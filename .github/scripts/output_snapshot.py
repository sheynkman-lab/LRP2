#!/usr/bin/env python3
"""Snapshot LRP2 pipeline outputs and compare later runs against the snapshot.

A snapshot is a small JSON file listing every published output file (paths relative to the
pipeline's --outdir) with a SHA-256 of its normalised contents. We leverage snapshots in this 
context for our CI/CD pipeline to ensure output consistency across runs. For specific rules on
how to normalise and compare files, see the accompanying tests/output_snapshots/rules.txt file.

How to use: 

  # Record the expected outputs (e.g. on the HPC, from a known-good run):
  output_snapshot.py snapshot RESULTS_DIR -o tests/output_snapshots/NAME.json

  # Check a new run against it (exit code 1 if anything changed):
  output_snapshot.py compare tests/output_snapshots/NAME.json RESULTS_DIR

  # Compare two snapshots (e.g. the expected one and one saved by a CI run):
  output_snapshot.py compare expected.json actual.json

Normalisation for 'content' checks:
  * gzip files (including BGZF) are compared decompressed, so compression level and the
    timestamp inside the gzip header don't matter;
  * in text files, Nextflow work-dir paths (.../work/ab/cdef.../, .../work/stage-<uuid>/) become
    <WORK>, the results directory's own path becomes <OUTDIR>, and Windows line endings become
    Unix ones;
  * per-rule options can round decimal numbers, sort lines, or drop lines (see rules.txt).
"""

import argparse
import gzip
import hashlib
import json
import os
import re
import shlex
import sys
from datetime import datetime, timezone

SNAPSHOT_FORMAT = 1
CHUNK = 1024 * 1024
BINARY_SNIFF = 8192

WORK_PATH_RE = re.compile(
    r"(?:/[^/\s'\"=,;:|<>()\[\]{}]+)*/work/(?:[0-9a-f]{2}/[0-9a-f]{30}|stage-[0-9a-f-]{36})"
)
# Decimal numbers only (a '.' or an exponent), not part of a longer token such as a version
# string ("1.2.3"), an ID ("ENST00000456328.2") or a word.
FLOAT_RE = re.compile(r"(?<![\w.])[-+]?(?:\d+\.\d*|\.\d+|\d+(?=[eE]))(?:[eE][-+]?\d+)?(?![\w.])")


# ----------------------------------------------------------------------------- rules
class Rule:
    def __init__(self, pattern, check, rounding=None, sort=False, drop=None, lineno=0):
        if check not in ("ignore", "exists", "content"):
            raise ValueError(f"rules line {lineno}: unknown check '{check}' "
                             "(use ignore, exists or content)")
        if check != "content" and (rounding is not None or sort or drop):
            raise ValueError(f"rules line {lineno}: round/sort/drop only apply to 'content'")
        self.pattern = pattern
        self.check = check
        self.rounding = rounding
        self.sort = sort
        self.drop = drop or []
        self.drop_re = [re.compile(d) for d in self.drop]
        self.regex = glob_to_regex(pattern)

    def describe(self):
        """Canonical text stored in snapshots, so a rule change is noticed on compare."""
        parts = [self.check]
        if self.rounding is not None:
            parts.append(f"round={self.rounding}")
        if self.sort:
            parts.append("sort")
        parts += [f"drop={d}" for d in self.drop]
        return " ".join(parts)


def glob_to_regex(pattern):
    """'**/' = any number of folders (including none), '**' = anything, '*' = anything within
    one path component, '?' = one character within a component."""
    out, i = [], 0
    while i < len(pattern):
        if pattern.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        elif pattern[i] == "*":
            out.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            out.append("[^/]")
            i += 1
        else:
            out.append(re.escape(pattern[i]))
            i += 1
    return re.compile("^" + "".join(out) + "$")


def load_rules(path):
    rules = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            try:
                tokens = shlex.split(line, comments=True)
            except ValueError as err:
                raise ValueError(f"{path} line {lineno}: {err}") from None
            if not tokens:
                continue
            if len(tokens) < 2:
                raise ValueError(f"{path} line {lineno}: expected '<pattern> <check> [options]'")
            pattern, check, opts = tokens[0], tokens[1], tokens[2:]
            rounding, sort, drop = None, False, []
            for opt in opts:
                if opt.startswith("round="):
                    rounding = int(opt.split("=", 1)[1])
                    if not 1 <= rounding <= 17:
                        raise ValueError(f"{path} line {lineno}: round= must be 1-17")
                elif opt == "sort":
                    sort = True
                elif opt.startswith("drop="):
                    drop.append(opt.split("=", 1)[1])
                else:
                    raise ValueError(f"{path} line {lineno}: unknown option '{opt}'")
            rules.append(Rule(pattern, check, rounding, sort, drop, lineno))
    if not rules:
        raise ValueError(f"{path}: no rules found")
    return rules


def rule_for(relpath, rules):
    for rule in rules:
        if rule.regex.match(relpath):
            return rule
    return None


# ----------------------------------------------------------------------------- fingerprints
def open_maybe_gzip(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return open(path, "rb")


def round_floats(line, digits):
    def repl(m):
        value = float(m.group(0))
        text = f"{value:.{digits}g}"
        return "0" if text in ("-0", "+0") else text
    return FLOAT_RE.sub(repl, line)


def fingerprint(path, rule, replacements):
    size = os.path.getsize(path)
    if rule.check == "exists":
        return {"rule": rule.describe(), "bytes": size}

    digest = hashlib.sha256()
    with open_maybe_gzip(path) as fh:
        head = fh.read(BINARY_SNIFF)
        if b"\x00" in head:
            # Binary content: hash the (decompressed) bytes as they are
            digest.update(head)
            for chunk in iter(lambda: fh.read(CHUNK), b""):
                digest.update(chunk)
            return {"rule": rule.describe(), "sha256": digest.hexdigest(),
                    "kind": "binary", "bytes": size}

        def lines():
            rest = fh.read(CHUNK)
            buffer = head + rest
            pending = b""
            while buffer:
                buffer = pending + buffer
                parts = buffer.split(b"\n")
                pending = parts.pop()
                for raw in parts:
                    yield raw
                buffer = fh.read(CHUNK)
            if pending:
                yield pending

        def normalise(raw):
            text = raw.decode("utf-8", errors="replace").rstrip("\r")
            if rule.drop_re and any(r.search(text) for r in rule.drop_re):
                return None
            if "/work/" in text:
                text = WORK_PATH_RE.sub("<WORK>", text)
            for old, new in replacements:
                if old in text:
                    text = text.replace(old, new)
            if rule.rounding is not None:
                text = round_floats(text, rule.rounding)
            return text

        count = 0
        if rule.sort:
            kept = [t for t in (normalise(r) for r in lines()) if t is not None]
            kept.sort()
            for text in kept:
                digest.update(text.encode("utf-8") + b"\n")
            count = len(kept)
        else:
            for raw in lines():
                text = normalise(raw)
                if text is None:
                    continue
                digest.update(text.encode("utf-8") + b"\n")
                count += 1
    return {"rule": rule.describe(), "sha256": digest.hexdigest(),
            "kind": "text", "lines": count, "bytes": size}


def take_snapshot(results_dir, rules, extra_replacements=()):
    results_dir = os.path.abspath(results_dir)
    if not os.path.isdir(results_dir):
        raise SystemExit(f"ERROR: results directory not found: {results_dir}")
    # Longest first, so a path is replaced before any of its prefixes
    replacements = sorted({(results_dir, "<OUTDIR>"),
                           (os.path.realpath(results_dir), "<OUTDIR>"),
                           *extra_replacements}, key=lambda p: -len(p[0]))
    files = {}
    for root, dirs, names in os.walk(results_dir):
        dirs.sort()
        for name in sorted(names):
            full = os.path.join(root, name)
            rel = os.path.relpath(full, results_dir).replace(os.sep, "/")
            rule = rule_for(rel, rules)
            if rule is None:
                raise SystemExit(f"ERROR: no rule matches '{rel}'; end rules.txt with a '**' line")
            if rule.check == "ignore":
                continue
            if not os.path.isfile(full):  # broken symlink, socket, ...
                files[rel] = {"rule": rule.describe(), "error": "not a readable file"}
                continue
            files[rel] = fingerprint(full, rule, replacements)
    return {"format": SNAPSHOT_FORMAT,
            "source": os.path.basename(results_dir),
            "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "files": files}


# ----------------------------------------------------------------------------- compare
def load_snapshot(path):
    with open(path, encoding="utf-8") as fh:
        snap = json.load(fh)
    if snap.get("format") != SNAPSHOT_FORMAT or "files" not in snap:
        raise SystemExit(f"ERROR: {path} is not a snapshot made by this script (format {SNAPSHOT_FORMAT})")
    return snap


def compare(expected, actual, rules):
    exp, act = expected["files"], actual["files"]
    problems, ok = [], 0
    for rel in sorted(set(exp) | set(act)):
        e, a = exp.get(rel), act.get(rel)
        if e is None:
            problems.append(("new file", rel, "not in the expected snapshot"))
            continue
        if a is None:
            rule = rule_for(rel, rules) if rules else None
            if rule is not None and rule.check == "ignore":
                continue  # ignored since the snapshot was made
            problems.append(("missing", rel, "expected but not produced"))
            continue
        if e.get("rule") != a.get("rule"):
            problems.append(("rule changed", rel,
                             f"snapshot used '{e.get('rule')}', rules.txt now says "
                             f"'{a.get('rule')}': re-create the snapshot"))
            continue
        if "error" in a or "error" in e:
            problems.append(("unreadable", rel, a.get("error") or e.get("error")))
            continue
        if a["rule"] == "exists":
            if (e.get("bytes", 0) > 0) and a.get("bytes", 0) == 0:
                problems.append(("empty", rel, "expected a non-empty file"))
            else:
                ok += 1
            continue
        if e.get("sha256") != a.get("sha256"):
            detail = "contents differ"
            if e.get("kind") == "text" and a.get("kind") == "text":
                detail += f" ({e.get('lines')} lines expected, {a.get('lines')} now)"
            problems.append(("changed", rel, detail))
        else:
            ok += 1
    return problems, ok


def report(problems, ok, label, summary_path, github):
    kinds = {}
    for kind, _, _ in problems:
        kinds[kind] = kinds.get(kind, 0) + 1
    status = "match" if not problems else "DIFFER"
    print(f"{label}: outputs {status}. {ok} file(s) identical"
          + "".join(f", {n} {k}" for k, n in sorted(kinds.items())))
    for kind, rel, detail in problems:
        print(f"  {kind:12} {rel}  [{detail}]")
    if github:
        for kind, rel, detail in problems[:50]:
            print(f"::error title={label}: {kind}::{rel} ({detail})")
        if len(problems) > 50:
            print(f"::error title={label}::...and {len(problems) - 50} more differences (see job log)")
    if summary_path:
        with open(summary_path, "a", encoding="utf-8") as fh:
            icon = "✅" if not problems else "❌"
            fh.write(f"\n## {icon} Output comparison: {label}\n\n")
            fh.write(f"{ok} file(s) identical to the expected snapshot"
                     + "".join(f", **{n} {k}**" for k, n in sorted(kinds.items())) + ".\n")
            if problems:
                fh.write("\n| Difference | File | Detail |\n|---|---|---|\n")
                for kind, rel, detail in problems[:200]:
                    fh.write(f"| {kind} | `{rel}` | {detail} |\n")
                if len(problems) > 200:
                    fh.write(f"\n...and {len(problems) - 200} more (see the job log).\n")


# ----------------------------------------------------------------------------- CLI
def parse_replacements(values):
    pairs = []
    for value in values or []:
        if "=" not in value:
            raise SystemExit(f"ERROR: --replace expects PATH=TOKEN, got '{value}'")
        old, new = value.split("=", 1)
        pairs.append((old.rstrip("/"), new))
    return pairs


def main(argv=None):
    here = os.path.dirname(os.path.abspath(__file__))
    default_rules = os.path.normpath(os.path.join(here, "..", "..", "tests", "output_snapshots", "rules.txt"))

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--rules", default=default_rules, help=f"rules file (default: {default_rules})")
    common.add_argument("--replace", action="append", metavar="PATH=TOKEN",
                        help="also replace PATH with TOKEN in text files (repeatable), e.g. the "
                             "pipeline's launch directory, which differs between machines")

    p_snap = sub.add_parser("snapshot", parents=[common], help="record a results directory")
    p_snap.add_argument("results_dir")
    p_snap.add_argument("-o", "--output", required=True, help="snapshot JSON file to write")

    p_cmp = sub.add_parser("compare", parents=[common], help="compare against a snapshot")
    p_cmp.add_argument("expected", help="expected snapshot JSON")
    p_cmp.add_argument("actual", help="results directory, or a snapshot JSON")
    p_cmp.add_argument("--label", default="outputs", help="name shown in messages")
    p_cmp.add_argument("--write-actual", metavar="FILE",
                       help="also save the actual snapshot here (to replace the expected one "
                            "if the change is intended)")
    p_cmp.add_argument("--summary", metavar="FILE", help="append a Markdown report (GitHub step summary)")
    p_cmp.add_argument("--github", action="store_true", help="print GitHub ::error annotations")

    args = parser.parse_args(argv)
    try:
        rules = load_rules(args.rules)
    except (OSError, ValueError) as err:
        raise SystemExit(f"ERROR: {err}")
    replacements = parse_replacements(args.replace)

    if args.command == "snapshot":
        snap = take_snapshot(args.results_dir, rules, replacements)
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as fh:
            json.dump(snap, fh, indent=1, sort_keys=True)
            fh.write("\n")
        n = len(snap["files"])
        print(f"Wrote {args.output}: {n} file(s) from {args.results_dir}")
        return 0

    if not os.path.isfile(args.expected):
        raise SystemExit(f"ERROR: expected snapshot not found: {args.expected}")
    expected = load_snapshot(args.expected)
    if os.path.isfile(args.actual) and args.actual.endswith(".json"):
        actual = load_snapshot(args.actual)
    else:
        actual = take_snapshot(args.actual, rules, replacements)
    if args.write_actual:
        os.makedirs(os.path.dirname(os.path.abspath(args.write_actual)), exist_ok=True)
        with open(args.write_actual, "w", encoding="utf-8") as fh:
            json.dump(actual, fh, indent=1, sort_keys=True)
            fh.write("\n")
    problems, ok = compare(expected, actual, rules)
    report(problems, ok, args.label, args.summary, args.github)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
