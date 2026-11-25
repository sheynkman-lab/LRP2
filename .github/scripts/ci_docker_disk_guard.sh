#!/usr/bin/env bash
# Keeps Docker's image store from filling the runner's disk during a Nextflow run.
#
# Docker keeps every image it has pulled until the job ends. One full pipeline run pulls about ten
# images (two SQANTI3 versions, ProteoWizard, FragPipe, PyTorch, ...) that together need more disk
# than the runner has left, so the last pull fails with "failed to register layer: ... no space
# left on device". Most of those images are only needed by one early step.
#
# This loop watches the disk Docker stores images on. While free space is below MIN_FREE_GB it
# removes images that no container is using. Docker refuses to remove an image a container still
# uses, so running tasks are never affected, and an image that is needed again later is simply
# pulled again. Above the threshold it does nothing.
#
# Usage: ci_docker_disk_guard.sh [MIN_FREE_GB] [INTERVAL_SECONDS] [GRACE_SECONDS]   (run in background)
set -u

MIN_FREE_GB="${1:-20}"     # start removing idle images below this much free space
INTERVAL="${2:-15}"        # seconds between checks
GRACE="${3:-300}"          # leave images pulled less than this many seconds ago alone: 'docker run'
                           # pulls an image and creates its container a moment later
REPORT_EVERY=600           # seconds between routine "free space" lines in the log

log() { echo "$(date -u '+%H:%M:%S') $*"; }
free_gb() { df -BG --output=avail "$1" 2> /dev/null | tail -n 1 | tr -dc '0-9'; }

log "watching Docker's disk: remove idle images when free space < ${MIN_FREE_GB} GB (check every ${INTERVAL}s)"
last_report=0
last_idle_note=0
while sleep "$INTERVAL"; do
  root=$(docker info --format '{{.DockerRootDir}}' 2> /dev/null) || continue
  free=$(free_gb "$root")
  [[ -n "$free" ]] || continue
  now=$(date +%s)
  if (( now - last_report >= REPORT_EVERY )); then
    log "free on ${root}: ${free} GB; images stored: $(docker image ls --quiet 2> /dev/null | sort -u | wc -l)"
    last_report=$now
  fi
  (( free < MIN_FREE_GB )) || continue

  removed=0
  while IFS= read -r ref; do
    [[ -z "$ref" || "$ref" == *"<none>"* ]] && continue
    pulled=$(docker image inspect --format '{{json .Metadata.LastTagTime}}' "$ref" 2> /dev/null | tr -d '"')
    pulled_s=$(date -d "$pulled" +%s 2> /dev/null || echo 0)
    (( now - pulled_s < GRACE )) && continue
    # Fails (and is ignored) when a container is using the image
    if docker image rm "$ref" > /dev/null 2>&1; then
      log "removed idle image ${ref}"
      removed=$((removed + 1))
    fi
  done < <(docker image ls --format '{{.Repository}}:{{.Tag}}' 2> /dev/null | sort -u)
  # Untagged leftovers (e.g. layers of an image whose tag moved)
  docker image prune --force > /dev/null 2>&1 || true

  if (( removed > 0 )); then
    log "free space was ${free} GB; removed ${removed} idle image(s); now $(free_gb "$root") GB free"
  elif (( now - last_idle_note >= 300 )); then
    log "free space is ${free} GB but every stored image is in use or was just pulled; nothing to remove"
    last_idle_note=$now
  fi
done
