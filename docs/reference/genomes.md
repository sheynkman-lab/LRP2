# Reference Genomes

The pipeline supports human and mouse data using GENCODE reference genomes across multiple versions.

## Supported Genomes

### Human (GRCh38)

- `GRCh38.p14.v49`
- `GRCh38.p14.v48`
- `GRCh38.p14.v47`
- `GRCh38.p14.v46`
- `GRCh38.p14.v45`
- `GRCh38.p14.v44`
- `GRCh38.p13.v43`
- `GRCh38.p13.v42`
- `GRCh38.p13.v41`
- `GRCh38.p13.v40`

### Human (GRCh37)

- `GRCh37.p13.v19` (and other GRCh37 versions)

### Mouse (GRCm39)

- `GRCm39.vM38`
- `GRCm39.vM37`
- `GRCm39.vM36`
- `GRCm39.vM35`
- `GRCm39.vM34`

## How It Works

The pipeline automatically downloads the appropriate FASTA and GTF files based on your `--genome` selection. Species is auto-detected from `--genome` and determines which CPAT model (human or mouse) is used for ORF prediction.

## Configuration

See `conf/gencode_references.config` for the full list of supported versions and their download URLs.

## Usage

Specify the genome version using the `--genome` parameter:

```bash
--genome GRCh38.p14.v49
```

## Future Support

Support for RefSeq / igenomes and custom references is under active development.

!!! note "Custom Genomes"
    If you need to use a custom reference genome not listed here, please contact the development team or check for updates in future releases.
