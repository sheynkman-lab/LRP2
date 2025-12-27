# LRP2-Lite: Long-Read Proteogenomics Lite Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) [![run with conda](https://img.shields.io/badge/run%20with-conda-43b02a.svg?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**LRP2_Lite** is a bioinformatics pipeline for comprehensive long-read proteogenomics analysis of PacBio Iso-Seq data. It takes full-length non-chimeric (FLNC) BAM files as input, and performs isoform discovery, quality control, and protein prediction, outputting high-confidence proteomes and detailed classification reports.

The pipeline combines state-of-the-art tools for long-read RNA sequencing analysis with custom filtering and classification scripts following best practices to identify novel protein isoforms from PacBio data.

## Pipeline Summary

The LRP2_Lite pipeline consists of three major stages:

### 1. PacBio Iso-Seq Processing (`01_pacbio_isoseq`)
- Merge FLNC BAM files per sample (**pbtk pbmerge**)
- Cluster merged reads into consensus isoforms (**isoseq cluster**)
- Align consensus isoforms to reference genome (**pbmm2**)
- Collapse redundant isoforms based on alignment (**isoseq collapse**)

### 2. Transcript Quality Control and Filtering (`02_sqanti_transcript`)
- Perform comprehensive quality control and classification (**sqanti_qc**)
- Filter transcripts (**filter_sqanti**) by:
  - Protein-coding status
  - Internal priming artifacts
  - Template switching artifacts
  - Structural categories
- Generate CPM-normalized hash IDs for protein mapping (**filter_sqanti**)

### 3. Protein Prediction and Classification (`03_predicted_proteome`)
- Predict open reading frames (**cpat_orf**)
- Filter and select best ORF predictions (**filter_cpat**)
- Classify predicted proteins (**sqanti_protein_classification**)
- Generate high-confidence protein sets with custom UTR classification (**protein_utr_classification**)

### Key Features:
- **Full isoform resolution**: Leverages PacBio long reads for complete transcript characterization
- **Quality-based filtering**: Multi-stage artifact removal and quality control
- **Protein-level analysis**: ORF prediction and protein classification
- **Reference-based annotation**: Integration with GENCODE annotations
- **Species support**: Human and mouse genomes (via iGenomes and GENCODE)
- **Reproducible**: Fully containerized with Docker/Singularity support

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to test your setup with `-profile test` before running the workflow on actual data.

### Preparing Input Data

First, prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv:**

```csv
sample_name,bam,condition,replicate
control_chr22,/path/to/sample_data/230801_pacbio_rbfox2_control_chr22.flnc.bam,control,rep1
treatment_chr22,/path/to/sample_data/230801_pacbio_rbfox2_RB-G5_chr22.flnc.bam,treatment,rep1
```

Each row represents a PacBio Iso-Seq FLNC (Full-Length Non-Chimeric) BAM file. The required columns are:

- `sample_name`: Unique sample identifier (no spaces)
- `bam`: Absolute path to the FLNC BAM file
- `condition`: Sample condition or group (optional, used for downstream analysis)
- `replicate`: Replicate identifier (optional)

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters.

### Running the Pipeline

Now, you can run the pipeline using:

```bash
nextflow run /path/to/LRP2_lite \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --genome GRCh38 \
    -profile singularity
```

This will launch the pipeline with the default parameters using Singularity containers. See below for more details about available profiles and parameters. 

If you are running on an HPC environment, you may also run using SLURM, e.g.:

```bash
nextflow run /path/to/LRP2_lite \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --genome GRCh38 \
    -profile singularity,slurm
```

### Example Command with Additional Parameters

```bash
nextflow run /path/to/LRP2_lite \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --species human \
    --protein_coding_filter true \
    --cpat_coding_threshold 0.364 \
    -profile singularity \
    -resume
```

### Available Profiles

The pipeline supports multiple execution profiles:

- `docker`: Run using Docker containers (recommended)
- `singularity`: Run using Singularity containers
- `conda`: Run using Conda environments
- `test`: Run minimal test dataset
- `<institution>`: Institutional configs (if configured)

### Key Parameters

**Input/Output:**
- `--input`: Path to comma-separated samplesheet file (required)
- `--outdir`: Path to output directory (required)

**Reference Genome:**
- `--genome`: Reference genome version (default: `GRCh38`)
- `--fasta`: Path to reference FASTA (auto-set from iGenomes)
- `--gencode_gtf`: Path to GENCODE GTF annotation (auto-set from iGenomes)

**Filtering Options:**
- `--protein_coding_filter`: Keep only protein-coding genes (default: `true`)
- `--internal_priming_filter`: Remove internal priming artifacts (default: `true`)
- `--template_switching_filter`: Remove template switching artifacts (default: `true`)
- `--structure_filter`: Filtering level - `strict` or `all` (default: `strict`)

**Protein Prediction:**
- `--species`: Species for CPAT models - `human` or `mouse` (default: `human`)
- `--min_orf`: Minimum ORF length in nucleotides (default: `75`)
- `--top_orf`: Number of ORF candidates to report (default: `5`)
- `--cpat_coding_threshold`: Coding probability threshold (default: human=0.364, mouse=0.44)
- `--protein_class_keep`: Protein categories to retain (default: `FPM,IPM,NPC,NPE`)

**Iso-Seq Collapse Parameters:**
- `--max_fuzzy_junction`: Maximum junction position difference (default: `0`)
- `--max_5p_diff`: Max base-pair difference at 5' end (default: `100`)
- `--max_3p_diff`: Max base-pair difference at 3' end (default: `200`)

For a complete list of parameters, run:

```bash
nextflow run /path/to/LRP2_lite --help
```

## Sample Data

### Downloading Sample Data

Example sample data for testing the pipeline is available on Zenodo:

**Zenodo DOI:** [10.5281/zenodo.18065306](https://zenodo.org/records/18065306)

The sample dataset includes chromosome 22 data from a PacBio Iso-Seq experiment studying RBFOX2 splicing regulation:
- `230801_pacbio_rbfox2_control_chr22.flnc.bam` - Control sample (chr22 only)
- `230801_pacbio_rbfox2_RB-G5_chr22.flnc.bam` - RBFOX2 knockdown sample (chr22 only)

### Using Sample Data

An example samplesheet for the sample data is included in `sample_data/samplesheet.csv`. After downloading the BAM files from Zenodo, you can run the pipeline on the sample data:

```bash
# Download sample data from Zenodo
wget https://zenodo.org/records/18065306/files/230801_pacbio_rbfox2_control_chr22.flnc.bam
wget https://zenodo.org/records/18065306/files/230801_pacbio_rbfox2_RB-G5_chr22.flnc.bam

# Move files to sample_data directory
mv *.flnc.bam sample_data/

# Run pipeline with sample data
nextflow run /path/to/LRP2_lite \
    --input sample_data/samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    -profile singularity
```

> [!NOTE]
> The sample data is restricted to chromosome 22 to reduce file size and computation time, making it ideal for testing the pipeline and becoming familiar with the workflow. We highly recommend trying out the sample data prior to running with your own data.

## Pipeline Output

The pipeline generates comprehensive outputs organized in three main directories:

### Output Directory Structure

```
<outdir>/
├── 01_pacbio_isoseq/          # Iso-Seq processing results
│   ├── isoseq_merge/          # Merged FLNC BAM files
│   ├── isoseq_cluster/        # Clustered consensus reads
│   ├── isoseq_align/          # Aligned consensus isoforms
│   └── isoseq_collapse/       # Collapsed isoform GFF files
├── 02_sqanti_transcript/      # Transcript QC and filtering
│   ├── sqanti_qc/             # SQANTI3 classification reports
│   └── custom_filtered/       # Filtered transcript sets
├── 03_predicted_proteome/     # Protein predictions
│   ├── cpat_orf/              # ORF predictions
│   ├── filter_cpat/           # Filtered ORFs and CDS
│   ├── sqanti_protein/        # Protein classifications
│   └── protein_utr/           # High-confidence proteins
└── pipeline_info/             # Execution reports and logs
    ├── execution_report.html
    ├── execution_timeline.html
    └── lrp2_software_versions.yml
```

### Key Output Files

**Transcript-level:**
- Collapsed isoform annotations (GFF/GTF format)
- SQANTI3 classification reports
- Filtered transcript sequences (FASTA)
- Hash ID mappings with CPM values

**Protein-level:**
- Predicted protein sequences (FASTA)
- Protein classification tables
- CDS annotations (GTF)
- Protein-to-transcript mappings
- High-confidence proteome sets

**Quality Control:**
- SQANTI3 QC reports and plots
- FLNC read counts
- Gene discovery statistics
- Junction analysis

For detailed information about output files, please refer to the [output documentation](docs/output.md).

## Credits

### Development Team

The LRP2_Lite pipeline was developed by the Sheynkman Lab and Knowles Lab:

- **Megan Schertzer**, Sheynkman Lab - Pipeline and module code development
- **Julia Lewandowski**, Knowles Lab - Nextflow workflow implementation

We thank the following people for their extensive assistance in the development of this pipeline...

## License

This pipeline is released under the MIT License.

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use sheynkmanlab/lrp2 for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
