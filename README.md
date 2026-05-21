![LRP2 Logo](assets/LRP2-Workflow.lightmode.drawio.png#gh-light-mode-only)
![LRP2 Logo](assets/LRP2-Workflow.darkmode.drawio.png#gh-dark-mode-only)

# LRP2: Long-Read Proteogenomics Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) [![run with conda](https://img.shields.io/badge/run%20with-conda-43b02a.svg?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**LRP2** is a Nextflow pipeline for long-read proteogenomics analysis. It takes PacBio full-length non-chimeric (FLNC) reads and/or raw mass spectrometry files as input, performs transcript discovery and quality control, ORF prediction and quality control, differential analysis, and proteomics database generation and search to validate protein isoforms.

## Pipeline Overview

The LRP2 Nextflow pipeline consists of five subworkflows:

| Subworkflow | Description |
|-------------|-------------|
| 1. PacBio Isocall | Align FLNC reads and collapse to isoforms with PacBio Isocall |
| 2. Transcriptome | Classify transcripts with SQANTI3, filter artifacts, assign deterministic hash-based isoform IDs |
| 3. Predicted proteome | Predict ORFs with CPAT, classify proteins with SQANTI protein |
| 4. Multi-sample analysis *(optional)* | Differential splicing (LeafCutter), expression and usage (edgeR, DRIMSeq) |
| 5. Proteomics *(optional)* | Build custom reference database, convert raw MS files, search with FragPipe or MetaMorpheus, map peptides to isoforms |

> **Note**: Differential usage analyses (DTU/DPU with DRIMSeq) require biological replicates. Datasets without replicates may lead to unexpected results with the pipeline currently.

> **Note**: The PROTEOMICS subworkflow can only run when protein samples (sample_type='protein') are provided in the samplesheet. When both RNA and protein samples are present, it searches against a concatenated database of the predicted proteome from Stage 3 plus the reference protein FASTA. For protein-only samples, it uses only the reference protein database via `--protein_fasta` (or auto-detected from the GENCODE genome if specified). Select the search engine with `--protein_search` (options: `fragpipe` or `metamorpheus`, default: `fragpipe`). 

## Quick Start

### Prerequisites

- **Nextflow** ≥ 24.04.2 ([install guide](https://www.nextflow.io/docs/latest/install.html))
- **Singularity/Apptainer** or **Docker** for containerized dependencies

### System Requirements and HPC Recommendations

**LRP2 is primarily designed for High-Performance Computing (HPC) environments.** Although the test dataset can run on a local machine (minimum of 4 CPUs and 32GB RAM recommended), real-world datasets require substantial computational resources to be processed efficiently.

#### A Primer on Containers

LRP2 uses containers to manage software dependencies. Containers allow for the packaging of libraries, code, and configurations such that each tool in the pipeline can be reliably run in any computing environment without compatibility issues. 

To run LRP2, you **must have one** of the following installed:

1. **Singularity/Apptainer** (required for HPC): Most HPC systems have Singularity or Apptainer pre-installed as a module. You can check this by running `module avail singularity` or `module avail apptainer`. If not available, contact your HPC administrator or see [Apptainer installation guide](https://apptainer.org/docs/admin/main/installation.html), 

2. **Docker** (may be used for local systems):
    - Installation guides: [Docker Desktop](https://docs.docker.com/get-docker/) (Mac/Windows) or [Docker Engine](https://docs.docker.com/engine/install/) (Linux)
> **Note**: Docker typically requires root/admin privileges, and may not be available on shared HPC systems. We therefore strongly recommend the use of Singularity/Apptainer. 

The pipeline will automatically pull and cache container images on first run. Singularity images are cached in `work/singularity/` by default.

### Clone the repository
```bash
git clone https://github.com/sheynkman-lab/LRP2.git
cd LRP2
```

### On an HPC (e.g., UVA Rivanna)

Start a persistent terminal session so the pipeline keeps running if you lose your connection:
```bash
screen -S lrp2
```
> **Tip**: To detach from screen, press `Ctrl+A` then `D`. To reattach later: `screen -r lrp2`

> **Tip**: Certain HPC systems (e.g. UVA Rivanna) only support `screen`, but you can use terminal multiplexer by running `tmux new -s lrp2` on other systems if supported/preferred.

Request an interactive job with enough resources for the test dataset:
```bash
ijob -c 4 --mem=64G -p your_slurm_partition -A your_allocation --time=4:00:00
```
> **Note**: Adjust for your HPC system. Replace `your_slurm_partition` with your SLURM partition and `your_allocation` with your SLURM allocation group. The `-c` (CPUs), `--mem` (memory), and `--time` values above are sufficient for the test dataset, but should be increased for larger datasets.

Load the required modules:
```bash
module load nextflow apptainer
```

### Run the RNA-only test dataset

From the `LRP2` directory, you can run the RNA-only test two ways:

To run locally:
```bash
nextflow run . -profile test_rna,singularity --outdir test_results
```

To submit individual tasks to the SLURM scheduler (recommended):
```bash
nextflow run . -profile test_rna,singularity,slurm --outdir test_results
```

### Run the RNA + DDA proteomics test dataset

1. To use FragPipe, you must be an academic user and accept the license agreement for MSFragger, IonQuant, and diaTracer. Before using FragPipe in LRP2 for the first time, please review the license agreement here: https://msfragger.arsci.com/upgrader/LICENSE-ACADEMIC.pdf. Once you are ready to accept the license agreement, you may run the following curl command in your terminal with your information substituted for FIRST_NAME, LAST_NAME, EMAIL, and ACADEMIC_INSTITUTION_NAME:

```bash
curl --location --request POST \
    'https://msfragger-upgrader.nesvilab.org/upgrader/upgrade_download.php' \
    --form 'transfer="academic"' \
    --form 'agreement2="true"' \
    --form 'agreement3="true"' \
    --form "first_name=YOUR_FIRST_NAME" \
    --form "last_name=YOUR_LAST_NAME" \
    --form "email=YOUR_EMAIL" \
    --form "organization=YOUR_INSTITUTION" \
    --form "download=4.4.1\$zip" \
    --form 'is_fragpipe="true"' \
    > /dev/null 2>&1
```
> **Note**: Tokens expire quickly. You will need a new token for each run.

2. Check your email for a 6-digit token.

3. From the `LRP2` directory, you can run the RNA + DDA proteomics test with your token two ways:

To run locally:
```bash
nextflow run . -profile test_dda,singularity --outdir test_results_dda \
    --fragpipe_token "YOUR_TOKEN"
```

To submit individual tasks to the SLURM scheduler (recommended):
```bash
nextflow run . -profile test_dda,singularity,slurm --outdir test_results_dda \
    --fragpipe_token "YOUR_TOKEN"
```

> **Note**: The `test_dda` profile automatically sets `--protein_search fragpipe` and `--fragpipe_license_accept true`.

## Preparing Input Data

### Input File Requirements

**RNA samples** must be provided as **PacBio full-length non-chimeric (FLNC) reads** as outputted by PacBio Isoseq refine, in either BAM or FASTQ format. It is assumed that input files are **post-processed** and have already undergone deconcatenation, demultiplexing, and primer removal. **Do NOT** provide raw subreads or CCS reads directly from the sequencer

**Protein samples** may be either **DDA** or **DIA**, and can be provided in `.mzML` format or any other non-mzML format (e.g. `.raw`).

### Samplesheet Structure

Prepare a samplesheet CSV describing your input data:

![Samplesheet structure](assets/samplesheet_structure.png)

> **Note**: Each RNA sample must have a unique `sample_name`, as these are used by Isocall to label count matrix columns. For protein samples, all raw files originating from the same biological sample (e.g., multiple fractions or injection replicates) should share the same `sample_name` so they are combined and searched together in FragPipe. This `sample_name` should match the corresponding RNA sample to link the RNA and protein data. The predicted proteome from that sample will be included in the proteomics search database. If a protein sample has no matched RNA sample, you can assign it a unique `sample_name` that does not match any RNA sample. In this case, only the GENCODE reference proteome will be used as the proteomics search database.

**Required columns:**

- `sample_name`: Each RNA sample must have a distinct value. Do not include any spaces in this value.
- `sample_path`: Absolute path to the input file. 
  - RNA samples should be PacBio FLNC `.bam` or `.fastq` files
  - Protein samples should be `.raw` or `.mzML` files
- `condition`: Sample group (e.g., "control", "treatment"). Used for differential analysis, which performs pairwise comparisons between groups. Two or more groups are supported. If you do not want differential analysis, assign the same condition to all samples. Do not include any spaces in this value.
- `sample_type`: Must be either `RNA` or `protein`.
- `mass_spec_type`: Must be either `DDA` or `DIA`. Required for protein samples. For RNA samples, specify `none` for this column. 


## Running the Pipeline

### Basic run with default parameters (RNA only)

You can run the pipeline from any directory by providing the path to the LRP2 repository.

Using tmux and an interactive node without SLURM (tasks run locally on the node — suitable for small tests but not large datasets): 
```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v49 \
    -profile singularity
```

On an HPC with SLURM (recommended — tasks are submitted as individual jobs for parallel execution):
```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v49 \
    -profile singularity,slurm
```
> **Note**: Running on a non-SLURM scheduler (e.g., LSF, PBS)? Contact us for support: cwp5au@virginia.edu
> **Note**: Differential analysis will run automatically when two or more conditions are specified in the samplesheet.

### With proteomics (FragPipe)

Requires protein samples in the samplesheet and a FragPipe academic license token (see [Run the RNA + DDA proteomics test dataset](#run-the-rna--dda-proteomics-test-dataset) above for obtaining a token).
```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v49 \
    --protein_search fragpipe \
    --fragpipe_token YOUR_TOKEN \
    -profile singularity,slurm
```
> **Note**: When both RNA and protein samples are provided, the pipeline searches against the predicted proteome combined with a reference protein FASTA (auto-detected from the genome or provided via `--protein_fasta`). For protein-only samples, only the reference database is used.

## Profile Options

Nextflow profiles control how the pipeline executes. Multiple profiles can be combined using commas (e.g., `-profile test_rna,singularity,slurm`).

### Container Profiles (choose ONE)

| Profile | Description | Best For |
|---------|-------------|----------|
| `singularity` | Use Singularity/Apptainer containers | **HPC systems** (most common) |
| `docker` | Use Docker containers | Local machines with Docker installed |
| `conda` | Use Conda environments | Systems without container support (slower) |

> **Note**: We have extensively tested the pipeline with `singularity` locally and on HPC systems, and recommend its usage. You may use `docker` on local machines. We do not recommend the use of `conda` except as a last resort due to it lacking the same reproducibility as containers. 

### Executor Profiles (optional)

| Profile | Description | When to Use |
|---------|-------------|-------------|
| `slurm` | Submit jobs to SLURM scheduler on an HPC | HPC environment |
| `lsf` | Submit jobs to LSF scheduler on an HPC | HPC environment |

If you do not specify an executor, LRP2 will run locally on the current node, which is suitable for an interactive session in an HPC environment or on a local machine.

> **Important**: When using `slurm` or `lsf`, Nextflow submits individual pipeline tasks as separate jobs. Without a scheduler profile, all tasks run on the node where you launch Nextflow (requires sufficient resources). If you intend to run locally, you may need to lower resource requirements in `conf/base.config`. 

### Test Profiles

| Profile | Description | Dataset |
|---------|-------------|---------|
| `test_rna` | RNA-only test dataset | Runs RNA subworkflows (S1 - S4) |
| `test_dda` | RNA + DDA proteomics test | Runs all subworkflows (S1 - S5) with FragPipe DDA search |
| `test_dia` | RNA + DIA proteomics test | Runs all subworkflows (S1 - S5) with FragPipe DIA search |

### Example Profile Combinations

**HPC with SLURM (recommended for production):**
```bash
-profile singularity,slurm
```

**Quick RNA test on HPC:**
```bash
-profile test_rna,singularity,slurm
```

**Local machine with Docker:**
```bash
-profile docker
```

## Reference Genome Support

The pipeline supports multiple multi-species (human, mouse) reference genome sources: 

1. **GENCODE genomes** (recommended) — High-quality annotations with multiple release versions
   - Human: `GRCh38.p14.v49`, `GRCh38.p14.v48`, ..., `GRCh38.p14.v44`, `GRCh38.p13.v43`, ..., `GRCh37.p13.v19`
   - Mouse: `GRCm39.vM38`, `GRCm39.vM37`, `GRCm39.vM36`, `GRCm39.vM35`, `GRCm39.vM34`
2. **Custom references** — Provide your own FASTA and GTF files by setting the `--fasta` and `--gtf` parameters
3. **RefSeq genomes (under development)** — Standard genome builds using NCBI/Ensembl annotations (e.g., `GRCh38`, `GRCm38`)

Currently, it is recommended to use GENCODE reference, with support for any version. The pipeline automatically downloads the appropriate FASTA and GTF files based on your `--genome` selection.

Support for RefSeq / igenomes is under active development.

### Viewing Available Genomes

To see all available preconfigured `--genome` options, run this command from the LRP2 directory:

```bash
grep -oE "^        '[^']+'" conf/gencode_references.config | tr -d "'" | tr -d ' ' | sort
```

### Multi-Species Support

The pipeline supports both human and mouse data analysis. The species being analyzed is automatically detected from the `--genome` parameter when using predefined GENCODE references:

**Human data:**
```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v49 \
    -profile singularity,slurm
```

**Mouse data:**
```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCm39.vM34 \
    -profile singularity,slurm
```

Species determines which CPAT model files (human vs mouse-specific models) are used for ORF prediction.

### Using Custom References

When providing custom FASTA and GTF files instead of using a predefined `--genome`, you **must explicitly specify** the `--species` parameter:

```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --fasta /path/to/custom_genome.fa \
    --gtf /path/to/custom_annotation.gtf \
    --species mouse \
    -profile singularity,slurm
```

> **Important**: The `--species` parameter is required when using custom references to ensure the correct CPAT models are used for ORF prediction. Valid values are `human` or `mouse`.

## Parameters

### Input/Output

- `--input` — Path to samplesheet CSV (required)
- `--outdir` — Path to output directory (required)
- `--genome` — Reference genome, specify version (default: `GRCh38.p14.v49`)

### PacBio Isocall Subworkflow

- `--min_read_support` — Minimum read support for transcripts (default: `3`)
- `--isocall_config`: Path to custom Isocall configuration TOML file (default: `bin/isocall_config.toml`)
  
### Transcriptome Subworkflow

- `--protein_coding_filter` — Keep only protein-coding genes (default: `true`)
- `--internal_priming_filter` — Remove internal priming artifacts (default: `true`)
- `--template_switching_filter` — Remove template switching artifacts (default: `true`)
- `--transcript_class_keep` — Structural categories to retain: FSM, ISM, NIC, NNC, or ALL (default: `FSM,ISM,NIC,NNC`)

### Predicted Proteome Subworkflow

- `--species` — Species: `human` or `mouse`. Auto-detected from `--genome` when using predefined GENCODE references. **Required** when using custom FASTA/GTF files.
- `--min_orf` — Minimum ORF length in nucleotides (default: `75`)
- `--cpat_coding_threshold` — Coding probability threshold (default: human=0.364, mouse=0.44)
- `--protein_class_keep` — Protein categories to retain (default: `FPM,NPC,NPE`)

### Differential Analysis

- `--min_samples_per_intron`: Minimum samples per intron for leafcutter (default: `2`)
- `--min_samples_per_group`: Minimum samples per group for leafcutter (default: `1`)
- `--min_usage_ratio`: Minimum junction usage ratio for filtering (default: `0.01`)
  
### Proteomics Subworkflow

- `--protein_search` — Search engine: `fragpipe` or `metamorpheus` (default: `fragpipe`)

For a complete list of parameters:
```bash
nextflow run /path/to/LRP2 --help
```

> [!WARNING]
> Please provide pipeline parameters via the CLI as shown or using the Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters.

## Pipeline Output

Each subworkflow outputs to numbered module directories. The final module in each subworkflow typically contains the key results, while earlier modules contain intermediate files.

```
<outdir>/
├── S1_PACBIO_ISOCALL/                   
│   ├── M1_ISOCALL_ALIGN/               
│   ├── M2_ISOCALL_PROFILE/             
│   ├── M3_ISOCALL_PREP/                 
│   ├── M4_ISOCALL_MERGE/                
│   └── M5_ISOCALL_CALL/                 # GTF of transcript structures and count matrix
├── S2_TRANSCRIPTOME/                   
│   ├── M1_SQANTI_QC/
│   ├── M2_GENERATE_HASHIDS/                
│   └── M2_FILTER_TRANSCRIPTOME/         # GTF, BED12, DNA FASTA, count matrix of the refined transcriptome (technical artifacts removed)
├── S3_PREDICTED_PROTEOME/               
│   ├── M1_CPAT_ORF/                     
│   ├── M2_FILTER_CPAT/                  # GTF with exon and CDS type columns for single best ORF identified per transcript 
│   ├── M3_SQANTI_PROTEIN/
│   └── M4_PROTEIN_CLASSIFICATION/       # GTF, BED12, protein FASTA, count matrix collpased to distinct ORFs 
├── S4_MULTISAMPLE_ANALYSIS/             # (optional)
│   ├── M1_LEAFCUTTER_LONGREAD/          # Differential splicing results
│   └── M2_DIFFERENTIAL_EXPRESSION/      # Differential expression/usage
│       ├── differential_gene_expression/
│       │   ├── *_DGE_edgeR_results.txt
│       │   ├── *_DGE_edgeR_raw_CPM_matrix.txt
│       │   ├── *_DGE_edgeR_normalized_CPM_matrix.txt
│       │   └── *_DGE_MD_plot.pdf
│       ├── differential_transcript_expression/
│       │   ├── *_DTE_edgeR_results.txt
│       │   └── *_DTE_MD_plot.pdf
│       ├── differential_ORF_expression/
│       │   ├── *_DE_ORF_edgeR_results.txt
│       │   └── *_DE_ORF_MD_plot.pdf
│       ├── differential_transcript_usage/
│       │   └── *_DTU_transcript_DRIMSeq_summary.txt
│       └── differential_ORF_usage/
│           └── *_DU_ORF_DRIMSeq_summary.txt
├── S5_PROTEOMICS/                       # (optional)
│   ├── M1_BUILD_PROTEOME_REFERENCE
│   ├── M2_MSCONVERT_MZML/         
│   ├── M3_FRAGPIPE/                     
│   └── M4_NOVEL_PEPTIDES/               # BED12 of peptides mapped to genome, summary table of novel and annotated peptides mapped to isoforms
└── pipeline_info/                       # Execution reports and logs
    ├── execution_report.html
    ├── execution_timeline.html
    └── lrp2_software_versions.yml
```

For detailed information about output files, please refer to the [output documentation](docs/output.md).

## Credits

### Development Team

The LRP2 pipeline was developed through a collaboration by the Sheynkman Lab and Knowles Lab:

- **Megan D. Schertzer**, Sheynkman Lab - Module code development
- **Julia T. Lewandowski**, Knowles Lab - Pipeline implementation

We thank the following people for their extensive assistance in the development of this pipeline: 
- **Emily F. Watts**, Sheynkman Lab - Contributions to LRP and conception of multi-sample analysis subworkflow.
- **Madison M. Mehlferber**, Sheynkman Lab - Contributor to the original LRP pipeline. Continued pipeline testing and feedback.
- **Will Rosenow**, Sheynkman Lab - Pipeline testing and feedback.
- **Scott I. Adamson**, Knowles Lab - Development of leafcutter-py. 
- **Jocelyne Bruand**, Pacific Biosciences - Development of Isocall. 
- **Elizabeth Tseng**, Pacific Biosciences - Development of Isocall. 
- **Egor Dolzhenko**, Pacific Bioscience - Lead Developer of Isocall. 

We especially thank the PIs associated with this project: 
- **David A. Knowles**, Development of Leafcutter and project support / funding 
- **Gloria Sheynkman**, Development/conceptualization of LRP and project support / funding

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
