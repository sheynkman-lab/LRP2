![LRP2 Logo](assets/LRP2-Workflow.lightmode.drawio.png#gh-light-mode-only)
![LRP2 Logo](assets/LRP2-Workflow.darkmode.drawio.png#gh-dark-mode-only)

# LRP2: Long-Read Proteogenomics Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) [![run with conda](https://img.shields.io/badge/run%20with-conda-43b02a.svg?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**LRP2** is a scalable, end-to-end long-read proteogenomics pipeline built in Nextflow. It identifies and validates protein isoforms by integrating PacBio long-read RNA-seq with mass spectrometry. Starting from full-length non-chimeric (FLNC) reads and/or raw MS files, LRP2 performs transcript discovery, ORF prediction, differential analysis, and mass spectrometry-based protein identification. 

## Pipeline Overview

The LRP2 Nextflow pipeline consists of five subworkflows:

| Subworkflow | Description |
|-------------|-------------|
| 1. PacBio Isocall | Align FLNC reads and collapse to isoforms with PacBio Isocall |
| 2. Transcriptome | Classify transcripts with SQANTI3, filter artifacts, assign deterministic hash-based isoform IDs |
| 3. Predicted proteome | Predict ORFs with CPAT, classify proteins with SQANTI protein |
| 4. Multi-sample analysis | Differential expression and usage (edgeR, DRIMSeq), differential splicing (LR LeafCutter, preliminary implementation) |
| 5. Proteomics | Build custom reference database, convert raw MS files, search with FragPipe or MetaMorpheus, map peptides to isoforms |
 
## Quick Start

> **Note**: The Quick Start uses SLURM on a UVA Rivanna-style HPC setup, which is what the pipeline has been most extensively tested with. LSF is also supported via the `lsf` profile. For other schedulers, see [Support and customization](#support-and-customization).
>
> Specifics vary by cluster — account flags (`-A`), partition names (`-p`), and module names (`nextflow`, `apptainer`) may differ on your system. Check with your HPC documentation or admin.

### Prerequisites

- **Nextflow** ≥ 24.04.2 ([install guide](https://www.nextflow.io/docs/latest/install.html))
- **Singularity/Apptainer** or **Docker** for containerized dependencies

### System Requirements and HPC Recommendations

**LRP2 is primarily designed for High-Performance Computing (HPC) environments.** Although the test dataset can run on a local machine (minimum of 4 CPUs and 32GB RAM recommended), real-world datasets require substantial computational resources to be processed efficiently.

### Containers

LRP2 uses containers to manage software dependencies. Containers allow for the packaging of libraries, code, and configurations such that each tool in the pipeline can be reliably run in any computing environment without compatibility issues. 

To run LRP2, you **must have one** of the following installed:

1. **Singularity/Apptainer** (required for HPC): Most HPC systems have Singularity or Apptainer pre-installed as a module. You can check this by running `module avail singularity` or `module avail apptainer`. If not available, contact your HPC administrator or see [Apptainer installation guide](https://apptainer.org/docs/admin/main/installation.html). 

2. **Docker** (may be used for local systems):
    - Installation guides: [Docker Desktop](https://docs.docker.com/get-docker/) (Mac/Windows) or [Docker Engine](https://docs.docker.com/engine/install/) (Linux)
      
> **Note**: Docker requires root/admin privileges, and is typically not permitted on shared HPC systems. We therefore strongly recommend the use of Singularity/Apptainer.

The pipeline will automatically pull and cache container images on first run. Singularity images are cached in `work/singularity/` by default.

### Clone the repository
```bash
git clone https://github.com/sheynkman-lab/LRP2.git
cd LRP2
```

### On an HPC (SLURM example)

Start a persistent terminal session so the pipeline keeps running if you lose your connection:
```bash
screen -S lrp2
```
> **Tip**: To detach from screen, press `Ctrl+A` then `D`. To reattach later: `screen -r lrp2`.
>
> Certain HPC systems (e.g. UVA Rivanna) only support `screen`. On systems that support it, you can use the tmux terminal multiplexer instead by running `tmux new -s lrp2`.

Request an interactive job with enough resources for the test dataset:
```bash
salloc -c 4 --mem=64G -p your_slurm_partition -A your_allocation --time=4:00:00
```

> **Note**: Adjust for your HPC system. Replace `your_slurm_partition` with your SLURM partition and `your_allocation` with your SLURM allocation group. UVA Rivanna users can substitute `ijob` for `salloc`. The `-c` (CPUs), `--mem` (memory), and `--time` values above are sufficient for the test dataset, but should be increased for larger datasets.

Load the required modules:
```bash
module load nextflow apptainer
```

### Run the RNA-only test dataset

From the `LRP2` directory:

```bash
nextflow run . \
    -profile test_rna,singularity,slurm \
    --outdir test_rna_results \
    --hpc_partition your_partition \
    --hpc_cluster_options "--account=your_allocation"
```
> **Note**: Replace `your_partition` with your cluster's partition name, and adjust `--hpc_cluster_options` for your cluster (e.g., `--account=` for SLURM, `-P` for LSF). LSF users: also swap `slurm` for `lsf` in the profile. See [HPC Scheduler Options](#hpc-scheduler-options).

> **Note**: To run locally on your current node instead of submitting to SLURM, drop the `slurm` profile: `-profile test_rna,singularity`.

### Run the RNA + DDA proteomics test dataset

**Step 1: Get a FragPipe academic license token.**

FragPipe requires an academic license for MSFragger, IonQuant, and diaTracer. Before using FragPipe in LRP2 for the first time, review the [academic license agreement](https://msfragger.arsci.com/upgrader/LICENSE-ACADEMIC.pdf). To accept the license and request a one-time token, run the following curl command in your terminal, substituting your information for `YOUR_FIRST_NAME`, `YOUR_LAST_NAME`, `YOUR_EMAIL`, and `YOUR_INSTITUTION`:

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

> **Note**: Non-academic users: see [Support and customization](#support-and-customization).

**Step 2: Check your email for a 6-digit token.**

**Step 3: Run the test.**

From the `LRP2` directory:

```bash
nextflow run . \
    -profile test_dda,singularity,slurm \
    --outdir test_results_dda \
    --fragpipe_token "YOUR_TOKEN" \
    --hpc_partition your_partition \
    --hpc_cluster_options "--account=your_allocation"
```
> **Note**: The `test_dda` profile automatically sets `--protein_search fragpipe` and `--fragpipe_license_accept true`.

> **Note**: To run locally instead of submitting to SLURM, drop the `slurm` profile: `-profile test_dda,singularity`.

## Preparing Input Data

### Input File Requirements

**RNA samples** must be provided as **PacBio full-length non-chimeric (FLNC) reads** as outputted by PacBio Isoseq refine, in either BAM or FASTQ format. It is assumed that input files are **post-processed** and have already undergone deconcatenation, demultiplexing, and primer removal. **Do NOT** provide raw subreads or CCS reads directly from the sequencer.

**Protein samples** may be either **DDA** or **DIA**, and can be provided in `.mzML` format or vendor-specific raw formats.

### Samplesheet Structure

Prepare a comma-delimited samplesheet (.csv) describing your input data:

![Samplesheet structure](assets/samplesheet_structure.png)

> **Note**: Sample naming conventions in the samplesheet:
>
> - **RNA samples**: Each RNA sample must have a unique `sample_name`. These are used by Isocall to label count matrix columns.
> - **Protein samples**: All raw files from the same biological sample (e.g., multiple fractions or injection replicates) must share the same `sample_name` so they are combined and searched together in FragPipe.
> - **Matched RNA + protein samples**: Use the same `sample_name` for the RNA and protein entries. The predicted proteome from that sample will be included in the proteomics search database.
> - **Unmatched protein samples**: If a protein sample has no matched RNA sample, assign it a `sample_name` that does not match any RNA sample. In this case, only the GENCODE reference proteome will be used as the proteomics search database.

**Samplesheet columns:**

- `sample_name`: Each RNA sample must have a distinct value. Do not include any spaces in this value.
- `sample_path`: Absolute path to the input file. 
  - RNA samples should be PacBio FLNC `.bam` or `.fastq` files
  - Protein samples should be `.raw` or `.mzML` files
- `condition`: Sample group (e.g., "control", "treatment"). Used for differential analysis, which performs pairwise comparisons between groups. Two or more groups are supported. If you do not want differential analysis, assign the same condition to all samples. Do not include any spaces in this value.
- `sample_type`: Must be either `RNA` or `protein`.
- `mass_spec_type`: Must be either `DDA` or `DIA`. Required for protein samples. For RNA samples, specify `none` for this column. 

## Running the Pipeline

For most datasets, we recommend running the pipeline from a driver shell script submitted to SLURM/LSF. This is more robust than running interactively (as shown in the [Quick Start](#quick-start) for test data), since large datasets may exceed interactive job time limits or resource quotas. 

### Recommended: Driver shell script (SLURM example)

Create a `run_lrp2.sh` script with your SLURM/LSF run parameters:

```bash
#!/bin/bash

#SBATCH --job-name=lrp2_driver
#SBATCH --partition=your_partition
#SBATCH --account=your_allocation
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G

module load nextflow apptainer

nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --dataset_name my_dataset \
    --genome GRCh38.p14.v49 \
    --protein_search fragpipe \
    --fragpipe_token "YOUR_TOKEN" \
    --hpc_partition standard \
    --hpc_cluster_options "--account=my_allocation" \
    -profile singularity,slurm
```

Submit with:

```bash
sbatch run_lrp2.sh
```
> **Note**: Customize the above template for your HPC. This includes `#SBATCH` directives (partition, account) and module names (nextflow, apptainer), and the `--hpc_partition` and `--hpc_cluster_options` pipeline parameters. For LSF, replace `#SBATCH` directives with `#BSUB` equivalents and use `-profile singularity,lsf`. For other schedulers, see [Support and customization](#support-and-customization).
>
> **Resource allocation** works on two levels:
> - **The driver job** (`#SBATCH` directives in the shell script): modest resources are sufficient — Nextflow itself only orchestrates submissions and doesn't run the heavy work.
> - **Individual pipeline tasks** (CPUs, memory, time per process): handled automatically by LRP2's internal configuration. You do **not** need to specify these on the command line. To customize, edit `conf/base.config`.
>
> Include `--fragpipe_token` only if running the proteomics subworkflow (see [Run the RNA + DDA proteomics test dataset](#run-the-rna--dda-proteomics-test-dataset) for obtaining a token). Differential analysis runs automatically when two or more conditions are present in the samplesheet.

## Profile Options

Nextflow profiles control how the pipeline executes. Multiple profiles can be combined using commas (e.g., `-profile test_rna,singularity,slurm`).

### Container Profiles (choose ONE)

| Profile | Description | Best For |
|---------|-------------|----------|
| `singularity` | Use Singularity/Apptainer containers | **HPC systems** (most common) |
| `docker` | Use Docker containers | Local machines with Docker installed |
| `conda` | Use Conda environments | Systems without container support (slower) |

> **Note**: We have extensively tested the pipeline with `singularity` locally and on HPC systems, and recommend its usage. You may use `docker` on local machines. We do not recommend the use of `conda` except as a last resort due to it lacking the same reproducibility as containers. 

### Executor Profiles (recommended)

| Profile | Description | When to Use |
|---------|-------------|-------------|
| `slurm` | Submit jobs to SLURM scheduler on an HPC | HPC environment |
| `lsf` | Submit jobs to LSF scheduler on an HPC | HPC environment |

> **Note**: When using `slurm` or `lsf`, Nextflow submits individual pipeline tasks as separate jobs. Without a scheduler profile, all tasks run on the node where you launch Nextflow (requires sufficient resources). If you intend to run locally, you may need to lower resource requirements in `conf/base.config`. 

### Test Profiles

| Profile | Description | Dataset |
|---------|-------------|---------|
| `test_rna` | RNA-only test dataset | Runs RNA subworkflows (S1 - S4) |
| `test_dda` | RNA + DDA proteomics test | Runs all subworkflows (S1 - S5) with FragPipe DDA search |

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

## Reference Genomes

The pipeline supports human and mouse data using GENCODE reference genomes across multiple versions:

- **Human**: `GRCh38.p14.v49`, `GRCh38.p14.v48`, ..., `GRCh38.p14.v44`, `GRCh38.p13.v43`, ..., `GRCh37.p13.v19`
- **Mouse**: `GRCm39.vM38`, `GRCm39.vM37`, `GRCm39.vM36`, `GRCm39.vM35`, `GRCm39.vM34`

The pipeline automatically downloads the appropriate FASTA and GTF files based on your `--genome` selection. Species is auto-detected from `--genome` and determines which CPAT model (human or mouse) is used for ORF prediction. See `conf/gencode_references.config` for the full list of supported versions.

Support for RefSeq / igenomes and custom references is under active development.

## Parameters

For a complete list of parameters:
```bash
nextflow run /path/to/LRP2 --help
```
> [!WARNING]
> Please provide pipeline parameters via the CLI as shown or using the Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration other than parameters.

### Input/Output

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to samplesheet CSV (required) | — |
| `--outdir` | Path to output directory (required) | — |
| `--dataset_name` | Run identifier used for output prefixes | `merged` |
| `--genome` | Reference genome version | `GRCh38.p14.v49` |

### HPC Scheduler Options

Use `--hpc_partition` for SLURM clusters and `--hpc_queue` for LSF clusters. SLURM defaults to Rivanna conventions — customize for your cluster.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--hpc_partition` | SLURM partition | `standard` |
| `--hpc_queue` | LSF queue name | — |
| `--hpc_cluster_options` | Additional scheduler-specific options (e.g., `--account=my_alloc` for SLURM, `-P my_project` for LSF) | — |

### S1 PacBio Isocall

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_read_support` | Minimum read support for transcripts | `3` |
| `--isocall_config` | Path to custom Isocall configuration TOML file | `bin/isocall_config.toml` |
  
### S2 Transcriptome

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--protein_coding_filter` | Keep only protein-coding genes | `true` |
| `--internal_priming_filter` | Remove internal priming artifacts | `true` |
| `--template_switching_filter` | Remove template switching artifacts | `true` |
| `--transcript_class_keep` | Structural categories to retain (FSM, ISM, NIC, NNC, ALL) | `FSM,ISM,NIC,NNC` |

### S3 Predicted Proteome

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_orf` | Minimum ORF length in nucleotides | `75` |
| `--cpat_coding_threshold` | Coding probability threshold | `0.364` (human), `0.44` (mouse) |
| `--protein_class_keep` | Protein categories to retain | `FPM,NPC,NPE` |

### S4 Multisample Analysis

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_samples_per_intron` | Minimum samples per intron for leafcutter | `2` |
| `--min_samples_per_group` | Minimum samples per group for leafcutter | `1` |
| `--min_usage_ratio` | Minimum junction usage ratio for filtering | `0.01` |
  
### S5 Proteomics

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--protein_search` | Search engine: `fragpipe` (required) | - |
| `--fragpipe_token` | Single-use academic license token for FragPipe (required if `--protein_search fragpipe`). See [Run the RNA + DDA proteomics test dataset](#run-the-rna--dda-proteomics-test-dataset) for how to obtain one. | — |
| `--fragpipe_workflow` | Path to a custom FragPipe workflow file specifying search parameters (modifications, enzymes, etc.) | default is selected by `mass_spec_type` |

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
│   └── M3_FILTER_TRANSCRIPTOME/         # GTF, BED12, DNA FASTA, count matrix of the refined transcriptome (technical artifacts removed)
├── S3_PREDICTED_PROTEOME/               
│   ├── M1_CPAT_ORF/                     
│   ├── M2_FILTER_CPAT/                  # GTF with exon and CDS type columns for single best ORF identified per transcript 
│   ├── M3_SQANTI_PROTEIN/
│   └── M4_PROTEIN_CLASSIFICATION/       # GTF, BED12, protein FASTA, count matrix collapsed to distinct ORFs 
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
│   ├── M1_BUILD_PROTEOME_REFERENCE/
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

- **Megan D. Schertzer**, Sheynkman Lab - Lead developer
- **Julia T. Lewandowski**, Knowles Lab - Lead developer

We thank the following people for their extensive assistance in the development of this pipeline: 
- **Emily F. Watts**, Sheynkman Lab - Contributions to LRP and conception of multi-sample analysis subworkflow.
- **Madison M. Mehlferber**, Sheynkman Lab - Contributor to the original LRP pipeline. Continued pipeline testing and feedback.
- **Will Rosenow**, Sheynkman Lab - Pipeline testing and feedback.
- **Scott I. Adamson**, Knowles Lab - Development of leafcutter-py. 
- **Jocelyne Bruand**, Pacific Biosciences - Development of Isocall. 
- **Elizabeth Tseng**, Pacific Biosciences - Development of Isocall. 
- **Egor Dolzhenko**, Pacific Biosciences - Lead Developer of Isocall. 

We especially thank the PIs that contributed to this project: 
- **David A. Knowles**, Development of LR LeafCutter and project support / funding 
- **Gloria Sheynkman**, Development/conceptualization of LRP and project support / funding

## Support and Customization

LRP2 supports a range of customization:

- **HPC environment**: SLURM and LSF schedulers with configurable partition/queue and cluster options; Singularity, Docker, or Conda containers
- **Input flexibility**: RNA-only, protein-only, or paired RNA + protein samples; DDA or DIA mass spec data
- **Reference genomes**: GENCODE human and mouse across multiple versions
- **Proteomics**: alternative search engines (FragPipe or MetaMorpheus) and customizable FragPipe workflows

We welcome input from the community — please reach out if you have a use case not covered by the defaults.

- **Issues and bug reports**: [GitHub Issues](https://github.com/sheynkman-lab/LRP2/issues)
- **Direct contact**: Megan Schertzer, cwp5au@virginia.edu

## License

This pipeline is released under the MIT License.

## Citations

### Citing LRP2

If you use LRP2 in your work, please cite:

> Schertzer MD, Lewandowski JT, et al. LRP2: A proteogenomics pipeline for long-read informed protein isoform analysis and discovery. *Manuscript in preparation.*

LRP2 builds on the original LRP framework:

> Miller, R. M., Jordan, B. T., Mehlferber, M. M., et al. 2022. "Enhanced protein isoform characterization through long-read proteogenomics." *Genome Biology* 23(1): 69. doi: [10.1186/s13059-022-02624-y](https://doi.org/10.1186/s13059-022-02624-y)

### Tool references

Please also cite the tools used by the pipeline:

- **Isocall** (PacBio) — [github.com/PacificBiosciences/isocall](https://github.com/PacificBiosciences/isocall)
  
- **SQANTI3**
  > Pardo-Palacios, F. J., Arzalluz-Luque, A., Kondratova, L., et al. 2024. "SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms." *Nature Methods* 21(5): 793–797. doi: [10.1038/s41592-024-02229-2](https://doi.org/10.1038/s41592-024-02229-2)

- **CPAT**
  > Wang, L., et al. 2013. "CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model." *Nucleic Acids Research* 41(6): e74. doi: [10.1093/nar/gkt006](https://doi.org/10.1093/nar/gkt006)

- **SQANTI Protein** (originally introduced in LRP)
  > Miller, R. M., Jordan, B. T., Mehlferber, M. M., et al. 2022. "Enhanced protein isoform characterization through long-read proteogenomics." *Genome Biology* 23(1): 69. doi: [10.1186/s13059-022-02624-y](https://doi.org/10.1186/s13059-022-02624-y)

- **edgeR**
  > Robinson, M. D., McCarthy, D. J., and Smyth, G. K. 2010. "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." *Bioinformatics* 26(1): 139–140. doi: [10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616)

- **DRIMSeq**
  > Nowicka, M., and Robinson, M. D. 2016. "DRIMSeq: a Dirichlet-multinomial framework for multivariate count outcomes in genomics." *F1000Research* 5: 1356. doi: [10.12688/f1000research.8900.2](https://doi.org/10.12688/f1000research.8900.2)

- **LeafCutter** (adapted as Long-read LeafCutter in LRP2; under active development)
  > Li, Y. I., Knowles, D. A., Humphrey, J., et al. 2017. "Annotation-free quantification of RNA splicing using LeafCutter." *Nature Genetics* 50(1): 151–158. doi: [10.1038/s41588-017-0004-9](https://doi.org/10.1038/s41588-017-0004-9)

- **FragPipe / MSFragger**
  > Kong, A. T., Leprevost, F. V., Avtonomov, D. M., Mellacheruvu, D., and Nesvizhskii, A. I. 2017. "MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry-based proteomics." *Nature Methods* 14(5): 513–520. doi: [10.1038/nmeth.4256](https://doi.org/10.1038/nmeth.4256)

- **IonQuant** (label-free quantification in FragPipe)
  > Yu, F., Haynes, S. E., and Nesvizhskii, A. I. 2021. "IonQuant enables accurate and sensitive label-free quantification with FDR-controlled match-between-runs." *Molecular & Cellular Proteomics* 20: 100077. doi: [10.1016/j.mcpro.2021.100077](https://doi.org/10.1016/j.mcpro.2021.100077)

- **MSFragger-DIA** (DIA proteomics in FragPipe)
  > Yu, F., Teo, G. C., Kong, A. T., et al. 2023. "Analysis of DIA proteomics data using MSFragger-DIA and FragPipe computational platform." *Nature Communications* 14(1): 4154. doi: [10.1038/s41467-023-39869-5](https://doi.org/10.1038/s41467-023-39869-5)

### Reference data

- **GENCODE**
  > Frankish, A., et al. 2023. "GENCODE: reference annotation for the human and mouse genomes in 2023." *Nucleic Acids Research* 51(D1): D942–D949. doi: [10.1093/nar/gkac1071](https://doi.org/10.1093/nar/gkac1071)

### Framework

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., and Nahnsen, S. 2020. "The nf-core framework for community-curated bioinformatics pipelines." *Nature Biotechnology* 38(3): 276–278. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)
