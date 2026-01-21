# LRP2-Lite: Long-Read Proteogenomics Lite Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) [![run with conda](https://img.shields.io/badge/run%20with-conda-43b02a.svg?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**LRP2_Lite** is a bioinformatics pipeline for comprehensive long-read proteogenomics analysis of PacBio Iso-Seq data and mass spectrometry proteomics data. It takes full-length non-chimeric (FLNC) BAM files and/or raw MS files as input, performing isoform discovery, quality control, protein prediction, and peptide-spectrum matching, outputting high-confidence proteomes and detailed classification reports.

The pipeline combines state-of-the-art tools for long-read RNA sequencing analysis and mass spectrometry proteomics with custom filtering and classification scripts following best practices to identify novel protein isoforms from PacBio data and validate them with experimental proteomics evidence.

## Pipeline Summary

The LRP2_Lite pipeline consists of five major stages:

### 1. PacBio Iso-Seq Processing (`01_pacbio_isoseq`)
- Merge FLNC BAM files per sample (**pbtk pbmerge**)
- Cluster merged reads into consensus isoforms (**isoseq cluster**)
- Align clustered consensus isoforms to reference genome (**pbmm2**)
- Collapse redundant isoforms (**isoseq collapse**)

### 2. Transcript Quality Control and Filtering (`02_transcriptome`)
- Perform comprehensive quality control and classification (**sqanti_qc**)
- Filter transcripts (**filter_transcriptome**) by:
  - Protein-coding status
  - Internal priming artifacts
  - Template switching artifacts
  - Structural categories
- Generate CPM-normalized hash IDs for protein mapping (**filter_transcriptome**)

### 3. Protein Prediction and Classification (`03_predicted_proteome`)
- Predict open reading frames (**cpat_orf**)
- Filter and select best ORF predictions (**filter_cpat**)
- Classify predicted proteins (**sqanti_protein_classification**)
- Generate high-confidence protein sets with custom UTR classification (**protein_utr_classification**)

### 4. Multi-Sample Differential Analysis (`04_multisample_analysis`) *(Optional)*
- Perform differential splicing analysis (**leafcutter_longread**):
  - Long-read junction clustering with minicutter
  - Differential intron usage testing
  - Subisoform cluster identification
- Conduct differential expression and usage analysis (**differential_expression**):
  - Differential gene expression (DGE) with edgeR
  - Differential transcript expression (DTE) with edgeR
  - Differential ORF expression with edgeR
  - Differential transcript usage (DTU) with DRIMSeq
  - Differential ORF usage (DPU) with DRIMSeq

> **Note**: Differential usage analyses (DTU/DPU with DRIMSeq) require biological replicates. Datasets without replicates may lead to unexpected results with the pipeline currently.

### 5. Proteomics Analysis (`05_proteomics`) *(Optional)*
- Convert raw MS data to mzML format (**msconvert**)
  - Supports .raw mass spectrometry files
  - Applies peak picking for centroid data
- Perform peptide-spectrum matching (**MetaMorpheus**):
  - Supports database search against predicted or reference proteome
  - Generates peptide-spectrum match (PSM) tables
  - Outputs comprehensive protein search results and statistics

> **Note**: The PROTEOMICS subworkflow can only run when protein samples (sample_type='protein') are provided in the samplesheet. It will search against the predicted proteome from Stage 3 if RNA samples were processed, and otherwise use a provided reference protein database via `--gencode_protein_fasta`. 

### Key Features:
- **Full isoform resolution**: Leverages PacBio long reads for complete transcript characterization
- **Quality-based filtering**: Multi-stage artifact removal and quality control
- **Protein-level analysis**: ORF prediction and protein classification
- **Proteomics integration**: Mass spectrometry protein search and analysis with MetaMorpheus
- **Multi-omics support**: Combined, integrated RNA-level and protein-level analysis
- **Species support**: Human and mouse genomes (via iGenomes)
- **Reproducible**: Fully containerized with Docker/Singularity support

## Usage
> **Note**: If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

### Preparing Input Data

First, prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv:**

```csv
sample_name,sample_path,condition,replicate,sample_type
control_chr22,/path/to/sample_data/230801_pacbio_rbfox2_control_chr22.flnc.bam,control,rep1,RNA
treatment_chr22,/path/to/sample_data/230801_pacbio_rbfox2_RB-G5_chr22.flnc.bam,treatment,rep1,RNA
```

Each row represents a PacBio Iso-Seq FLNC (Full-Length Non-Chimeric) BAM file, OR a Mass Spectrometry protein .raw or .mzML file. The columns are:

**Required:**
- `sample_name`: Unique sample identifier (no spaces)
- `sample_path`: Absolute path to the sample (*.bam if RNA, *.raw or .mzML if protein)
- `sanple_type`: Data type of the sample, which may be either 'RNA' or 'protein'

**Optional (required for differential analysis):**
- `condition`: Sample condition or group (e.g., "control", "treatment")
- `replicate`: Replicate identifier (e.g., "rep1", "rep2")

> **Note**: The `condition` column (or alternatively `group`) is required if you want to run differential analysis with `--run_differential_analysis`. The column can also be named `group` instead of `condition`.


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

### Example Command with Differential Analysis

To enable multi-sample differential analysis, use the `--run_differential_analysis` flag and specify the labels for you control and treatment groups, respectively:

```bash
nextflow run /path/to/LRP2_lite \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --run_differential_analysis \
    --control_group control \
    --experimental_group treatment \
    -profile singularity,slurm \
    -resume
```

> **Important**: For differential analysis, your samplesheet must include a `condition` column (or `group` column) with at least two groups. Biological replicates are strongly recommended for robust statistical analysis.

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

**Differential Analysis (Optional):**
- `--run_differential_analysis`: Enable multi-sample differential analysis (default: `false`)
- `--control_group`: Control/reference group name from samplesheet `condition` column (required if enabled)
- `--experimental_group`: Experimental/treatment group name from samplesheet `condition` column (required if enabled)
- `--sample_metadata`: Path to sample metadata CSV (default: uses `--input` samplesheet)
- `--min_samples_per_intron`: Minimum samples per intron for leafcutter (default: `2`)
- `--min_samples_per_group`: Minimum samples per group for leafcutter (default: `1`)
- `--min_usage_ratio`: Minimum junction usage ratio for filtering (default: `0.01`)

**Proteomics Analysis (Optional):**
- `--gencode_protein_fasta`: Path to reference protein database FASTA (required for protein-only analysis without RNA samples)
- `--metamorpheus_config`: Path to MetaMorpheus TOML configuration file (default: uses built-in ``sample_data/SearchTask.toml``)
- `--msconvert_peak_picking`: Enable peak picking during mzML conversion (default: `true`)

> **Note**: When both RNA and protein samples are provided, the pipeline automatically uses the predicted proteome from Stage 3 as the search database. For protein-only analysis (no RNA samples), you must provide `--gencode_protein_fasta`.

For a complete list of parameters, run:

```bash
nextflow run /path/to/LRP2_lite --help
```

> [!WARNING]
> Please provide pipeline parameters via the CLI as shown or using the Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters.

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

The pipeline generates comprehensive outputs organized in the following directories:

### Output Directory Structure

```
<outdir>/
├── S1_PACBIO_ISOSEQ/                    # Stage 1: Iso-Seq processing results
│   ├── M1_ISOSEQ_MERGE/                 # Merged FLNC BAM files
│   ├── M2_ISOSEQ_CLUSTER/               # Clustered consensus isoforms
│   ├── M3_ISOSEQ_ALIGN/                 # Aligned clustered isoforms
│   └── M4_ISOSEQ_COLLAPSE/              # Collapsed isoform GFF files
├── S2_TRANSCRIPTOME/                    # Stage 2: Transcript QC and filtering
│   ├── M1_SQANTI_QC/                    # SQANTI3 classification reports
│   └── M2_FILTER_TRANSCRIPTOME/         # Filtered transcript sets
├── S3_PREDICTED_PROTEOME/               # Stage 3: Protein predictions
│   ├── M1_CPAT_ORF/                     # ORF predictions
│   ├── M2_FILTER_CPAT/                  # Filtered ORFs and CDS
│   ├── M3_SQANTI_PROTEIN_CLASSIFICATION/# Protein classifications
│   └── M4_PROTEIN_UTR_CLASSIFICATION/   # High-confidence proteins
├── S4_MULTISAMPLE_ANALYSIS/             # Stage 4: Differential analysis (optional)
│   ├── M1_LEAFCUTTER_LONGREAD/          # Differential splicing results
│   │   ├── *_intron_coords.txt
│   │   ├── *_exon_coords.txt
│   │   ├── *_subisoform_clusters.txt
│   │   ├── *_cluster_significance.txt
│   │   └── *_effect_sizes.txt
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
├── S5_PROTEOMICS/                       # Stage 5: Proteomics analysis (optional)
│   ├── M1_MSCONVERT_MZML/               # .raw -> .mzML file conversion
│   │   └── *.mzML
│   └── M2_METAMORPHEUS/                 # MetaMorpheus search results
│       ├── *.psmtsv                     # Peptide-spectrum match tables
│       └── results/                     # Detailed search results
└── pipeline_info/                       # Execution reports and logs
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

**Proteomics (if protein samples provided):**
- Converted mzML files from raw MS data
- Peptide-spectrum match (PSM) tables
- MetaMorpheus search results and statistics
- Identified peptides and proteins

**Quality Control:**
- SQANTI3 QC reports and plots
- FLNC read counts
- Gene discovery statistics
- Junction analysis

**Differential Analysis (if enabled):**
- **Leafcutter long-read splicing:**
  - Intron and exon coordinate tables
  - Subisoform cluster assignments
  - Differential splicing significance and effect sizes
  - PSL alignment files
- **Differential expression (edgeR):**
  - Gene, transcript, and ORF-level differential expression results
  - Raw and normalized CPM matrices
  - MD plots for quality assessment
  - LogFC, p-values, and FDR-corrected values
- **Differential usage (DRIMSeq, requires replicates):**
  - Transcript and ORF-level differential usage summaries
  - Gene-level and feature-level significance
  - Proportional usage estimates per group

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
