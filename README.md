![LRP2 Logo](assets/LRP2-Workflow.lightmode.drawio.png#gh-light-mode-only)
![LRP2 Logo](assets/LRP2-Workflow.darkmode.drawio.png#gh-dark-mode-only)

# LRP2: Long-Read Proteogenomics Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) [![run with conda](https://img.shields.io/badge/run%20with-conda-43b02a.svg?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**LRP2** is a bioinformatics pipeline for comprehensive long-read proteogenomics analysis of PacBio Iso-Seq data and mass spectrometry proteomics data. It takes full-length non-chimeric (FLNC) BAM files and/or raw MS files as input, performing isoform discovery, quality control, protein prediction, and peptide-spectrum matching, outputting high-confidence proteomes and detailed classification reports.

The pipeline combines state-of-the-art tools for long-read RNA sequencing analysis and mass spectrometry proteomics with custom filtering and classification scripts following best practices to identify novel protein isoforms from PacBio data and validate them with experimental proteomics evidence.

## Pipeline Summary

The LRP2 pipeline consists of five major stages:

### 1. PacBio Long-Read Isoform Processing (`01_pacbio_isocall`)

**PacBio Isocall is used for long-read processing.** Isocall is PacBio's up-and-coming successor to Isoseq for isoform discovery and quantification. It offers excellent scalability to large datasets with drastically reduced runtime and improved accuracy in isoform detection.
- Align FLNC BAM files to reference genome (**pbmm2 align**)
- Profile aligned reads to create junction chain profiles per sample (**isocall profile**)
- Merge profiles from all samples (**isocall merge**)
- Prepare known isoforms database from reference GTF (**isocall prep-isoforms**)
- Call both known and novel isoforms using merged profiles (**isocall call**)

### 2. Transcript Quality Control and Filtering (`02_transcriptome`)
- Perform comprehensive quality control and classification (**sqanti_qc**)
- Generate unique hash identifiers for transcripts (**generate_hashids**):
  - Create hash IDs from junction chain coordinates
  - Output transcript-to-hash ID mapping files
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
- Perform peptide-spectrum matching with **FragPipe** or **MetaMorpheus**:
  - **FragPipe** (default, recommended and can run on DDA or DIA samples): Ultra-fast proteomics search pipeline
    - MSFragger for peptide identification
    - IonQuant for label-free quantification
    - DIA-NN integration for DIA data analysis
    - Supports both DDA and DIA acquisition modes
  - **MetaMorpheus** (alternative, only runnable on DDA samples): Comprehensive open-source proteomics tool
    - Database search against predicted or reference proteome
    - Generates peptide-spectrum match (PSM) tables
    - Outputs comprehensive protein search results and statistics

> **Note**: The PROTEOMICS subworkflow can only run when protein samples (sample_type='protein') are provided in the samplesheet. When both RNA and protein samples are present, it searches against a concatenated database of the predicted proteome from Stage 3 plus the reference protein FASTA. For protein-only samples, it uses only the reference protein database via `--protein_fasta` (or auto-detected from the GENCODE genome if specified). Select the search engine with `--protein_search` (options: `fragpipe` or `metamorpheus`, default: `fragpipe`). 

### Key Features:
- **Full isoform resolution**: Leverages PacBio long reads for complete transcript characterization
- **Quality-based filtering**: Multi-stage artifact removal and quality control
- **Protein-level analysis**: ORF prediction and protein classification
- **Proteomics integration**: Mass spectrometry protein search with FragPipe or MetaMorpheus
- **Multi-omics support**: Combined, integrated RNA-level and protein-level analysis
- **Flexible genome support**: RefSeq genomes (via iGenomes) and GENCODE genomes with multiple release versions
- **Reproducible**: Fully containerized with Docker/Singularity support

### Reference Genome Support

The pipeline supports three types of reference genome sources:

1. **RefSeq Genomes (via iGenomes)**: Standard genome builds using NCBI/Ensembl annotations from iGenomes (e.g. `--genome GRCh38`)
2. **GENCODE Genomes**: High-quality genome annotations from GENCODE, supporting multiple release versions (e.g. `--genome GRCh38.p14.v49`)
3. **Custom genome references**: You may pass paths to a custom, local FASTA and GTF file using the ``--fasta`` and ``--gencode_gtf`` parameters, respectively 

The pipeline automatically downloads and uses the appropriate FASTA and GTF files based on your `--genome` selection. GENCODE genomes are particularly useful when you need specific annotation versions or want the latest curated gene models.

## Usage
> **Note**: If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

### Preparing Input Data

First, prepare a samplesheet with your input data that looks as follows:

**samplesheet.csv:**

```csv
sample_name,sample_path,condition,sample_type,mass_spec_type
control_rep1,/path/to/control_rep1.flnc.bam,control,RNA,none
control_rep2,/path/to/control_rep2.flnc.bam,control,RNA,none
control_rep3,/path/to/control_rep3.flnc.bam,control,RNA,none
treatment_rep1,/path/to/treatment_rep1.flnc.bam,treatment,RNA,none
treatment_rep2,/path/to/treatment_rep2.flnc.bam,treatment,RNA,none
treatment_rep3,/path/to/treatment_rep3.flnc.bam,treatment,RNA,none
control_protein,/path/to/control_injection1.raw,control,protein,DDA
control_protein,/path/to/control_injection2.raw,control,protein,DDA
control_protein,/path/to/control_injection3.raw,control,protein,DDA
treatment_protein,/path/to/treatment_frac1.mzML,treatment,protein,DIA
treatment_protein,/path/to/treatment_frac2.mzML,treatment,protein,DIA
```

Each row with sample_type ``RNA`` represents a PacBio IsoSeq FLNC (Full-Length Non-Chimeric) ``.bam`` file.
Each row with sample_type ``protein`` represents a mass spectrometry sample containing protein data in either ``.raw`` or ``.mzML`` format.

**Required columns:**
- `sample_name`: Unique identifier for each biological replicate. **Each RNA sample must have a distinct `sample_name` value.** The sample names are used by Isocall to label count matrix columns. Protein samples with matching `sample_name` and `condition` will be grouped together for proteomics analysis. Do not include any spaces in this value.
- `sample_path`: Absolute or relative path to the sample file.
- `condition`: Sample condition or group (e.g., "control", "treatment"). Used for grouping samples and differential analysis. Do not include any spaces in this value.
- `sample_type`: Sample data type, which must be either ``RNA`` or ``protein``.

**Optional columns:**
- `mass_spec_type`: Mass spectrometry data acquisition type, which should be one of ``DDA`` or ``DIA``. This value is required for protein samples. For RNA samples, specify ``none``.


### Running the Pipeline

Now, you can run the pipeline using:

```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --genome GRCh38.p14.v46 \
    -profile singularity
```

This will launch the pipeline with the default parameters using Singularity containers. See below for more details about available profiles and parameters. 

If you are running on an HPC environment, you may also run using SLURM, e.g.:

```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --genome GRCh38.p14.v46 \
    -profile singularity,slurm
```

### Example Command with Additional Parameters

```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v46 \
    --species human \
    --protein_coding_filter true \
    --cpat_coding_threshold 0.364 \
    -profile singularity \
    -resume
```

### Example Command with Differential Analysis

To enable multi-sample differential analysis, use the `--run_differential_analysis` flag and specify the labels for you control and treatment groups, respectively:

```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v46 \
    --run_differential_analysis \
    --control_group control \
    --experimental_group treatment \
    -profile singularity,slurm \
    -resume
```

> **Important**: For differential analysis, your samplesheet must include a `condition` column (or `group` column) with at least two groups. Biological replicates are strongly recommended for robust statistical analysis.

### Example Command with Proteomics (FragPipe)

1. To use FragPipe, you must be an academic user and accept the license agreement for MSFragger, IonQuant, and diaTracer. Before using FragPipe in LRP2 for the first time, please review the license agreement here: https://msfragger.arsci.com/upgrader/LICENSE-ACADEMIC.pdf

Once you are ready to accept the license agreement, you may run the following curl command in your terminal with your information substituted for FIRST_NAME, LAST_NAME, EMAIL, and ACADEMIC_INSTITUTION_NAME: 

```bash
curl --location --request POST \
              'https://msfragger-upgrader.nesvilab.org/upgrader/upgrade_download.php' \
              --form 'transfer="academic"' \
              --form 'agreement2="true"' \
              --form 'agreement3="true"' \
              --form "first_name=FIRST_NAME" \
              --form "last_name=LAST_NAME" \
              --form "email=EMAIL" \
              --form "organization=ACADEMIC_INSTITUTION_NAME" \
              --form "download=4.4.1\$zip" \
              --form 'is_fragpipe="true"' \
              > /dev/null 2>&1
```

After submitting this command, you will receive a 6-digit token to your provided email, which is required to run command-line FragPipe in the pipeline in the next step. 

2) To then proceed with analyzing protein samples with FragPipe, you should include the required FragPipe registration parameters consisting of the same information (FIRST_NAME, LAST_NAME, EMAIL, ACADEMIC_INSTITUTION NAME), the 6-digit token you received to your email, and a confirmation that you accept the license agreement in your LRP2 run command:

```bash
nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v46 \
    --protein_search fragpipe \
    --fragpipe_first_name FIRST_NAME\
    --fragpipe_last_name LAST_NAME \
    --fragpipe_email EMAIL \
    --fragpipe_institution ACADEMIC_INSTITUTION_NAME \
    --fragpipe_token 123456 \
    --fragpipe_license_accept true \
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
- `--genome`: Reference genome version (default: `GRCh38.p14.v46`)
- `--fasta`: Path to reference FASTA (auto-set from iGenomes)
- `--gencode_gtf`: Path to GENCODE GTF annotation (auto-set from iGenomes)

**Filtering Options:**
- `--protein_coding_filter`: Keep only protein-coding genes (default: `true`)
- `--internal_priming_filter`: Remove internal priming artifacts (default: `true`)
- `--template_switching_filter`: Remove template switching artifacts (default: `true`)
- `--transcript_class_keep`: Comma-separated transcript structural categories to retain - FSM, ISM, NIC, NNC, or 'ALL' (default: `FSM,NIC,NNC`)

**Protein Prediction:**
- `--species`: Species for CPAT models - `human` or `mouse` (default: `human`)
- `--min_orf`: Minimum ORF length in nucleotides (default: `75`)
- `--top_orf`: Number of ORF candidates to report (default: `5`)
- `--cpat_coding_threshold`: Coding probability threshold (default: human=0.364, mouse=0.44)
- `--protein_class_keep`: Protein categories to retain (default: `FPM,IPM,NPC,NPE`)

**PacBio Long-Read Processing (ISOCALL):**
- `--min_read_support`: Minimum read support required for calling novel transcripts (default: `3`)
- `--max_bundles_per_gene`: Maximum bundles per gene during read sampling (default: `100`)
- `--isocall_config`: Path to Isocall configuration TOML file (default: `bin/isocall_config.toml`)

**Differential Analysis (Optional):**
- `--run_differential_analysis`: Enable multi-sample differential analysis (default: `false`)
- `--control_group`: Control/reference group name from samplesheet `condition` column (required if enabled)
- `--experimental_group`: Experimental/treatment group name from samplesheet `condition` column (required if enabled)
- `--sample_metadata`: Path to sample metadata CSV (default: uses `--input` samplesheet)
- `--min_samples_per_intron`: Minimum samples per intron for leafcutter (default: `2`)
- `--min_samples_per_group`: Minimum samples per group for leafcutter (default: `1`)
- `--min_usage_ratio`: Minimum junction usage ratio for filtering (default: `0.01`)

**Proteomics Analysis (Optional):**
- `--protein_search`: Proteomics search engine - `fragpipe` or `metamorpheus` (default: `fragpipe`)
- `--protein_fasta`: Path to reference protein database FASTA (optional - auto-detected from GENCODE genome if not provided)
- `--msconvert_peak_picking`: Enable peak picking during mzML conversion (default: `true`)

**FragPipe-specific parameters** (when `--protein_search fragpipe`):
- `--fragpipe_first_name`: User first name for FragPipe registration (required)
- `--fragpipe_last_name`: User last name for FragPipe registration (required)
- `--fragpipe_email`: User email for FragPipe registration (required)
- `--fragpipe_institution`: Institution name for FragPipe registration (required)
- `--fragpipe_token`: FragPipe license token (required - obtain from [fragpipe.nesvilab.org](https://fragpipe.nesvilab.org))
- `--fragpipe_license_accept`: Accept FragPipe academic license (required: `true`)

**MetaMorpheus-specific parameters** (when `--protein_search metamorpheus`):
- `--metamorpheus_config`: Path to MetaMorpheus TOML configuration file (default: uses built-in `sample_data/SearchTask.toml`)

> **Note**: When both RNA and protein samples are provided, the pipeline concatenates the predicted proteome from Stage 3 with the reference protein FASTA (either user-provided via `--protein_fasta` or auto-detected from GENCODE genome) to create a comprehensive search database. For protein-only analysis (no RNA samples), only the reference protein FASTA will be used.

For a complete list of parameters, run:

```bash
nextflow run /path/to/LRP2 --help
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
nextflow run /path/to/LRP2 \
    --input sample_data/samplesheet.csv \
    --outdir results \
    --genome GRCh38.p14.v46 \
    -profile singularity
```

> [!NOTE]
> The sample data is restricted to chromosome 22 for RNA samples and only a subset of peptides for mass spec samples to reduce file size and computation time, making it ideal for testing the pipeline and becoming familiar with the workflow. We highly recommend trying out the sample data prior to running with your own data!

## Test Data

The pipeline includes test data for quick validation of the installation and functionality. Test data is located in the `test_data/` directory and includes subsetted paired RNA and mass spec data for ENCODE4 K562 and HepG2 samples.

### Running RNA-only Test

To run the RNA-only test profile, you may run:

```bash
nextflow run sheynkmanlab/lrp2 \
    -profile test_rna,<docker/singularity> \
    --outdir test_results
```

This will process RNA samples from both K562 and HepG2 cell lines using the test configuration.

### Running RNA + DDA Mass Spec Test

To run the RNA + DDA mass spectrometry test profile with FragPipe, you must first obtain a FragPipe license token:

1. **Obtain FragPipe Token**: Run the following curl command with your information:

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

You will receive a 6-digit token via email.

2. **Run the test**: Use the token and your registration information:

```bash
nextflow run sheynkmanlab/lrp2 \
    -profile test_dda,<docker/singularity> \
    --outdir test_results_dda \
    --fragpipe_first_name YOUR_FIRST_NAME \
    --fragpipe_last_name YOUR_LAST_NAME \
    --fragpipe_email YOUR_EMAIL \
    --fragpipe_institution YOUR_INSTITUTION \
    --fragpipe_token 123456
```

Note: The `test_dda` profile automatically sets `--protein_search fragpipe` and `--fragpipe_license_accept true` so you do not need to specify these parameters.

## Pipeline Output

The pipeline generates comprehensive outputs organized in the following directories:

### Output Directory Structure

```
<outdir>/
├── S1_PACBIO_ISOCALL/                   # STAGE 1: Isocall processing results
│   ├── M1_ISOCALL_ALIGN/                # Align FLNC BAM files per sample in parallel
│   ├── M2_ISOCALL_PROFILE/              # Junction chain profiles per sample in parallel
│   ├── M3_ISOCALL_PREP/                 # Build known isoforms database from GTF
│   ├── M4_ISOCALL_MERGE/                # Merge profiles from all samples together
│   └── M5_ISOCALL_CALL/                 # Called isoforms using merged profiles and reference
├── S2_TRANSCRIPTOME/                    # STAGE 2: Transcript QC and filtering
│   ├── M1_SQANTI_QC/                    # SQANTI3 classification reports
│   └── M2_FILTER_TRANSCRIPTOME/         # Filtered transcript sets
├── S3_PREDICTED_PROTEOME/               # STAGE 3: Protein predictions
│   ├── M1_CPAT_ORF/                     # ORF predictions
│   ├── M2_FILTER_CPAT/                  # Filtered ORFs and CDS
│   ├── M3_SQANTI_PROTEIN_CLASSIFICATION/# Protein classifications
│   └── M4_PROTEIN_UTR_CLASSIFICATION/   # High-confidence proteins
├── S4_MULTISAMPLE_ANALYSIS/             # STAGE 4: Differential analysis (optional)
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
├── S5_PROTEOMICS/                       # STAGE 5: Proteomics analysis (optional)
│   ├── M1_MSCONVERT_MZML/               # .raw -> .mzML file conversion if needed
│   │   └── *.mzML
│   ├── M2_FRAGPIPE/                     # FragPipe search results (if --protein_search fragpipe)
│   │   ├── psm.tsv                      # Peptide-spectrum matches
│   │   ├── peptide.tsv                  # Identified peptides
│   │   ├── protein.tsv                  # Identified proteins
│   │   ├── ion.tsv                      # Ion-level quantification
│   │   └── combined_protein.tsv         # Combined protein quantification
│   └── M2_METAMORPHEUS/                 # MetaMorpheus search results (if --protein_search metamorpheus)
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

**Proteomics:**
- Converted mzML files from raw MS data
- **FragPipe outputs**:
  - Peptide-spectrum match (PSM) tables with confidence scores
  - Peptide-level identifications and quantification
  - Protein-level identifications and quantification
  - Ion-level quantification for label-free analysis
  - Combined protein reports across sample groups
- **MetaMorpheus outputs**:
  - Peptide-spectrum match (PSM) tables
  - Protein search results and statistics
  - Identified peptides and proteins

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

The LRP2 pipeline was developed through a collaboration by the Sheynkman Lab and Knowles Lab:

- **Megan D. Schertzer**, Sheynkman Lab - Module code development
- **Julia T. Lewandowski**, Knowles Lab - Pipeline implementation

We thank the following people for their extensive assistance in the development of this pipeline: 
- **Emily F. Watts**, Sheynkman Lab - Contributions to LRP and conception of multi-sample analysis paradigm
- **Madison M. Mehlferber**, Sheynkman Lab - Pipeline testing and feedback
- **Will Rosenow**, Sheynkman Lab - Pipeline testing and feedback
- **Scott I. Adamson**, Knowles Lab - Development of leafcutter-py 
- **Jocelyne Bruand**, Pacific Biosciences - Development of Isocall 
- **Elizabeth Tseng**, Pacific Biosciences - Development of Isocall 
- **Egor Dolzhenko**, Pacific Bioscience - Lead Developer of Isocall 

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
