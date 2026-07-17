# LRP2: Long-Read Proteogenomics Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/) [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) [![run with conda](https://img.shields.io/badge/run%20with-conda-43b02a.svg?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**LRP2** is a scalable, end-to-end long-read proteogenomics pipeline built in Nextflow. It identifies and validates protein isoforms by integrating PacBio long-read RNA-seq with mass spectrometry. Starting from full-length non-chimeric (FLNC) reads and/or raw MS files, LRP2 performs transcript discovery, ORF prediction, differential analysis, and mass spectrometry-based protein identification.

## Pipeline Overview
![LRP2 Logo](assets/LRP2-Workflow.lightmode.drawio.png#gh-light-mode-only)
![LRP2 Logo](assets/LRP2-Workflow.darkmode.drawio.png#gh-dark-mode-only)

The LRP2 Nextflow pipeline consists of five subworkflows:

| Subworkflow | Description |
|-------------|-------------|
| 1. PacBio Isocall | Align FLNC reads and collapse to isoforms with PacBio Isocall |
| 2. Transcriptome | Classify transcripts with SQANTI3, filter artifacts, assign deterministic hash-based isoform IDs |
| 3. Predicted proteome | Predict ORFs with CPAT, classify proteins with SQANTI protein |
| 4. Multi-sample analysis | Differential expression and usage (edgeR, DRIMSeq), differential splicing (LR LeafCutter, preliminary implementation) |
| 5. Proteomics | Build custom reference database, convert raw MS files, search with FragPipe or MetaMorpheus, map peptides to isoforms |

## Getting Started

- **[Installation](getting-started/installation.md)**: Install Nextflow and container runtime
- **[System Requirements](getting-started/requirements.md)**: Hardware and software requirements
- **[Quick Start](getting-started/quickstart.md)**: Run the test dataset in minutes

## Documentation

- **[Usage Guide](usage/index.md)**: Detailed instructions for running the pipeline
- **[Output Files](output/structure.md)**: Understand LRP2's outputs and results structure
- **[Parameters](usage/parameters.md)**: Understand LRP2's parameters
- **[Support](support.md)**: How to get help and report issues

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

## License

This pipeline is released under the MIT License.

## Citing LRP2

If you use LRP2 in your work, please cite our [preprint](https://www.biorxiv.org/content/10.64898/2026.05.27.728216v1):

> Schertzer MD, Lewandowski JT, et al. LRP2: A proteogenomics pipeline for long-read informed protein isoform analysis and discovery. *Manuscript in preparation.*

LRP2 builds on the original LRP framework:

> Miller, R. M., Jordan, B. T., Mehlferber, M. M., et al. 2022. "Enhanced protein isoform characterization through long-read proteogenomics." *Genome Biology* 23(1): 69. doi: [10.1186/s13059-022-02624-y](https://doi.org/10.1186/s13059-022-02624-y)

See the [full citation list](reference/citations.md) for all tools used by the pipeline.
