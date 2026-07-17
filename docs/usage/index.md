# Usage Overview

This section provides detailed instructions for running the LRP2 pipeline on your own data.

## Before You Start

Make sure you have completed the [installation](../getting-started/installation.md) and successfully run the [quick start](../getting-started/quickstart.md) test dataset.

## Running the Pipeline

There are two main approaches to running LRP2:

### 1. Interactive Job (Test Data Only)

Suitable for small test datasets. See the [Quick Start](../getting-started/quickstart.md) guide.

### 2. Driver Shell Script (Recommended for Production)

For most datasets, we recommend running the pipeline from a driver shell script submitted to SLURM/LSF. This is more robust than running interactively since large datasets may exceed interactive job time limits or resource quotas.

See [Running the Pipeline](running-pipeline.md) for detailed instructions.

## Key Topics

- **[Preparing Input Data](input-preparation.md)**: Format your samplesheet and prepare FLNC reads and mass spec files
- **[Running the Pipeline](running-pipeline.md)**: Execute LRP2 with production data
- **[Profile Options](profiles.md)**: Choose container and executor profiles
- **[Parameters](parameters.md)**: Complete reference for all pipeline parameters
- **[Multisample Re-runs](multisample-rerun.md)**: Re-run differential analysis with different settings

## Getting Help

If you encounter issues:

1. Check the [parameters documentation](parameters.md) to ensure correct configuration
2. Review the [output structure](../output/structure.md) to verify results
3. See the [support page](../support.md) for how to get help
