# Running the Pipeline

For most datasets, we recommend running the pipeline from a driver shell script submitted to SLURM/LSF. This is more robust than running interactively (as shown in the [Quick Start](../getting-started/quickstart.md) for test data), since large datasets may exceed interactive job time limits or resource quotas.

## Recommended: Driver Shell Script

### SLURM Example

Create a `run_lrp2.sh` script with your SLURM run parameters:

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
    --hpc_queue your_queue \
    -profile singularity,slurm
```

Submit with:

```bash
sbatch run_lrp2.sh
```

### LSF Example

For LSF systems, use `#BSUB` directives instead:

```bash
#!/bin/bash

#BSUB -J lrp2_driver
#BSUB -q your_queue
#BSUB -P your_project
#BSUB -W 72:00
#BSUB -n 1
#BSUB -M 40G

module load nextflow apptainer

nextflow run /path/to/LRP2 \
    --input samplesheet.csv \
    --outdir results \
    --dataset_name my_dataset \
    --genome GRCh38.p14.v49 \
    --protein_search fragpipe \
    --fragpipe_token "YOUR_TOKEN" \
    --hpc_queue your_queue \
    --hpc_cluster_options '-P your_project' \
    -profile singularity,lsf
```

Submit with:

```bash
bsub < run_lrp2.sh
```

## Important Notes

!!! note "Resource Allocation"
    Resource allocation works on two levels:

    1. **The driver job** (`#SBATCH` or `#BSUB` directives in the shell script): modest resources are sufficient — Nextflow itself only orchestrates submissions and doesn't run the heavy work.

    2. **Individual pipeline tasks** (CPUs, memory, time per process): handled automatically by LRP2's internal configuration. You do **not** need to specify these on the command line. To customize, edit `conf/base.config`.

!!! tip "FragPipe Token"
    Include `--fragpipe_token` only if running the proteomics subworkflow. See [Quick Start](../getting-started/quickstart.md#step-1-get-a-fragpipe-academic-license-token) for obtaining a token.

!!! info "Differential Analysis"
    Differential analysis runs automatically when two or more conditions are present in the samplesheet.

## Customization

Customize the above template for your HPC:

- **SLURM/LSF directives**: partition/queue, account/project
- **Module names**: `nextflow`, `apptainer` (may be `singularity` on some systems)
- **Pipeline parameters**: `--hpc_queue` and `--hpc_cluster_options`

For other schedulers, see the [Support](../support.md) page.
