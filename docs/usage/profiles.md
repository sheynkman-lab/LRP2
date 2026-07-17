# Profile Options

Nextflow profiles control how the pipeline executes. Multiple profiles can be combined using commas (e.g., `-profile test_rna,singularity,slurm`).

## Container Profiles (choose ONE)

| Profile | Description | Best For |
|---------|-------------|----------|
| `singularity` | Use Singularity/Apptainer containers | **HPC systems** (most common) |
| `docker` | Use Docker containers | Local machines with Docker installed |
| `conda` | Use Conda environments | Systems without container support (slower) |

!!! warning "Conda Not Recommended"
    We have extensively tested the pipeline with `singularity` locally and on HPC systems, and recommend its usage. You may use `docker` on local machines. We do not recommend the use of `conda` except as a last resort due to it lacking the same reproducibility as containers.

## Executor Profiles (recommended)

| Profile | Description | When to Use |
|---------|-------------|-------------|
| `slurm` | Submit jobs to SLURM scheduler on an HPC | HPC environment with SLURM |
| `lsf` | Submit jobs to LSF scheduler on an HPC | HPC environment with LSF |

!!! note "Local Execution"
    When using `slurm` or `lsf`, Nextflow submits individual pipeline tasks as separate jobs. Without a scheduler profile, all tasks run on the node where you launch Nextflow (requires sufficient resources). If you intend to run locally, you may need to lower resource requirements in `conf/base.config`.

## Test Profiles

| Profile | Description | Dataset |
|---------|-------------|---------|
| `test_rna` | RNA-only test dataset | Runs RNA subworkflows (S1 - S4) |
| `test_dda` | RNA + DDA proteomics test | Runs all subworkflows (S1 - S5) with FragPipe DDA search |

## Example Profile Combinations

### HPC with SLURM (recommended for production)

```bash
-profile singularity,slurm
```


### Quick RNA test on HPC

```bash
-profile test_rna,singularity,slurm
```

### Local machine with Docker

```bash
-profile docker
```
