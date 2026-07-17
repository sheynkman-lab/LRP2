# System Requirements and HPC Recommendations

**LRP2 is primarily designed for High-Performance Computing (HPC) environments.** Although the test dataset can run on a local machine (minimum of 4 CPUs and 32GB RAM recommended), real-world datasets require substantial computational resources to be processed efficiently.

## Minimum Requirements for Test Dataset

- **CPUs**: 4 cores
- **Memory**: 32 GB RAM
- **Storage**: At least 50 GB free space for test data and intermediate files

## HPC Recommendations

For production runs with real-world datasets:

- **Scheduler**: SLURM or LSF
- **Container runtime**: Singularity/Apptainer
- **Persistent session**: screen or tmux for long-running jobs

### Resource Allocation

Resource allocation works on two levels:

1. **The driver job**: Modest resources are sufficient for the driver script, as these are used by Nextflow for task submission and pipeline orchestration. Typically 1 CPU and 40GB memory is adequate.

2. **Individual pipeline tasks**: CPUs, memory, and time per process are handled automatically by LRP2's internal configuration. You do **not** need to specify these on the command line. To customize, edit `conf/base.config` directly to set values for a process label (`process_low`, `process_medium`, `process_high`, etc.)

## Software Requirements

- **Nextflow** ≥ 24.04.2
- **Singularity/Apptainer** or **Docker**
- Standard Unix utilities (bash, git, curl)

## Network Requirements

The pipeline will download reference genomes (assuming you are using an inbuilt GENCODE reference with the `--genome` parameter) and container images on first run. Please ensure that your system has internet access, or that you have pre-staged these resources prior to running. 
