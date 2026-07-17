# Quick Start

!!! note "HPC Setup"
    The Quick Start uses SLURM on a UVA Rivanna-style HPC setup, which is what the pipeline has been most extensively tested with. LSF is also supported via the `lsf` profile. For other schedulers, see [Support and customization](../support.md).

    Specifics vary by cluster — account flags (`-A`), partition names (`-p`), and module names (`nextflow`, `apptainer`) may differ on your system. Check with your HPC documentation or admin.

## On an HPC (SLURM example)

### Start a Persistent Session

Start a persistent terminal session so the pipeline keeps running if you lose your connection:

```bash
screen -S lrp2
```

!!! tip "Screen Usage"
    To detach from screen, press `Ctrl+A` then `D`. To reattach later: `screen -r lrp2`.

    Certain HPC systems (e.g. UVA Rivanna) only support `screen`. On systems that support it, you can use the tmux terminal multiplexer instead by running `tmux new -s lrp2`.

### Request an Interactive Job

Request an interactive job with enough resources for the test dataset:

```bash
salloc -c 4 --mem=64G -p your_slurm_partition -A your_allocation --time=4:00:00
```

!!! note "Adjust for Your System"
    Replace `your_slurm_partition` with your SLURM partition and `your_allocation` with your SLURM allocation group. UVA Rivanna users can substitute `ijob` for `salloc`. The `-c` (CPUs), `--mem` (memory), and `--time` values above are sufficient for the test dataset, but should be increased for larger datasets.

### Load Required Modules

```bash
module load nextflow apptainer
```

## Run the RNA-only Test Dataset

From the `LRP2` directory:

```bash
nextflow run . \
    -profile test_rna,singularity,slurm \
    --outdir test_rna_results \
    --hpc_queue your_queue
```

!!! note "Configuration Notes"
    - Replace `your_queue` with your cluster's partition (SLURM) or queue (LSF) name.
    - If your cluster requires additional scheduler options such as account strings or QOS flags, pass them via `--hpc_cluster_options` (e.g., `--hpc_cluster_options '-A your_allocation'` for SLURM, `--hpc_cluster_options '-P your_project'` for LSF).
    - LSF users: swap `slurm` for `lsf` in the profile.
    - See [HPC Scheduler Options](../usage/parameters.md#hpc-scheduler-options) for more details.

!!! tip "Running Locally"
    To run locally on your current node instead of submitting to SLURM, drop the `slurm` profile: `-profile test_rna,singularity`.

## Run the RNA + DDA Proteomics Test Dataset

### Step 1: Get a FragPipe Academic License Token

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

!!! warning "Token Expiration"
    Tokens expire quickly. You will need a new token for each run.

!!! note "Non-Academic Users"
    See [Support and customization](../support.md) for alternative options.

### Step 2: Check Your Email

Check your email for a 6-digit token.

### Step 3: Run the Test

From the `LRP2` directory:

```bash
nextflow run . \
    -profile test_dda,singularity,slurm \
    --outdir test_results_dda \
    --fragpipe_token "YOUR_TOKEN" \
    --hpc_queue your_queue
```

!!! note
    The `test_dda` profile automatically sets `--protein_search fragpipe` and `--fragpipe_license_accept true`, which are arguments you need to specify when you want to run FragPipe on your own data.

!!! tip "Running Locally"
    To run locally instead of submitting to SLURM, drop the `slurm` profile: `-profile test_dda,singularity`.

## Next Steps

After successfully running the test dataset, you are now ready to run the pipeline on your data! Please refer to the following documentation: 

- Learn about [preparing your own input data](../usage/input-preparation.md)
- Understand the [pipeline parameters](../usage/parameters.md)
- Explore and understand LRP2's [output files](../output/structure.md)
