# Parameters

For a complete list of parameters:

```bash
nextflow run /path/to/LRP2 --help
```

!!! warning "Parameter Configuration"
    Please provide pipeline parameters via the CLI as shown or using the Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration other than parameters.

## Input/Output

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to samplesheet CSV (required) | — |
| `--outdir` | Path to output directory (required) | — |
| `--dataset_name` | Run identifier used for output prefixes | `merged` |
| `--genome` | Reference genome version | `GRCh38.p14.v49` |

## HPC Scheduler Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--hpc_queue` | SLURM/LSF partition/queue name | — |
| `--hpc_cluster_options` | Additional scheduler-specific options (e.g., '-A my_alloc' for SLURM, '-P my_project' for LSF) | — |

### Example Usage

**SLURM with account requirement:**

```bash
--hpc_queue standard --hpc_cluster_options '-A my_allocation'
```

**LSF with project requirement:**

```bash
--hpc_queue normal --hpc_cluster_options '-P my_project'
```

## S1 PacBio Isocall

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_read_support` | Minimum read support for transcripts | `3` |
| `--isocall_config` | Path to custom Isocall configuration TOML file | `bin/isocall_config.toml` |

## S2 Transcriptome

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--protein_coding_filter` | Keep only protein-coding genes | `true` |
| `--internal_priming_filter` | Remove internal priming artifacts | `true` |
| `--template_switching_filter` | Remove template switching artifacts | `true` |
| `--transcript_class_keep` | Structural categories to retain (FSM, ISM, NIC, NNC, ALL) | `FSM,ISM,NIC,NNC` |

## S3 Predicted Proteome

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_orf` | Minimum ORF length in nucleotides | `75` |
| `--cpat_coding_threshold` | Coding probability threshold | `0.364` (human), `0.44` (mouse) |
| `--protein_class_keep` | Protein categories to retain | `FPM,NPC,NPE` |

## S4 Multisample Analysis

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_samples_per_intron` | Minimum samples per intron for leafcutter | `2` |
| `--min_samples_per_group` | Minimum samples per group for leafcutter | `1` |
| `--min_usage_ratio` | Minimum junction usage ratio for filtering | `0.01` |

## S5 Proteomics

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--protein_search` | Search engine: `fragpipe` (required) | - |
| `--fragpipe_token` | Single-use academic license token for FragPipe (required if `--protein_search fragpipe`). See [Quick Start](../getting-started/quickstart.md#step-1-get-a-fragpipe-academic-license-token) for how to obtain one. | — |
| `--fragpipe_workflow` | Path to a custom FragPipe workflow file specifying search parameters (modifications, enzymes, etc.) | default is selected by `mass_spec_type` |
