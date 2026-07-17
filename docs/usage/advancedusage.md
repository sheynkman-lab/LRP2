# Advanced Usage

This page covers advanced usage modes that allow you to re-run specific subworkflows of the LRP2 pipeline independently, using outputs from previous runs. These modes are useful for iterative analyses and parameter optimization.

---

## Re-running S4 Multisample Analysis Only

If you have already completed a full LRP2 run and want to re-run just the differential analysis modules (S4 MULTISAMPLE ANALYSIS) with different parameters, you can use **multisample-only mode** by passing a samplesheet through `--multisample_metadata`, as well as paths to three required input transcript and ORF-level output files from a previous run.

This allows you to skip subworkflows S1-S3 (PacBio Isocall, Transcriptome, Predicted Proteome) and rerun only S4 (Multisample Analysis), which is useful for testing different statistical thresholds, filtering parameters, and/or condition subgroupings of your samples. 

### Usage

```bash
nextflow run /path/to/LRP2 \
    --multisample_metadata samplesheet.csv \
    --transcripts_gtf results/S2_TRANSCRIPTOME/M3_FILTER_TRANSCRIPTOME/merged.transcriptome.filtered.gtf \
    --transcript_counts results/S2_TRANSCRIPTOME/M3_FILTER_TRANSCRIPTOME/merged.transcriptome.filtered_hashids_with_cpm.txt \
    --orf_counts results/S3_PREDICTED_PROTEOME/M4_PROTEIN_CLASSIFICATION/merged.predicted_proteome.collapsed_high_confidence_ORF_hashids_with_cpm.txt \
    --outdir results_reanalysis \
    --min_samples_per_intron 1 \
    --min_usage_ratio 0.05 \
    -profile singularity,slurm
```

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--multisample_metadata` | Path to samplesheet CSV. Must have columns: `sample_name`, `sample_path`, `condition` or `group`, and `sample_type`. |
| `--transcripts_gtf` | Path to filtered transcriptome GTF from a previous run |
| `--transcript_counts` | Path to transcript counts file from a previous run |
| `--orf_counts` | Path to ORF counts file from a previous run |
| `--outdir` | Output directory for the reanalysis |

### Metadata CSV Format

The metadata CSV passed to `--multisample_metadata` **must** have columns named `sample_name`, `sample_path`, `condition` or `group`, and `sample_type`.

!!! tip "Reusing Samplesheets"
    The standard samplesheet format will work for `--multisample_metadata`. However, if you are rerunning this for the same data with only different condition group labels for samples or minor parameter value changes, we recommend creating a unique samplesheet for each one and choosing a corresponding name for your output results directory to help keep your results organized and straightforward to differentiate between.

---

## Re-running S5 Proteomics Only

If you have previously completed the LRP2 pipeline through S3 (PREDICTED PROTEOME) and want to re-run only the proteomics analysis (S5 PROTEOMICS), you can use **proteomics-only mode** by providing the LRP2-generated protein FASTA and GTF from your previous run. This is useful for re-analyzing proteomics data with different search parameters and configurations,or if you would like to add new mass spectrometry samples while using the same LRP2 predicted proteome reference. 

### Usage

```bash
nextflow run /path/to/LRP2 \
    --input samplesheet_proteomics_only.csv \
    --lrp_protein_fasta results/S3_PREDICTED_PROTEOME/M4_PROTEIN_CLASSIFICATION/merged.predicted_proteome.best_ORF.fa \
    --lrp_gtf results/S3_PREDICTED_PROTEOME/M2_FILTER_CPAT/merged.predicted_proteome.best_ORF.gtf \
    --gencode_fasta gencode.v46.pc_translations.fa \
    --gencode_gtf gencode.v46.annotation.gtf \
    --outdir results_proteomics_reanalysis \
    -profile singularity,slurm
```

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to samplesheet CSV with **only proteomics (protein) samples** |
| `--lrp_protein_fasta` | Path to LRP2 protein FASTA from a previous S3 run |
| `--lrp_gtf` | Path to LRP2 CDS GTF from a previous S3 run  |
| `--gencode_fasta` | Path to GENCODE protein FASTA for reference proteome. |
| `--gencode_gtf` | Path to GENCODE annotation GTF. Required for novel peptide mapping if `--gencode_fasta` is provided. |
| `--outdir` | Output directory for the proteomics reanalysis |

!!! warning "Protein Samples Only"
    For proteomics-only mode to work correctly, your samplesheet must contain **only** `sample_type=protein` entries. Do not include any RNA samples. The pipeline will automatically detect the absence of RNA samples and skip S1-S4.

### Important Considerations

If using an LRP2 protein FASTA from S3 containing ALL predicted ORFs, keep in mind that proteomic search spaces will be larger, which affects false discovery rates. For more targeted proteomics searches, you may want to pre-filter the LRP2 FASTA. 

## Next Steps

For a list of all available parameters and their descriptions, see the [Parameters](parameters.md) documentation.
