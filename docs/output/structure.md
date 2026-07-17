# Pipeline Output Structure

Each subworkflow outputs to numbered module directories. The final module in each subworkflow typically contains the key results, while earlier modules contain intermediate files.

## Directory Structure

```
<outdir>/
├── S1_PACBIO_ISOCALL/
│   ├── M1_ISOCALL_ALIGN/
│   ├── M2_ISOCALL_PROFILE/
│   ├── M3_ISOCALL_PREP/
│   ├── M4_ISOCALL_MERGE/
│   └── M5_ISOCALL_CALL/                 # GTF of transcript structures and count matrix
├── S2_TRANSCRIPTOME/
│   ├── M1_SQANTI_QC/
│   ├── M2_GENERATE_HASHIDS/
│   └── M3_FILTER_TRANSCRIPTOME/         # GTF, BED12, DNA FASTA, count matrix of the refined transcriptome (technical artifacts removed)
├── S3_PREDICTED_PROTEOME/
│   ├── M1_CPAT_ORF/
│   ├── M2_FILTER_CPAT/                  # GTF with exon and CDS type columns for single best ORF identified per transcript
│   ├── M3_SQANTI_PROTEIN/
│   └── M4_PROTEIN_CLASSIFICATION/       # GTF, BED12, protein FASTA, count matrix collapsed to distinct ORFs
├── S4_MULTISAMPLE_ANALYSIS/             # (optional)
│   ├── M1_LEAFCUTTER_LONGREAD/          # Differential splicing results
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
├── S5_PROTEOMICS/                       # (optional)
│   ├── M1_BUILD_PROTEOME_REFERENCE/
│   ├── M2_MSCONVERT_MZML/
│   ├── M3_FRAGPIPE/
│   └── M4_NOVEL_PEPTIDES/               # BED12 of peptides mapped to genome, summary table of novel and annotated peptides mapped to isoforms
└── pipeline_info/                       # Execution reports and logs
    ├── execution_report.html
    ├── execution_timeline.html
    └── lrp2_software_versions.yml
```

## Key Output Locations

### S1: PacBio Isocall

**Location**: `S1_PACBIO_ISOCALL/M5_ISOCALL_CALL/`

- GTF of transcript structures
- Count matrix

### S2: Transcriptome

**Location**: `S2_TRANSCRIPTOME/M3_FILTER_TRANSCRIPTOME/`

- GTF, BED12, DNA FASTA
- Count matrix of the refined transcriptome (technical artifacts removed)

### S3: Predicted Proteome

**Intermediate**: `S3_PREDICTED_PROTEOME/M2_FILTER_CPAT/`
- GTF with exon and CDS type columns for single best ORF identified per transcript

**Final**: `S3_PREDICTED_PROTEOME/M4_PROTEIN_CLASSIFICATION/`
- GTF, BED12, protein FASTA
- Count matrix collapsed to distinct ORFs

### S4: Multisample Analysis (Optional)

**Location**: `S4_MULTISAMPLE_ANALYSIS/`

Includes differential expression/usage results:
- Gene expression (DGE)
- Transcript expression (DTE)
- ORF expression (DE_ORF)
- Transcript usage (DTU)
- ORF usage (DU_ORF)
- Differential splicing (LeafCutter)

### S5: Proteomics (Optional)

**Location**: `S5_PROTEOMICS/M4_NOVEL_PEPTIDES/`

- BED12 of peptides mapped to genome
- Summary table of novel and annotated peptides mapped to isoforms

## Pipeline Reports

**Location**: `pipeline_info/`

- `execution_report.html`: Resource usage and task statistics
- `execution_timeline.html`: Timeline visualization of task execution
- `lrp2_software_versions.yml`: Versions of all software used

For detailed information about specific output files, see [Output Files](files.md).
