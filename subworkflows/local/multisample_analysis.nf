/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { LEAFCUTTER_LONGREAD    } from '../../modules/local/leafcutter_longread/main'
include { DIFFERENTIAL_EXPRESSION   } from '../../modules/local/differential_expression/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow MULTISAMPLE_ANALYSIS {
    take:
    ch_transcripts              // channel: [meta, filtered_gtf, transcript_counts]
    ch_orfs                     // channel: [meta, orf_counts]
    sample_metadata             // path: sample metadata CSV
    lr_leafcutter_script        // path: R script for leafcutter clustering
    leafcutter_ds_script        // path: Python script for leafcutter differential splicing
    multisample_script          // path: R script for multisample analysis
    min_samples_per_intron      // val: minimum samples per intron
    min_samples_per_group       // val: minimum samples per group
    min_usage_ratio             // val: minimum usage ratio
    drimseq_min_gene_expr       // val: minimum gene expression for DRIMSeq
    drimseq_min_isoform_prop    // val: minimum isoform proportion for DRIMSeq

    main:
    ch_versions = channel.empty()

    // Parse metadata to generate all pairwise comparisons
    def metadata_file = sample_metadata instanceof java.nio.file.Path ? sample_metadata.toFile() : new File(sample_metadata.toString())
    def metadata_content = metadata_file.text
    def lines = metadata_content.split('\n')
    def header_parts = lines[0].split(',')
    def group_idx = header_parts.findIndexOf { it.trim() == 'group' }
    def unique_conditions = [] as Set
    lines.drop(1).each { line ->
        if (line.trim()) {
            def parts = line.split(',')
            if (parts.size() > group_idx) {
                unique_conditions.add(parts[group_idx].trim())
            }
        }
    }

    def condition_pairs = []
    def conditions_list = unique_conditions.toList().sort()
    conditions_list.eachWithIndex { cond1, i ->
        conditions_list.eachWithIndex { cond2, j ->
            if (j > i) {
                condition_pairs << [cond1, cond2]
            }
        }
    }

    ch_comparison_pairs = channel.fromList(condition_pairs)
        .map { pair ->
            [pair[0], pair[1], "${pair[0]}_vs_${pair[1]}"]
        }

    ch_transcripts_for_comparisons = ch_transcripts
        .map { meta, gtf, counts -> [gtf, counts] }
        .combine(ch_comparison_pairs)
        .map { gtf, counts, control, experimental, comparison_id ->
            [[id: comparison_id], gtf, counts, control, experimental]
        }

    ch_transcripts_comparisons_split = ch_transcripts_for_comparisons
        .multiMap { meta, gtf, counts, control, experimental ->
            data: [meta, gtf, counts]
            control: control
            experimental: experimental
        }

    //
    // MODULE: Run long-read leafcutter clustering and differential splicing for each comparison
    //
    LEAFCUTTER_LONGREAD (
        ch_transcripts_comparisons_split.data,
        sample_metadata,
        lr_leafcutter_script,
        leafcutter_ds_script,
        ch_transcripts_comparisons_split.control,
        ch_transcripts_comparisons_split.experimental,
        min_samples_per_intron,
        min_samples_per_group,
        min_usage_ratio
    )
    ch_versions = ch_versions.mix(LEAFCUTTER_LONGREAD.out.versions)

    //
    // MODULE: Run differential expression and usage analysis for each comparison
    //
    ch_orfs_for_comparisons = ch_orfs
        .map { meta, orf_counts -> orf_counts }
        .combine(ch_comparison_pairs)
        .map { orf_counts, control, experimental, comparison_id ->
            [[id: comparison_id], orf_counts, control, experimental]
        }

    ch_counts_for_comparisons = ch_transcripts_for_comparisons
        .map { meta, gtf, counts, control, experimental -> [meta.id, meta, counts, control, experimental] }
        .join(
            ch_orfs_for_comparisons.map { meta, orf_counts, control, experimental -> [meta.id, orf_counts] },
            by: 0
        )
        .map { comparison_id, meta, transcript_counts, control, experimental, orf_counts ->
            [meta, transcript_counts, orf_counts, control, experimental]
        }

    ch_counts_comparisons_split = ch_counts_for_comparisons
        .multiMap { meta, transcript_counts, orf_counts, control, experimental ->
            data: [meta, transcript_counts, orf_counts]
            control: control
            experimental: experimental
        }

    DIFFERENTIAL_EXPRESSION (
        ch_counts_comparisons_split.data,
        sample_metadata,
        multisample_script,
        ch_counts_comparisons_split.control,
        ch_counts_comparisons_split.experimental,
        drimseq_min_gene_expr,
        drimseq_min_isoform_prop
    )
    ch_versions = ch_versions.mix(DIFFERENTIAL_EXPRESSION.out.versions)

    emit:
    // Leafcutter outputs
    psl                     = LEAFCUTTER_LONGREAD.out.psl                      // [meta, *.psl]
    intron_coords           = LEAFCUTTER_LONGREAD.out.intron_coords            // [meta, *_intron_coords.txt]
    exon_coords             = LEAFCUTTER_LONGREAD.out.exon_coords              // [meta, *_exon_coords.txt]
    subisoform_clusters     = LEAFCUTTER_LONGREAD.out.subisoform_clusters      // [meta, *_subisoform_clusters.txt]
    leafcutter_counts       = LEAFCUTTER_LONGREAD.out.counts_matrix            // [meta, *.counts.gz]
    groups_file             = LEAFCUTTER_LONGREAD.out.groups_file              // [meta, *_groups_file.txt]
    rplots                  = LEAFCUTTER_LONGREAD.out.rplots                   // [meta, Rplots.pdf]
    cluster_significance    = LEAFCUTTER_LONGREAD.out.cluster_significance     // [meta, *_cluster_significance.txt]
    effect_sizes            = LEAFCUTTER_LONGREAD.out.effect_sizes             // [meta, *_effect_sizes.txt]

    // Differential gene expression
    dge_results             = DIFFERENTIAL_EXPRESSION.out.dge_results             // [meta, *_DGE_edgeR_results.txt]
    dge_raw_cpm             = DIFFERENTIAL_EXPRESSION.out.dge_raw_cpm             // [meta, *_DGE_edgeR_raw_CPM_matrix.txt]
    dge_norm_cpm            = DIFFERENTIAL_EXPRESSION.out.dge_norm_cpm            // [meta, *_DGE_edgeR_normalized_CPM_matrix.txt]
    dge_plot                = DIFFERENTIAL_EXPRESSION.out.dge_plot                // [meta, *_DGE_MD_plot.pdf]

    // Differential transcript expression
    dte_results             = DIFFERENTIAL_EXPRESSION.out.dte_results             // [meta, *_DTE_edgeR_results.txt]
    dte_raw_cpm             = DIFFERENTIAL_EXPRESSION.out.dte_raw_cpm             // [meta, *_DTE_edgeR_raw_CPM_matrix.txt]
    dte_norm_cpm            = DIFFERENTIAL_EXPRESSION.out.dte_norm_cpm            // [meta, *_DTE_edgeR_normalized_CPM_matrix.txt]
    dte_plot                = DIFFERENTIAL_EXPRESSION.out.dte_plot                // [meta, *_DTE_MD_plot.pdf]

    // Differential ORF expression
    de_orf_results          = DIFFERENTIAL_EXPRESSION.out.de_orf_results          // [meta, *_DE_ORF_edgeR_results.txt]
    de_orf_raw_cpm          = DIFFERENTIAL_EXPRESSION.out.de_orf_raw_cpm          // [meta, *_DE_ORF_edgeR_raw_CPM_matrix.txt]
    de_orf_norm_cpm         = DIFFERENTIAL_EXPRESSION.out.de_orf_norm_cpm         // [meta, *_DE_ORF_edgeR_normalized_CPM_matrix.txt]
    de_orf_plot             = DIFFERENTIAL_EXPRESSION.out.de_orf_plot             // [meta, *_DE_ORF_MD_plot.pdf]

    // Differential usage
    dtu_summary             = DIFFERENTIAL_EXPRESSION.out.dtu_summary             // [meta, *_DTU_transcript_DRIMSeq_summary.txt]
    du_orf_summary          = DIFFERENTIAL_EXPRESSION.out.du_orf_summary          // [meta, *_DU_ORF_DRIMSeq_summary.txt]

    versions                = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
