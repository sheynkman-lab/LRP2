/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PBTK_PBMERGE as ISOSEQ_MERGE_ALIGNED    } from '../../modules/nf-core/pbtk/pbmerge/main'
include { PBTK_PBMERGE as ISOSEQ_MERGE_FLNC } from '../../modules/nf-core/pbtk/pbmerge/main'
include { ISOSEQ_ALIGN                          } from '../../modules/local/isoseq_align/main'
include { ISOSEQ_COLLAPSE                       } from '../../modules/local/isoseq_collapse/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PACBIO_ISOSEQ {
    take:
    ch_samples              // channel: [meta, [flnc_bams]]
    reference_fasta         // path: reference genome FASTA

    main:
    ch_versions = channel.empty()

    //
    // Flatten individual FLNC BAMs for parallel alignment with unique meta IDs
    // [meta, [bam1, bam2]] -> [meta_sample1, bam1], [meta_sample2, bam2]
    // note: each gets unique ID from sample_names array to avoid filename collisions
    //
    ch_individual_flnc = ch_samples
        .flatMap { meta, bams ->
            // Create individual [meta, bam] items with unique IDs from sample_names array
            bams.withIndex().collect { bam, index ->
                def individual_meta = [
                    id: meta.sample_names[index],
                    dataset_id: meta.id,
                    sample_names: meta.sample_names,
                    bam_ids: meta.bam_ids,
                    conditions: meta.conditions,
                    replicates: meta.replicates,
                    sample_types: meta.sample_types
                ]
                [individual_meta, bam]
            }
        }

    //
    // MODULE 1: Align each individual FLNC BAM to reference genome in parallel
    //
    ISOSEQ_ALIGN (
        ch_individual_flnc,
        reference_fasta
    )
    ch_versions = ch_versions.mix(ISOSEQ_ALIGN.out.versions)

    //
    // Collect aligned BAMs back together by dataset_id for merging
    // [meta_sample1, aligned1], [meta_sample2, aligned2] -> [dataset_id, [aligned1, aligned2], meta_for_dataset]
    //
    ch_aligned_grouped = ISOSEQ_ALIGN.out.bam
        .map { meta, bam ->
            [meta.dataset_id, bam, meta]  // [dataset_id, bam, meta]
        }
        .groupTuple(by: 0)  // group by dataset_id
        .map { dataset_id, bams, metas ->
            // Recreate the dataset-level meta using info from first sample
            def dataset_meta = [
                id: dataset_id,  // e.g., "merged"
                sample_names: metas[0].sample_names,
                bam_ids: metas[0].bam_ids,
                conditions: metas[0].conditions,
                replicates: metas[0].replicates,
                sample_types: metas[0].sample_types
            ]
            [dataset_meta, bams]  // [meta, [aligned1, aligned2, ...]]
        }

    //
    // MODULE 2: Merge aligned BAMs
    //
    ISOSEQ_MERGE_ALIGNED (
        ch_aligned_grouped
    )
    ch_versions = ch_versions.mix(ISOSEQ_MERGE_ALIGNED.out.versions)

    //
    // Merge original FLNC BAMs (needed internally for isoseq collapse flnc_count.txt output file to ave correct FLNC counts)
    //
    ISOSEQ_MERGE_FLNC (
        ch_samples
    )
    ch_versions = ch_versions.mix(ISOSEQ_MERGE_FLNC.out.versions)

    //
    // MODULE 3: Collapse redundant isoforms based on merged alignment
    //
    ISOSEQ_COLLAPSE (
        ISOSEQ_MERGE_ALIGNED.out.bam
            .join(ISOSEQ_MERGE_FLNC.out.bam, by: 0)
            .join(ISOSEQ_MERGE_ALIGNED.out.pbi, by: 0)
            .join(ISOSEQ_MERGE_FLNC.out.pbi, by: 0)
            .map { meta, aligned_bam, flnc_bam, aligned_pbi, flnc_pbi ->
                [meta, aligned_bam, flnc_bam, aligned_pbi, flnc_pbi] }
    )
    ch_versions = ch_versions.mix(ISOSEQ_COLLAPSE.out.versions)

    emit:
    // ISOSEQ_MERGE_ALIGNED outputs (M2)
    merged_aligned_bam = ISOSEQ_MERGE_ALIGNED.out.bam         // [meta, *.bam] - merged aligned BAMs
    merged_aligned_pbi = ISOSEQ_MERGE_ALIGNED.out.pbi         // [meta, *.pbi]

    // ISOSEQ_COLLAPSE outputs (M3)
    collapsed_gff    = ISOSEQ_COLLAPSE.out.gff              // [meta, *.collapsed.gff]
    collapsed_fasta  = ISOSEQ_COLLAPSE.out.fasta            // [meta, *.collapsed.fasta]
    collapsed_abundance = ISOSEQ_COLLAPSE.out.abundance     // [meta, *.collapsed.abundance.txt]
    collapsed_count  = ISOSEQ_COLLAPSE.out.flnc_count       // [meta, *.collapsed.flnc_count.txt]
    collapsed_group  = ISOSEQ_COLLAPSE.out.group_txt        // [meta, *.collapsed.group.txt]
    collapsed_read_stat = ISOSEQ_COLLAPSE.out.read_stat     // [meta, *.collapsed.read_stat.txt]
    collapsed_report = ISOSEQ_COLLAPSE.out.report           // [meta, *.collapsed.report.json]

    versions         = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

