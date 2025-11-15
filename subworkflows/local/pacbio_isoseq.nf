/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PBTK_PBMERGE        } from '../../modules/nf-core/pbtk/pbmerge/main'
include { ISOSEQ_CLUSTER      } from '../../modules/nf-core/isoseq/cluster/main'
include { ISOSEQ_ALIGN        } from '../../modules/local/isoseq_align/main'
include { ISOSEQ_COLLAPSE     } from '../../modules/local/isoseq_collapse/main'

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
    // MODULE: Merge all FLNC BAMs per sample
    //
    PBTK_PBMERGE (
        ch_samples
    )
    ch_versions = ch_versions.mix(PBTK_PBMERGE.out.versions)

    //
    // MODULE: Cluster merged FLNC reads into consensus isoforms
    //
    ISOSEQ_CLUSTER (
        PBTK_PBMERGE.out.bam
    )
    ch_versions = ch_versions.mix(ISOSEQ_CLUSTER.out.versions)

    //
    // MODULE: Align clustered consensus isoforms to reference genome
    //
    ISOSEQ_ALIGN (
        ISOSEQ_CLUSTER.out.bam,
        reference_fasta
    )
    ch_versions = ch_versions.mix(ISOSEQ_ALIGN.out.versions)

    //
    // MODULE: Collapse redundant isoforms based on alignment
    //
    ISOSEQ_COLLAPSE (
        ISOSEQ_ALIGN.out.bam
            .join(PBTK_PBMERGE.out.bam, by: 0)
            .join(PBTK_PBMERGE.out.pbi, by: 0)
            .join(ISOSEQ_CLUSTER.out.pbi, by: 0)
            .map { meta, aligned_bam, flnc_bam, flnc_pbi, clustered_pbi ->
                [meta, aligned_bam, flnc_bam, clustered_pbi, flnc_pbi] }
    )
    ch_versions = ch_versions.mix(ISOSEQ_COLLAPSE.out.versions)

    emit:
    merged_bam       = PBTK_PBMERGE.out.bam                 // [meta, *.bam]
    clustered_bam    = ISOSEQ_CLUSTER.out.bam               // [meta, *.transcripts.bam]
    aligned_bam      = ISOSEQ_ALIGN.out.bam                 // [meta, *.aligned.bam]
    collapsed_gff    = ISOSEQ_COLLAPSE.out.gff              // [meta, *.collapsed.gff]
    collapsed_count  = ISOSEQ_COLLAPSE.out.flnc_count       // [meta, *.collapsed.flnc_count.txt]
    collapsed_group  = ISOSEQ_COLLAPSE.out.group_txt        // [meta, *.collapsed.group.txt]
    versions         = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

