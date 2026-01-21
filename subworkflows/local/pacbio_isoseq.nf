/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PBTK_PBMERGE as ISOSEQ_MERGE      } from '../../modules/nf-core/pbtk/pbmerge/main'
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
    ISOSEQ_MERGE (
        ch_samples
    )
    ch_versions = ch_versions.mix(ISOSEQ_MERGE.out.versions)

    //
    // MODULE: Cluster merged FLNC reads into consensus isoforms
    //
    ISOSEQ_CLUSTER (
        ISOSEQ_MERGE.out.bam
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
            .join(ISOSEQ_MERGE.out.bam, by: 0)
            .join(ISOSEQ_MERGE.out.pbi, by: 0)
            .join(ISOSEQ_CLUSTER.out.pbi, by: 0)
            .map { meta, aligned_bam, flnc_bam, flnc_pbi, clustered_pbi ->
                [meta, aligned_bam, flnc_bam, clustered_pbi, flnc_pbi] }
    )
    ch_versions = ch_versions.mix(ISOSEQ_COLLAPSE.out.versions)

    emit:
    // ISOSEQ_MERGE outputs
    merged_bam       = ISOSEQ_MERGE.out.bam                 // [meta, *.bam]
    merged_pbi       = ISOSEQ_MERGE.out.pbi                 // [meta, *.pbi]

    // ISOSEQ_CLUSTER outputs
    clustered_bam    = ISOSEQ_CLUSTER.out.bam               // [meta, *.transcripts.bam]
    clustered_pbi    = ISOSEQ_CLUSTER.out.pbi               // [meta, *.transcripts.bam.pbi]
    cluster          = ISOSEQ_CLUSTER.out.cluster           // [meta, *.transcripts.cluster]
    cluster_report   = ISOSEQ_CLUSTER.out.cluster_report    // [meta, *.transcripts.cluster_report.csv]
    transcriptset    = ISOSEQ_CLUSTER.out.transcriptset     // [meta, *.transcripts.transcriptset.xml]
    hq_bam           = ISOSEQ_CLUSTER.out.hq_bam            // [meta, *.transcripts.hq.bam]
    hq_pbi           = ISOSEQ_CLUSTER.out.hq_pbi            // [meta, *.transcripts.hq.bam.pbi]
    hq_fasta         = ISOSEQ_CLUSTER.out.hq_fasta          // [meta, *.transcripts.hq.fasta.gz]
    lq_bam           = ISOSEQ_CLUSTER.out.lq_bam            // [meta, *.transcripts.lq.bam]
    lq_pbi           = ISOSEQ_CLUSTER.out.lq_pbi            // [meta, *.transcripts.lq.bam.pbi]
    lq_fasta         = ISOSEQ_CLUSTER.out.lq_fasta          // [meta, *.transcripts.lq.fasta.gz]

    // ISOSEQ_ALIGN outputs
    aligned_bam      = ISOSEQ_ALIGN.out.bam                 // [meta, *.aligned.bam]
    aligned_bai      = ISOSEQ_ALIGN.out.bai                 // [meta, *.aligned.bam.bai]

    // ISOSEQ_COLLAPSE outputs
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

