/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SQANTI_QC              } from '../../modules/local/sqanti_qc/main'
include { FILTER_TRANSCRIPTOME   } from '../../modules/local/filter_transcriptome/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow TRANSCRIPTOME {
    take:
    ch_samples                  // channel: [meta, [collapsed_gff, flnc_counts]]
    reference_gtf               // path: reference GENCODE GTF
    reference_fasta             // path: reference GENCODE FASTA
    sample_metadata             // path: sample metadata CSV
    filter_script               // path: R filtering script
    hashlib_script              // path: Python hashlib_id_generator script

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Run SQANTI3 QC on collapsed transcripts
    //
    SQANTI_QC (
        ch_samples,
        reference_gtf,
        reference_fasta
    )
    ch_versions = ch_versions.mix(SQANTI_QC.out.versions)

    //
    // MODULE: Filter SQANTI3 output
    //
    FILTER_TRANSCRIPTOME (
        SQANTI_QC.out.classification
            .join(SQANTI_QC.out.corrected_gtf, by: 0)
            .join(SQANTI_QC.out.corrected_fasta, by: 0)
            .map { meta, classification, gtf, fasta ->
                [meta, classification, gtf, fasta] },
        reference_gtf,
        reference_fasta,
        sample_metadata,
        filter_script,
        hashlib_script
    )
    ch_versions = ch_versions.mix(FILTER_TRANSCRIPTOME.out.versions)

    emit:
    classification          = SQANTI_QC.out.classification                    // [meta, *_classification.txt]
    corrected_gtf           = SQANTI_QC.out.corrected_gtf                     // [meta, *_corrected.gtf]
    corrected_fasta         = SQANTI_QC.out.corrected_fasta                   // [meta, *_corrected.fasta]
    junctions               = SQANTI_QC.out.junctions                         // [meta, *_junctions.txt]

    classification_filtered = FILTER_TRANSCRIPTOME.out.classification_filtered       // [meta, *_classification_filtered.txt]
    corrected_gtf_filtered  = FILTER_TRANSCRIPTOME.out.corrected_gtf_filtered        // [meta, *_corrected_filtered.gtf]
    corrected_fasta_filtered= FILTER_TRANSCRIPTOME.out.corrected_fasta_filtered      // [meta, *_corrected_filtered.fasta]
    hashids_filtered        = FILTER_TRANSCRIPTOME.out.hashids_filtered              // [meta, *_hashids_with_cpm_filtered.txt]
    hashids_all             = FILTER_TRANSCRIPTOME.out.hashids_all                   // [meta, *_all_hashids_with_cpm.txt]
    dropout_transcripts     = FILTER_TRANSCRIPTOME.out.dropout_transcripts           // [meta, dropout/*_dropout_transcripts.tsv]
    corrected_dropout_fasta = FILTER_TRANSCRIPTOME.out.corrected_dropout_fasta       // [meta, dropout/*_corrected_dropout.fasta]
    corrected_dropout_gtf   = FILTER_TRANSCRIPTOME.out.corrected_dropout_gtf         // [meta, dropout/*_corrected_dropout.gtf]

    versions                = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
