/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SQANTI_QC              } from '../../modules/local/sqanti_qc/main'
include { GENERATE_HASHIDS       } from '../../modules/local/generate_hashids/main'
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
    filter_script               // path: R filtering script (02_filter_sqanti_transcripts.R)
    hashlib_script              // path: Python hashlib_id_generator script
    generate_hashids_script     // path: R script for generating hashids (00_generate_hashids.R)

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
    // MODULE: Generate hashids for all transcripts
    //
    GENERATE_HASHIDS (
        SQANTI_QC.out.corrected_gtf
            .join(SQANTI_QC.out.classification, by: 0),
        reference_gtf,
        hashlib_script,
        generate_hashids_script
    )
    ch_versions = ch_versions.mix(GENERATE_HASHIDS.out.versions)

    //
    // MODULE: Filter SQANTI3 output (requires hashids_mapping from GENERATE_HASHIDS)
    //
    FILTER_TRANSCRIPTOME (
        SQANTI_QC.out.classification
            .join(SQANTI_QC.out.corrected_gtf, by: 0)
            .join(SQANTI_QC.out.corrected_fasta, by: 0)
            .join(GENERATE_HASHIDS.out.hashids_mapping, by: 0)
            .join(GENERATE_HASHIDS.out.corrected_psl, by: 0)
            .map { meta, classification, gtf, fasta, hashids, psl ->
                [meta, classification, gtf, fasta, hashids, psl] },
        reference_gtf,
        reference_fasta,
        sample_metadata,
        filter_script
    )
    ch_versions = ch_versions.mix(FILTER_TRANSCRIPTOME.out.versions)

    emit:
    // SQANTI_QC outputs
    classification          = SQANTI_QC.out.classification                    // [meta, *_classification.txt]
    corrected_gtf           = SQANTI_QC.out.corrected_gtf                     // [meta, *_corrected.gtf]
    corrected_fasta         = SQANTI_QC.out.corrected_fasta                   // [meta, *_corrected.fasta]
//    corrected_genepred      = SQANTI_QC.out.corrected_genepred                // [meta, *_corrected.genePred]
//    corrected_cds_gff       = SQANTI_QC.out.corrected_cds_gff                 // [meta, *_corrected.gtf.cds.gff]
//    isoforms_gtf            = SQANTI_QC.out.isoforms_gtf                      // [meta, *.isoforms.gtf]
    junctions               = SQANTI_QC.out.junctions                         // [meta, *_junctions.txt]
//    params                  = SQANTI_QC.out.params                            // [meta, *.params.txt]
//    refannotation_genepred  = SQANTI_QC.out.refannotation_genepred            // [meta, refAnnotation_*.genePred]
    
    // GENERATE_HASHIDS outputs
    hashids_mapping         = GENERATE_HASHIDS.out.hashids_mapping               // [meta, *_transcriptome_hashids_mapping.txt]
    
    // FILTER_TRANSCRIPTOME outputs
//    transcriptome_classification = FILTER_TRANSCRIPTOME.out.classification           // [meta, *_transcriptome_classification.txt]
    classification_filtered = FILTER_TRANSCRIPTOME.out.classification_filtered       // [meta, *_classification_filtered.txt]
//    transcriptome_corrected_gtf = FILTER_TRANSCRIPTOME.out.corrected_gtf             // [meta, *_transcriptome_corrected.gtf]
    corrected_gtf_filtered  = FILTER_TRANSCRIPTOME.out.corrected_gtf_filtered        // [meta, *_corrected_filtered.gtf]
    corrected_bed_filtered  = FILTER_TRANSCRIPTOME.out.corrected_bed_filtered        // [meta, *_corrected_filtered.bed]
//    transcriptome_corrected_fasta = FILTER_TRANSCRIPTOME.out.corrected_fasta         // [meta, *_transcriptome_corrected.fasta]
    corrected_fasta_filtered= FILTER_TRANSCRIPTOME.out.corrected_fasta_filtered      // [meta, *_corrected_filtered.fasta]
//    corrected_psl           = FILTER_TRANSCRIPTOME.out.corrected_psl                 // [meta, *_transcriptome_corrected.psl]
//    hashids_mapping         = FILTER_TRANSCRIPTOME.out.hashids_mapping               // [meta, *_transcriptome_hashids_mapping.txt]
    hashids_filtered        = FILTER_TRANSCRIPTOME.out.hashids_filtered              // [meta, *_hashids_with_cpm_filtered.txt]
    hashids_all             = FILTER_TRANSCRIPTOME.out.hashids_all                   // [meta, *_all_hashids_with_cpm.txt]

    versions                = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
