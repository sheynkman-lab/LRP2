/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CPAT_ORF                      } from '../../modules/local/cpat_orf/main'
include { FILTER_CPAT                   } from '../../modules/local/filter_cpat/main'
include { SQANTI_PROTEIN_CLASSIFICATION } from '../../modules/local/sqanti_protein_classification/main'
include { PROTEIN_UTR_CLASSIFICATION    } from '../../modules/local/protein_utr_classification/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREDICTED_PROTEOME {
    take:
    ch_samples                  // channel: [meta, [corrected_fasta, corrected_gtf, classification_filtered, hashids_filtered]]
    reference_gtf               // path: reference GENCODE GTF
    hexamer_file                // path: CPAT hexamer file (species-specific)
    logit_model                 // path: CPAT logit model (species-specific)
    filter_cpat_script          // path: R script for filtering CPAT output
    sqanti_protein_script       // path: Python script for SQANTI protein classification
    protein_class_script        // path: R script for 5'UTR protein classification

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Predict ORFs using CPAT
    //
    CPAT_ORF (
        ch_samples.map { meta, fasta, gtf, classification, hashids ->
            [meta, fasta] },
        hexamer_file,
        logit_model
    )
    ch_versions = ch_versions.mix(CPAT_ORF.out.versions)

    //
    // MODULE: Filter CPAT output and call best ORFs
    //
    FILTER_CPAT (
        CPAT_ORF.out.orf_prob
            .join(CPAT_ORF.out.orf_seqs, by: 0)
            .join(ch_samples, by: 0)
            .map { meta, orf_prob, orf_seqs, fasta, gtf, classification, hashids ->
                [meta, orf_prob, orf_seqs, fasta, gtf, classification] },
        reference_gtf,
        filter_cpat_script
    )
    ch_versions = ch_versions.mix(FILTER_CPAT.out.versions)

    //
    // MODULE: Run SQANTI3 protein classification
    //
    SQANTI_PROTEIN_CLASSIFICATION (
        FILTER_CPAT.out.cds_gtf
            .join(FILTER_CPAT.out.best_orfs_mapped, by: 0)
            .map { meta, cds_gtf, best_orfs ->
                [meta, cds_gtf, best_orfs] },
        reference_gtf,
        sqanti_protein_script
    )
    ch_versions = ch_versions.mix(SQANTI_PROTEIN_CLASSIFICATION.out.versions)

    //
    // MODULE: Custom protein and 5'UTR classification
    //
    PROTEIN_UTR_CLASSIFICATION (
        SQANTI_PROTEIN_CLASSIFICATION.out.protein_classification
            .join(FILTER_CPAT.out.cds_gtf, by: 0)
            .join(ch_samples, by: 0)
            .map { meta, protein_class, cds_gtf, fasta, gtf, classification, hashids ->
                [meta, protein_class, cds_gtf, fasta, hashids] },
        reference_gtf,
        protein_class_script
    )
    ch_versions = ch_versions.mix(PROTEIN_UTR_CLASSIFICATION.out.versions)

    emit:
    // CPAT outputs
    orf_prob                    = CPAT_ORF.out.orf_prob                                     // [meta, *_cpat.ORF_prob.tsv]
    orf_seqs                    = CPAT_ORF.out.orf_seqs                                     // [meta, *_cpat.ORF_seqs.fa]
    cpat_error                  = CPAT_ORF.out.error_log                                    // [meta, *_cpat.error]

    // FILTER_CPAT outputs
    all_orfs_mapped             = FILTER_CPAT.out.all_orfs_mapped                           // [meta, *_all_orfs_mapped.tsv]
    best_orfs_mapped            = FILTER_CPAT.out.best_orfs_mapped                          // [meta, *_best_orfs_mapped.tsv]
    cds_gtf                     = FILTER_CPAT.out.cds_gtf                                   // [meta, *_CDS.gtf]

    // SQANTI_PROTEIN outputs
    protein_classification      = SQANTI_PROTEIN_CLASSIFICATION.out.protein_classification  // [meta, *.sqanti_protein_classification.tsv]

    // PROTEIN_UTR_CLASSIFICATION outputs
    protein_all_isoforms        = PROTEIN_UTR_CLASSIFICATION.out.protein_all_isoforms       // [meta, *_protein_all_isoforms.txt]
    protein_high_confidence     = PROTEIN_UTR_CLASSIFICATION.out.protein_high_confidence    // [meta, *_protein_high_confidence.txt]
    protein_gtf                 = PROTEIN_UTR_CLASSIFICATION.out.protein_gtf                // [meta, *_protein_high_confidence.gtf]
    protein_fasta               = PROTEIN_UTR_CLASSIFICATION.out.protein_fasta              // [meta, *_protein_high_confidence.fa]
    hashids_orf                 = PROTEIN_UTR_CLASSIFICATION.out.hashids_orf                // [meta, *_hashids_with_cpm_ORF.txt]

    versions                    = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
