/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PACBIO_ISOSEQ          } from '../subworkflows/local/pacbio_isoseq'
include { TRANSCRIPTOME          } from '../subworkflows/local/transcriptome'
include { PREDICTED_PROTEOME     } from '../subworkflows/local/predicted_proteome'
include { MULTISAMPLE_ANALYSIS   } from '../subworkflows/local/multisample_analysis'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_lrp2_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LRP2 {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Run PacBio IsoSeq analysis (merge, cluster, align, collapse)
    //
    PACBIO_ISOSEQ (
        ch_samplesheet,
        params.fasta
    )
    ch_versions = ch_versions.mix(PACBIO_ISOSEQ.out.versions)

    //
    // SUBWORKFLOW: Run SQANTI3 QC and filtering
    //
    // Note: samplesheet is used as sample_metadata if not explicitly provided
    sample_metadata_file = params.sample_metadata ?: params.input

    TRANSCRIPTOME (
        PACBIO_ISOSEQ.out.collapsed_gff
            .join(PACBIO_ISOSEQ.out.collapsed_count, by: 0)
            .map { meta, gff, count ->
                [meta, gff, count] },
        file(params.gencode_gtf),
        file(params.gencode_fasta),
        file(sample_metadata_file),
        file(params.filter_script),
        file(params.hashlib_script)
    )
    ch_versions = ch_versions.mix(TRANSCRIPTOME.out.versions)

    //
    // SUBWORKFLOW: Run predicted proteome analysis (ORF calling, CPAT filtering, protein classification)
    //
    // Determine species-specific CPAT files
    def hexamer_file = params.species == 'human' ?
        file(params.human_hexamer) : file(params.mouse_hexamer)
    def logit_model = params.species == 'human' ?
        file(params.human_logit_model) : file(params.mouse_logit_model)

    PREDICTED_PROTEOME (
        TRANSCRIPTOME.out.corrected_fasta_filtered
            .join(TRANSCRIPTOME.out.corrected_gtf_filtered, by: 0)
            .join(TRANSCRIPTOME.out.classification_filtered, by: 0)
            .join(TRANSCRIPTOME.out.hashids_filtered, by: 0)
            .map { meta, fasta, gtf, classification, hashids ->
                [meta, fasta, gtf, classification, hashids] },
        file(params.gencode_gtf),
        hexamer_file,
        logit_model,
        file(params.filter_cpat_script),
        file(params.sqanti_protein_script),
        file(params.protein_class_script)
    )
    ch_versions = ch_versions.mix(PREDICTED_PROTEOME.out.versions)

    // //
    // // SUBWORKFLOW: Run differential analysis (optional)
    // //
    if (params.run_differential_analysis) {
         // Prepare transcript channel with GTF and counts
         ch_transcripts = TRANSCRIPTOME.out.corrected_gtf_filtered
             .join(TRANSCRIPTOME.out.hashids_filtered, by: 0)

         // Prepare ORF counts channel
         ch_orfs = PREDICTED_PROTEOME.out.hashids_orf

         MULTISAMPLE_ANALYSIS (
             ch_transcripts,
             ch_orfs,
             file(sample_metadata_file),
             file(params.lr_leafcutter_script),
             file(params.multisample_script),
             params.control_group,
             params.experimental_group,
             params.min_samples_per_intron,
             params.min_samples_per_group,
             params.min_usage_ratio
         )
         ch_versions = ch_versions.mix(MULTISAMPLE_ANALYSIS.out.versions)
     }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'lrp2_software_versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

