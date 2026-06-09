#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sheynkmanlab/lrp2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sheynkmanlab/lrp2
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LRP2  } from './workflows/lrp2'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_lrp2_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_lrp2_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_lrp2_pipeline'
include { logColours              } from './subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// NOTE: params.fasta and params.gtf are set in nextflow.config after gencode config is loaded
// Users can override with --fasta and --gtf if needed

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SHEYNKMANLAB_LRP2 {

    take:
    samplesheet         // channel: combined samplesheet with both RNA and protein samples

    main:

    //
    // WORKFLOW: Run pipeline
    //
    LRP2 (
        samplesheet
    )
    emit:
    versions = LRP2.out.versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // Auto-detect multisample-only mode based on parameter presence
    //
    def is_multisample_only_mode = params.transcripts_gtf &&
                                    params.transcript_counts &&
                                    params.orf_counts &&
                                    params.multisample_metadata

    //
    // Validate input parameters
    //
    if (is_multisample_only_mode) {
        log.info "${logColours(params.monochrome_logs).blue}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).purple}[sheynkmanlab/lrp2]${logColours(params.monochrome_logs).blue} Running in MULTISAMPLE ONLY mode${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).blue}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${logColours(params.monochrome_logs).reset}"

        // Validate all four parameters are provided
        def required_params = [
            'transcripts_gtf': params.transcripts_gtf,
            'transcript_counts': params.transcript_counts,
            'orf_counts': params.orf_counts,
            'multisample_metadata': params.multisample_metadata
        ]

        required_params.each { param_name, param_value ->
            if (!param_value) {
                error("ERROR: --${param_name} is required for multisample-only mode. All four parameters (--transcripts_gtf, --transcript_counts, --orf_counts, --multisample_metadata) must be provided together.")
            }
        }

        // Validate files exist
        def files_to_check = [
            'transcripts_gtf': file(params.transcripts_gtf),
            'transcript_counts': file(params.transcript_counts),
            'orf_counts': file(params.orf_counts),
            'multisample_metadata': file(params.multisample_metadata)
        ]

        files_to_check.each { param_name, file_obj ->
            if (!file_obj.exists()) {
                error("ERROR: File not found for --${param_name}: ${file_obj}")
            }
        }

        // Validate metadata file format
        def metadata_file = file(params.multisample_metadata)
        def metadata_lines = metadata_file.text.split('\n')
        if (metadata_lines.size() < 2) {
            error("ERROR: Metadata file is empty or has no data rows: ${params.multisample_metadata}")
        }

        def header = metadata_lines[0].toLowerCase()
        // Accept either 'name' or 'sample_name' for sample id, and 'group' or 'condition' for experimental condition
        def has_name_column = header.contains('name') || header.contains('sample_name')
        def has_group_column = header.contains('group') || header.contains('condition')

        if (!has_name_column || !has_group_column) {
            def missing = []
            if (!has_name_column) missing << "'name' or 'sample_name'"
            if (!has_group_column) missing << "'group' or 'condition'"
            error("ERROR: Metadata file must have ${missing.join(' and ')} column(s). Found header: ${metadata_lines[0]}\nFile: ${params.multisample_metadata}")
        }

        log.info "${logColours(params.monochrome_logs).green} ${logColours(params.monochrome_logs).reset} All multisample-only input files validated successfully!"
        log.info "${logColours(params.monochrome_logs).dim}  - Transcripts GTF: ${params.transcripts_gtf}${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).dim}  - Transcript counts: ${params.transcript_counts}${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).dim}  - ORF counts: ${params.orf_counts}${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).dim}  - Sample metadata: ${params.multisample_metadata}${logColours(params.monochrome_logs).reset}"
        log.info ""

    } else {
        // Normal mode: require samplesheet
        if (!params.input) {
            error("ERROR: --input (samplesheet) is required when not running in multisample-only mode.\n\nTo run in multisample-only mode, provide all four parameters:\n  --transcripts_gtf\n  --transcript_counts\n  --orf_counts\n  --multisample_metadata")
        }
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    SHEYNKMANLAB_LRP2 (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

