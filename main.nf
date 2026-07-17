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

    }

    //
    // Auto-detect skip-isocall mode based on parameter presence
    //
    def skip_isocall_mode = params.external_gtf || params.external_counts

    if (skip_isocall_mode) {
        log.info "${logColours(params.monochrome_logs).blue}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).purple}[sheynkmanlab/lrp2]${logColours(params.monochrome_logs).blue} Running in SKIP ISOCALL mode (for non-PacBio data)${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).blue}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${logColours(params.monochrome_logs).reset}"

        // Validate both parameters are provided together
        if (!params.external_gtf || !params.external_counts) {
            error("ERROR: Both --external_gtf and --external_counts must be provided together to skip PACBIO_ISOCALL.\n\nProvided:\n  --external_gtf: ${params.external_gtf ?: 'NOT PROVIDED'}\n  --external_counts: ${params.external_counts ?: 'NOT PROVIDED'}\n\nTo use skip-isocall mode, provide both parameters. To run normally, omit both parameters.")
        }

        // Validate samplesheet is still required
        if (!params.input) {
            error("ERROR: --input (samplesheet) is still required in skip-isocall mode.\n\nThe samplesheet contains important metadata needed for downstream analysis.\nProvide a samplesheet with your RNA samples even when using external GTF/counts.")
        }

        // Validate files exist
        def external_gtf_file = file(params.external_gtf)
        def external_counts_file = file(params.external_counts)

        if (!external_gtf_file.exists()) {
            error("ERROR: External GTF file not found: ${params.external_gtf}")
        }
        if (!external_counts_file.exists()) {
            error("ERROR: External count matrix file not found: ${params.external_counts}")
        }

        // Validate GTF format (basic check)
        def gtf_content = external_gtf_file.name.endsWith('.gz') ?
            new java.util.zip.GZIPInputStream(new FileInputStream(external_gtf_file)).text :
            external_gtf_file.text

        def gtf_lines = gtf_content.split('\n').findAll { line -> !line.startsWith('#') && line.trim() }
        if (gtf_lines.size() == 0) {
            error("ERROR: External GTF file appears to be empty or contains only comments: ${params.external_gtf}")
        }

        // Validate count matrix format and sample name matching
        def counts_content = external_counts_file.text
        def counts_lines = counts_content.split('\n').findAll { it.trim() }

        if (counts_lines.size() < 2) {
            error("ERROR: External count matrix appears to be empty or has no data rows: ${params.external_counts}")
        }

        def counts_header = counts_lines[0]
        def counts_header_parts = counts_header.split(/[,\t]/)

        if (counts_header_parts.size() < 2) {
            error("ERROR: External count matrix must have at least 2 columns (transcript_id + at least one sample). Found: ${counts_header_parts.size()} column(s)\nHeader: ${counts_header}\nFile: ${params.external_counts}")
        }

        // Extract sample names from count matrix (skip first column which is transcript_id)
        def count_matrix_samples = counts_header_parts.drop(1).collect { it.trim() }

        // Parse samplesheet to get RNA sample names
        def samplesheet_file = file(params.input)
        def samplesheet_content = samplesheet_file.text
        def samplesheet_lines = samplesheet_content.split('\n')
        def samplesheet_header_parts = samplesheet_lines[0].split(',')
        def sample_name_idx = samplesheet_header_parts.findIndexOf { it.trim() == 'sample_name' }
        def sample_type_idx = samplesheet_header_parts.findIndexOf { it.trim() == 'sample_type' }

        if (sample_name_idx == -1) {
            error("ERROR: Samplesheet must have 'sample_name' column. Found header: ${samplesheet_lines[0]}")
        }

        def rna_sample_names = []
        samplesheet_lines.drop(1).each { line ->
            if (line.trim()) {
                def parts = line.split(',')
                if (parts.size() > sample_type_idx && sample_type_idx != -1) {
                    def sample_type = parts[sample_type_idx].trim().toLowerCase()
                    if (sample_type == 'rna') {
                        rna_sample_names << parts[sample_name_idx].trim()
                    }
                }
            }
        }

        if (rna_sample_names.size() == 0) {
            log.warn "WARNING: No RNA samples found in samplesheet. Count matrix validation skipped."
        } else {
            // Check that all count matrix samples are in the samplesheet
            def missing_in_samplesheet = count_matrix_samples.findAll { !rna_sample_names.contains(it) }
            def missing_in_counts = rna_sample_names.findAll { !count_matrix_samples.contains(it) }

            if (missing_in_samplesheet.size() > 0) {
                log.warn "WARNING: The following samples in count matrix are not in samplesheet RNA samples: ${missing_in_samplesheet.join(', ')}"
            }
            if (missing_in_counts.size() > 0) {
                error("ERROR: The following RNA samples in samplesheet are missing from count matrix: ${missing_in_counts.join(', ')}\n\nCount matrix samples: ${count_matrix_samples.join(', ')}\nRNA samples in samplesheet: ${rna_sample_names.join(', ')}")
            }
        }

        log.info "${logColours(params.monochrome_logs).green} ${logColours(params.monochrome_logs).reset} All skip-isocall input files validated successfully!"
        log.info "${logColours(params.monochrome_logs).dim}  - External GTF: ${params.external_gtf}${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).dim}  - External counts: ${params.external_counts}${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).dim}  - Count matrix samples (${count_matrix_samples.size()}): ${count_matrix_samples.join(', ')}${logColours(params.monochrome_logs).reset}"
        log.info "${logColours(params.monochrome_logs).dim}  - RNA samples in samplesheet (${rna_sample_names.size()}): ${rna_sample_names.join(', ')}${logColours(params.monochrome_logs).reset}"
        log.info ""
    }

    // Normal mode: require samplesheet (unless in multisample-only mode)
    if (!is_multisample_only_mode && !params.input) {
        error("ERROR: --input (samplesheet) is required when not running in multisample-only mode.\n\nTo run in multisample-only mode, provide all four parameters:\n  --transcripts_gtf\n  --transcript_counts\n  --orf_counts\n  --multisample_metadata")
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

