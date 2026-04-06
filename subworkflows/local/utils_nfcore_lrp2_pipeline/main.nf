//
// Subworkflow with functionality specific to the sheynkmanlab/lrp2 pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Sanitize samplesheet to handle Excel-created files
    //
    def sanitized_input = sanitizeSamplesheet(params.input)

    //
    // Create channel from input file provided through params.input
    //

    // Determine input type based on file content
    def input_type = detectInputType(sanitized_input)
    def rna_channel_data = []
    def protein_channel_data = []
    def combined_channel_data = []

    if (input_type == "sample_path") {
        def samplesheet_list = samplesheetToList(sanitized_input, "${projectDir}/assets/schema_input_bam.json")
        // Expected structure: [meta, sample_path] where meta contains: sample_name, condition, sample_type, mass_spec_type

        def rna_samples = samplesheet_list.findAll { row ->
            def sample_type = row[0].containsKey('sample_type') ? row[0].sample_type : 'unknown'
            sample_type.toLowerCase() == 'rna'
        }
        def protein_samples = samplesheet_list.findAll { row ->
            def sample_type = row[0].containsKey('sample_type') ? row[0].sample_type : 'unknown'
            sample_type.toLowerCase() == 'protein'
        }

        // =====================================================================
        // VALIDATE RNA SAMPLES
        // =====================================================================

        // Validate unique RNA sample_name values
        def rna_sample_names = rna_samples.collect { row -> row[0].sample_name }
        def duplicate_rna_names = rna_sample_names.findAll { name ->
            rna_sample_names.count(name) > 1
        }.unique()

        if (duplicate_rna_names.size() > 0) {
            error "ERROR: Duplicate RNA sample_name values detected: ${duplicate_rna_names.join(', ')}. " +
                  "Each RNA sample must have a unique sample_name. " +
                  "Please ensure all RNA samples have distinct sample_name values in your samplesheet."
        }

        // =====================================================================
        // SAMPLESHEET QC SUMMARY
        // =====================================================================
        def total_samples = samplesheet_list.size()
        def rna_count = rna_samples.size()
        def protein_count = total_samples - rna_count

        if (total_samples == 0) {
            error "ERROR: No samples found in samplesheet. Please provide at least one sample."
        }

        // Warn if no RNA samples (transcriptome/proteome prediction will be skipped)
        if (rna_count == 0 && protein_count > 0) {
            log.warn "WARNING: No RNA samples found in samplesheet. Only proteomics modules will run."
        } else if (rna_count == 0) {
            error "ERROR: No samples found in samplesheet. Please provide at least one RNA or protein sample."
        }

        // Count unique biological replicates (unique sample_names across all samples)
        def unique_sample_names = samplesheet_list.collect { row -> row[0].sample_name }.unique()

        // Count unique conditions
        def all_conditions = samplesheet_list.collect { row ->
            row[0].containsKey('condition') ? row[0].condition : 'unknown'
        }
        def unique_conditions = all_conditions.unique()

        // Calculate protein-only grouping statistics (protein samples grouped by sample_name + condition)
        def protein_groups = protein_samples.groupBy { row ->
            def sample_name = row[0].sample_name
            def condition = row[0].condition
            return "${sample_name}_${condition}"
        }

        // Display QC Summary with colors
        def colors = logColours(params.monochrome_logs)
        log.info ""
        log.info "${colors.blue}======================================================================${colors.reset}"
        log.info "${colors.blue}               SAMPLESHEET QC SUMMARY${colors.reset}"
        log.info "${colors.blue}======================================================================${colors.reset}"
        log.info ""
        log.info "${colors.cyan}Sample Overview:${colors.reset}"
        log.info "  ${colors.dim}Total files detected              :${colors.reset} ${colors.green}${total_samples}${colors.reset}"
        log.info "  ${colors.dim}Total RNA files                   :${colors.reset} ${colors.green}${rna_count}${colors.reset}"
        log.info "  ${colors.dim}Total protein files               :${colors.reset} ${protein_count > 0 ? colors.yellow : colors.dim}${protein_count}${colors.reset}"
        log.info "  ${colors.dim}Unique biological replicates      :${colors.reset} ${colors.green}${unique_sample_names.size()}${colors.reset}"
        log.info ""
        log.info "${colors.cyan}Protein Sample Overview:${colors.reset}"
        log.info "  ${colors.dim}Total unique sample groups       :${colors.reset} ${colors.green}${protein_groups.size()}${colors.reset}"
        protein_groups.each { group_key, group_samples ->
            def count = group_samples.size()
            log.info "    ${colors.purple}-${colors.reset} ${colors.bold}${group_key.padRight(30)}${colors.reset} : ${colors.green}${count}${colors.reset} file(s)"
        }
        log.info ""
        log.info "${colors.cyan}Sample Conditions Summary:${colors.reset}"
        log.info "  ${colors.dim}Total unique conditions detected  :${colors.reset} ${colors.green}${unique_conditions.size()}${colors.reset}"
        unique_conditions.each { condition ->
            log.info "    ${colors.purple}-${colors.reset} ${colors.bold}${condition}${colors.reset}"
        }
        log.info ""
        log.info "${colors.blue}======================================================================${colors.reset}"
        log.info ""

        // Create RNA channel data - keep individual samples separate for PACBIO_ISOCALL
        // Individual samples are needed for per-sample alignment and profiling
        if (rna_samples.size() > 0) {
            // Each sample needs its own meta and path for parallel processing in ISOCALL_ALIGN/PROFILE
            rna_channel_data = rna_samples.collect { row ->
                def sample_name = row[0].sample_name
                def condition = row[0].containsKey('condition') ? row[0].condition : 'unknown'
                def sample_path = row[1]

                def meta = [
                    id: sample_name,
                    sample_name: sample_name,
                    condition: condition,
                    sample_type: 'RNA'
                ]

                return [meta, sample_path]
            }
        }

        // Populate protein channel data - keep as individual samples with sample_type = 'protein'
        // Add group_key to metadata for matching with RNA samples later
        protein_channel_data = protein_samples.collect { row ->
            def meta = row[0]
            def sample_name = meta.sample_name
            def condition = meta.condition
            def group_key = "${sample_name}_${condition}"

            // Add group_key to metadata for later (need for matching with RNA groups)
            meta = meta + [
                sample_type: 'protein',
                group_key: group_key,
                id: "${meta.sample_name}_${new File(row[1].toString()).name.replaceAll(/\.(raw|mzML|mzml)$/, '')}"
            ]
            return [meta, row[1]]
        }

        // Combine RNA and protein channel data into a single list
        if (rna_channel_data) {
            combined_channel_data.addAll(rna_channel_data)
        }
        if (protein_channel_data) {
            combined_channel_data.addAll(protein_channel_data)
        }

    } else {
        // For FASTQ input, process as before but store in combined array
        combined_channel_data = samplesheetToList(sanitized_input, "${projectDir}/assets/schema_input.json")
            .collect { meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
            }
            .groupBy { item -> item[0] }
            .collect { entry ->
                def sample_id = entry.key
                def runs = entry.value
                def metas = runs.collect { run -> run[1] }
                def fastqs = runs.collect { run -> run[2] }
                validateInputSamplesheet([sample_id, metas, fastqs])
            }
            .collect { meta, fastqs ->
                return [ meta, fastqs.flatten() ]
            }
        // For FastQ input, all samples are RNA
    }

    emit:
    samplesheet = combined_channel_data.isEmpty() ? channel.empty() : channel.fromList(combined_channel_data)
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                null  // No MultiQC report - pipeline doesn't use MultiQC
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Sanitize samplesheet to handle Excel-created files (cleanup Windows line endings (\r) and BOM chars that get populated)
//
def sanitizeSamplesheet(input_file) {
    def file = new File(input_file)
    def content = file.getText('UTF-8')

    def has_bom = content.startsWith('\uFEFF')
    def has_windows_endings = content.contains('\r')

    // If the file is already clean, return original path
    if (!has_bom && !has_windows_endings) {
        return input_file
    }

    // otherwise file needs sanitization - create sanitized version
    log.info "Sanitizing samplesheet: Detected encoding issues (BOM: ${has_bom}, Windows line endings: ${has_windows_endings})"

    // Remove BOM if present (UTF-8 BOM: EF BB BF)
    if (has_bom) {
        content = content.substring(1)
    }
    content = content.replaceAll('\r\n', '\n')  // Remove Windows line endings (\r\n -> \n)
    content = content.replaceAll('\r', '\n')    // Remove any remaining carriage returns

    // Write sanitized content to new samplesheet (need this because sampleseetToList() expects a filepath and not raw content)
    def sanitized_file = new File("${input_file}.sanitized")
    sanitized_file.write(content, 'UTF-8')
    log.info "Sanitized samplesheet created at: ${sanitized_file}"

    return sanitized_file.toString()
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet (FastQ)
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Detect input type based on samplesheet headers
//
def detectInputType(input_file) {
    def file = new File(input_file)
    def firstLine = file.readLines()[0]

    if (firstLine.contains("sample_path") && firstLine.contains("condition")) {
        return "sample_path"
    } else if (firstLine.contains("fastq")) {
        return "fastq"
    } else {
        error("Could not determine input type from samplesheet headers. Expected either 'sample_path' or 'fastq' columns.")
    }
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// ANSI colours used for terminal logging
//
def logColours(monochrome_logs=true) {
    def colorcodes = [:] as Map

    // Reset / Meta
    colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']  = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']    = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']  = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow'] = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']   = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple'] = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']   = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']  = monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']  = monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']    = monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']  = monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow'] = monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']   = monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple'] = monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']   = monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']  = monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']  = monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']    = monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']  = monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow'] = monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']   = monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple'] = monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']   = monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']  = monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']  = monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']    = monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']  = monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow'] = monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']   = monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple'] = monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']   = monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']  = monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']  = monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']    = monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']  = monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow'] = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']   = monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple'] = monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']   = monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']  = monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}
