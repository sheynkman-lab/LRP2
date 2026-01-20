/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PACBIO_ISOSEQ          } from '../subworkflows/local/pacbio_isoseq'
include { TRANSCRIPTOME          } from '../subworkflows/local/transcriptome'
include { PREDICTED_PROTEOME     } from '../subworkflows/local/predicted_proteome'
include { PROTEOMICS             } from '../subworkflows/local/proteomics'
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
    ch_samplesheet         // channel: combined samplesheet with both RNA and protein samples

    main:

    ch_versions = channel.empty()

    //
    // Separate RNA and protein samples based on sample_type metadata
    //
    def isRnaSample = { meta ->
        if (meta.containsKey('sample_types')) {
            return meta.sample_types.any { type -> type.toLowerCase() == 'rna' }
        }
        else if (meta.containsKey('sample_type')) {
            return meta.sample_type.toLowerCase() == 'rna'
        }
        else {
            return true
        }
    }

    // Create separate RNA and protein channels using filter
    ch_rna_samples = ch_samplesheet
        .filter { meta, data -> isRnaSample(meta) }

    ch_protein_samples = ch_samplesheet
        .filter { meta, data ->
            meta.containsKey('sample_type') && meta.sample_type.toLowerCase() == 'protein'
        }

    //
    // Count RNA samples to determine if RNA subworkflows should run
    //
    ch_rna_samples
        .toList()
        .map { samples ->
            def count = samples.size()
            if (count > 0) {
                log.info "Detected ${count} RNA sample(s) - RNA analysis subworkflows will run"
            } else {
                log.info "No RNA samples detected - skipping RNA analysis subworkflows (PACBIO_ISOSEQ, TRANSCRIPTOME, PREDICTED_PROTEOME)"
            }
            return tuple(count, samples)
        }
        .set { ch_rna_data }

    ch_rna_count = ch_rna_data.map { count, _samples -> count }
    ch_rna_samples_filtered = ch_rna_data
        .filter { count, _samples -> count > 0 }
        .flatMap { _count, samples -> samples }

    //
    // SUBWORKFLOW: Run PacBio IsoSeq analysis (only if RNA samples present)
    //
    PACBIO_ISOSEQ (
        ch_rna_samples_filtered,
        params.fasta
    )
    ch_versions = ch_versions.mix(PACBIO_ISOSEQ.out.versions.ifEmpty([]))

    //
    // SUBWORKFLOW: Run SQANTI3 QC and filtering (only if RNA samples present)
    //
    // Note: samplesheet is used as sample_metadata
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
    ch_versions = ch_versions.mix(TRANSCRIPTOME.out.versions.ifEmpty([]))

    //
    // SUBWORKFLOW: Run predicted proteome analysis (only if RNA samples present)
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
    ch_versions = ch_versions.mix(PREDICTED_PROTEOME.out.versions.ifEmpty([]))

    //
    // SUBWORKFLOW: Run proteomics analysis (only if protein samples are present)
    //
    // Count protein samples to determine if PROTEOMICS should run
    ch_protein_samples
        .toList()
        .map { samples ->
            def count = samples.size()
            if (count > 0) {
                log.info "Detected ${count} protein sample(s) - PROTEOMICS subworkflow will run"
            } else {
                log.info "No protein samples detected - skipping PROTEOMICS subworkflow"
            }
            return tuple(count, samples)
        }
        .set { ch_protein_data }

    ch_protein_count = ch_protein_data.map { count, _samples -> count }
    ch_protein_samples_filtered = ch_protein_data
        .filter { count, _samples -> count > 0 }
        .flatMap { _count, samples -> samples }

    // Determine protein database source
    // If RNA samples were processed, use predicted proteome; otherwise use provided reference
    ch_rna_count
        .map { count ->
            if (count == 0) {
                // No RNA samples - must use provided protein database
                if (params.gencode_protein_fasta) {
                    log.info "No RNA samples detected - using provided protein database: ${params.gencode_protein_fasta}"
                    return file(params.gencode_protein_fasta)
                } else {
                    error "ERROR: No protein database available. For protein-only analysis, you must specify --gencode_protein_fasta"
                }
            } else {
                return null
            }
        }
        .filter { db -> db != null }
        .mix(
            PREDICTED_PROTEOME.out.protein_fasta
                .map { _meta, fasta -> fasta }
                .ifEmpty([])
        )
        .first()
        .set { ch_protein_db }

    // Prepare MetaMorpheus config file as channel
    ch_metamorpheus_config = channel.value(
        params.metamorpheus_config ?
            file(params.metamorpheus_config) :
            file("${projectDir}/sample_data/SearchTask.toml")
    )
    // Prepare mm_writable directory as a channel
    ch_mm_writable = channel.value(file("${projectDir}/assets/mm_writable_placeholder"))

    PROTEOMICS (
        ch_protein_samples_filtered,
        ch_protein_db,
        ch_metamorpheus_config,
        ch_mm_writable
    )
    ch_versions = ch_versions.mix(PROTEOMICS.out.versions.ifEmpty([]))

    //
    // SUBWORKFLOW: Run differential analysis (optional - requires RNA samples)
    //
    // Only run if run_differential_analysis is enabled AND RNA samples are present
    ch_rna_count
        .subscribe { count ->
            if (params.run_differential_analysis && count > 0) {
                log.info "Differential analysis enabled and RNA samples detected - will run MULTISAMPLE_ANALYSIS"
            } else if (params.run_differential_analysis && count == 0) {
                log.warn "Differential analysis enabled but no RNA samples detected - skipping MULTISAMPLE_ANALYSIS"
            }
        }

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

