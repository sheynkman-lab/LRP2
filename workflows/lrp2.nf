/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP as GUNZIP_FASTA         } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF           } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENCODE_FASTA } from '../modules/nf-core/gunzip/main'
include { PACBIO_ISOSEQ                  } from '../subworkflows/local/pacbio_isoseq'
include { TRANSCRIPTOME                  } from '../subworkflows/local/transcriptome'
include { PREDICTED_PROTEOME             } from '../subworkflows/local/predicted_proteome'
include { PROTEOMICS                     } from '../subworkflows/local/proteomics'
include { MULTISAMPLE_ANALYSIS           } from '../subworkflows/local/multisample_analysis'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText         } from '../subworkflows/local/utils_nfcore_lrp2_pipeline'

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
    // Decompress reference files if they are gzipped (e.g. GENCODE files gzipped by default)
    //
    def fasta_file = params.fasta ? file(params.fasta) : null
    def is_fasta_gzipped = fasta_file && fasta_file.name.endsWith('.gz')

    if (is_fasta_gzipped) {
        ch_fasta_input = channel.of([[ id: 'genome_fasta' ], fasta_file])
        GUNZIP_FASTA(ch_fasta_input)
        ch_fasta = GUNZIP_FASTA.out.gunzip.map { _meta, file -> file }
    } else {
        ch_fasta = fasta_file ? channel.value(fasta_file) : channel.empty()
    }

    // Handle GTF decompression
    def gtf_file = params.gencode_gtf ? file(params.gencode_gtf) : null
    def is_gtf_gzipped = gtf_file && gtf_file.name.endsWith('.gz')

    if (is_gtf_gzipped) {
        ch_gtf_input = channel.of([[ id: 'gencode_gtf' ], gtf_file])
        GUNZIP_GTF(ch_gtf_input)
        ch_gtf = GUNZIP_GTF.out.gunzip.map { _meta, file -> file }
    } else {
        ch_gtf = gtf_file ? channel.value(gtf_file) : channel.empty()
    }

    // Handle gencode_fasta decompression (used by TRANSCRIPTOME)
    def gencode_fasta_file = params.gencode_fasta ? file(params.gencode_fasta) : null
    def is_gencode_fasta_gzipped = gencode_fasta_file && gencode_fasta_file.name.endsWith('.gz')

    if (is_gencode_fasta_gzipped) {
        ch_gencode_fasta_input = channel.of([[ id: 'gencode_fasta' ], gencode_fasta_file])
        GUNZIP_GENCODE_FASTA(ch_gencode_fasta_input)
        ch_gencode_fasta = GUNZIP_GENCODE_FASTA.out.gunzip.map { _meta, file -> file }
    } else {
        ch_gencode_fasta = gencode_fasta_file ? channel.value(gencode_fasta_file) : channel.empty()
    }

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
        ch_fasta
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
        ch_gtf,
        ch_gencode_fasta,
        file(sample_metadata_file),
        file(params.filter_script),
        file(params.hashlib_script),
        file(params.generate_hashids_script)
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
        ch_gtf,
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

    //
    // SUBWORKFLOW: Run proteomics analysis
    // Only runs if --gencode_protein_fasta parameter is provided
    // This ensures users explicitly opt-in to proteomics analysis
    //
    ch_protein_count
        .subscribe { count ->
            if (params.gencode_protein_fasta && count > 0) {
                log.info "Detected ${count} protein sample(s) and --gencode_protein_fasta provided - PROTEOMICS subworkflow will run"
            } else if (params.gencode_protein_fasta && count == 0) {
                log.warn "--gencode_protein_fasta provided but no protein samples detected - skipping PROTEOMICS subworkflow"
            } else if (!params.gencode_protein_fasta && count > 0) {
                log.warn "Protein samples detected but --gencode_protein_fasta not provided - skipping PROTEOMICS subworkflow"
            } else {
                log.info "No protein samples detected - skipping PROTEOMICS subworkflow"
            }
        }

    // note: pipeline will only execute PROTEOMICS if gencode_protein_fasta parameter is provided by user!
    if (params.gencode_protein_fasta) {
        // If RNA samples were processed, use predicted proteome; otherwise use provided reference
        ch_protein_db = ch_rna_count
            .map { count ->
                if (count == 0) {
                    // No RNA samples - must use provided protein database
                    log.info "No RNA samples detected - using provided protein database: ${params.gencode_protein_fasta}"
                    return file(params.gencode_protein_fasta)
                } else {
                    // Signal to use predicted proteome from RNA analysis
                    return null
                }
            }
            .filter { db -> db != null }
            .mix(
                PREDICTED_PROTEOME.out.protein_all_orfs_fasta
                    .map { _meta, fasta -> fasta }
            )
            .first()

        ch_metamorpheus_config = channel.value(
            params.metamorpheus_config ?
                file(params.metamorpheus_config) :
                file("${projectDir}/sample_data/SearchTask.toml")
        )
        ch_mm_writable = channel.value(file("${projectDir}/assets/mm_writable_placeholder"))

        PROTEOMICS (
            ch_protein_samples_filtered,
            ch_protein_db,
            ch_metamorpheus_config,
            ch_mm_writable
        )
        ch_versions = ch_versions.mix(PROTEOMICS.out.versions.ifEmpty([]))
    }

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

