/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP as GUNZIP_FASTA         } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF           } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENCODE_FASTA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PROTEIN_FASTA } from '../modules/nf-core/gunzip/main'
include { CAT_CAT                        } from '../modules/nf-core/cat/cat/main'
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
    // Create separate RNA and protein channels using filter
    ch_rna_samples = ch_samplesheet
        .filter { meta, data ->
            meta.containsKey('sample_type') && meta.sample_type.toLowerCase() == 'rna'
        }

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
    // Determine protein FASTA for proteomics analysis
    // Priority: 1) User-provided --protein_fasta, 2) Auto-detect from GENCODE genome, 3) Skip if neither available
    //
    def protein_fasta_path = params.protein_fasta

    // If not provided by user, get from GENCODE genome default
    if (!protein_fasta_path && params.gencode_refs?.containsKey(params.genome)) {
        protein_fasta_path = params.gencode_refs[params.genome].protein_fasta
        log.info "Auto-detected protein FASTA from GENCODE genome ${params.genome}: ${protein_fasta_path}"
    }

    def protein_fasta_file = protein_fasta_path ? file(protein_fasta_path) : null
    def is_protein_fasta_gzipped = protein_fasta_file && protein_fasta_file.name.endsWith('.gz')

    if (is_protein_fasta_gzipped) {
        ch_protein_fasta_input = channel.of([[ id: 'protein_fasta' ], protein_fasta_file])
        GUNZIP_PROTEIN_FASTA(ch_protein_fasta_input)
        ch_protein_fasta = GUNZIP_PROTEIN_FASTA.out.gunzip.map { _meta, file -> file }
    } else {
        ch_protein_fasta = protein_fasta_file ? channel.value(protein_fasta_file) : channel.empty()
    }

    //
    // SUBWORKFLOW: Run proteomics analysis
    // Runs if protein_fasta is available (either user-provided or auto-detected) AND protein samples are present
    //
    ch_protein_count
        .subscribe { count ->
            if (protein_fasta_path && count > 0) {
                log.info "Detected ${count} protein sample(s) and protein FASTA available - PROTEOMICS subworkflow will run"
            } else if (protein_fasta_path && count == 0) {
                log.warn "Protein FASTA available but no protein samples detected - skipping PROTEOMICS subworkflow"
            } else if (!protein_fasta_path && count > 0) {
                log.warn "Protein samples detected but no protein FASTA available (neither --protein_fasta provided nor auto-detected from GENCODE genome) - skipping PROTEOMICS subworkflow"
            } else {
                log.info "No protein samples detected - skipping PROTEOMICS subworkflow"
            }
        }

    // note: pipeline will only execute PROTEOMICS if protein_fasta is available (user-provided or auto-detected)
    if (protein_fasta_path) {
        ch_metamorpheus_config = channel.value(
            params.metamorpheus_config ?
                file(params.metamorpheus_config) :
                file("${projectDir}/sample_data/SearchTask.toml")
        )
        ch_mm_writable = channel.value(file("${projectDir}/assets/mm_writable_placeholder"))

        // Prepare protein samples - check if matching RNA samples exist:
        // 1) For RNA+protein samplesheets with matching sample names, concatenate predicted proteome with reference
        // 2) For protein-only samples, just use reference FASTA

        // Create a channel from predicted proteome FASTAs keyed by sample_name
        // This will be empty if no RNA samples were processed
        ch_predicted_proteomes = PREDICTED_PROTEOME.out.protein_all_orfs_fasta
            .map { meta, fasta -> [meta.sample_name, fasta] }
            .ifEmpty([])

        // Join protein samples with predicted proteomes by sample_name
        // This will only match samples that have both RNA and protein data
        ch_protein_with_predicted = ch_protein_samples_filtered
            .map { meta, file -> [meta.sample_name, meta, file] }
            .join(ch_predicted_proteomes, by: 0, remainder: true)
            .branch { _sample_name, meta, file, predicted_fasta ->
                has_rna: predicted_fasta != null
                    return [meta, file, predicted_fasta]
                no_rna: true
                    return [meta, file]
            }

        // For samples WITH matching RNA data: concatenate predicted proteome + reference FASTA
        ch_concat_input = ch_protein_with_predicted.has_rna
            .combine(ch_protein_fasta)
            .map { meta, _file, predicted_fasta, reference_fasta ->
                // Create input for CAT_CAT: [meta, [file1, file2]]
                def cat_meta = meta + [id: "${meta.sample_name}_combined_proteome"]
                log.info "Protein sample ${meta.id} has matching RNA sample - concatenating predicted proteome with reference FASTA"
                return [cat_meta, [reference_fasta, predicted_fasta]]
            }

        // Concatenate FASTAs for samples with RNA data (only if there are any)
        CAT_CAT(ch_concat_input)

        // Combine concatenated FASTAs back with original protein sample metadata
        ch_protein_with_concat = ch_protein_with_predicted.has_rna
            .map { meta, file, _predicted_fasta -> [meta.sample_name, meta, file] }
            .join(
                CAT_CAT.out.file_out.map { meta, fasta -> [meta.id.replaceAll('_combined_proteome$', ''), fasta] },
                by: 0
            )
            .map { _sample_name, meta, file, concatenated_fasta ->
                [meta, file, concatenated_fasta]
            }

        // For samples WITHOUT matching RNA data: use reference FASTA
        ch_protein_no_rna = ch_protein_with_predicted.no_rna
            .combine(ch_protein_fasta)
            .map { meta, file, protein_fasta ->
                log.info "Protein sample ${meta.id} has no matching RNA sample - using reference protein FASTA"
                return [meta, file, protein_fasta]
            }

        // Combine both branches
        ch_protein_for_proteomics = ch_protein_with_concat.mix(ch_protein_no_rna)

        // Pass protein samples with their corresponding protein databases to PROTEOMICS
        PROTEOMICS (
            ch_protein_for_proteomics,  // [meta, file, protein_db]
            ch_metamorpheus_config,
            ch_mm_writable,
            params.protein_search,
            params.fragpipe_workflow,
            // FragPipe authentication parameters
            params.fragpipe_first_name,
            params.fragpipe_last_name,
            params.fragpipe_email,
            params.fragpipe_institution,
            params.fragpipe_token,
            params.fragpipe_license_accept
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

