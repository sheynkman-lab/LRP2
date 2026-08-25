/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP as GUNZIP_GTF           } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_FASTA         } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PROTEIN_FASTA } from '../modules/nf-core/gunzip/main'
include { GZIP as GZIP_GTF               } from '../modules/local/gzip/main'
include { SANITIZE_PAR_IDS               } from '../modules/local/sanitize_par_ids/main'
include { BUILD_PROTEOME_REFERENCE       } from '../modules/local/build_proteome_reference/main'
include { PACBIO_ISOCALL                 } from '../subworkflows/local/pacbio_isocall'
include { TRANSCRIPTOME                  } from '../subworkflows/local/transcriptome'
include { PREDICTED_PROTEOME             } from '../subworkflows/local/predicted_proteome'
include { PROTEOMICS                     } from '../subworkflows/local/proteomics'
include { MULTISAMPLE_ANALYSIS           } from '../subworkflows/local/multisample_analysis'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { logColours                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
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

    // Setup color codes for logging
    def colors = logColours(params.monochrome_logs)

    //
    // Detect multisample-only mode
    //
    def is_multisample_only = params.transcripts_gtf &&
                               params.transcript_counts &&
                               params.orf_counts &&
                               params.multisample_metadata

    //
    // Detect skip-isocall mode (for non-PacBio data like ONT)
    //
    def skip_isocall_mode = params.external_gtf && params.external_counts

    //
    // Separate RNA and protein samples based on sample_type metadata
    // Note: Genome FASTA decompression is deferred until after we determine if RNA samples are present, and skip samplesheet parsing if in multisample-only mode
    //
    def samplesheet_file = !is_multisample_only ? file(params.input) : null
    def has_rna_samples_sync = false
    def has_protein_samples_sync = false

    if (!is_multisample_only) {
        // Parse samplesheet synchronously to determine if RNA samples exist (for conditional subworkflow execution)
        def samplesheet_content = samplesheet_file.text
        def samplesheet_lines = samplesheet_content.split('\n')
        def samplesheet_header_parts = samplesheet_lines[0].split(',')
        def sample_type_idx = samplesheet_header_parts.findIndexOf { it.trim() == 'sample_type' }

        samplesheet_lines.drop(1).each { line ->
            if (line.trim()) {
                def parts = line.split(',')
                if (parts.size() > sample_type_idx) {
                    def sample_type = parts[sample_type_idx].trim().toLowerCase()
                    if (sample_type == 'rna') {
                        has_rna_samples_sync = true
                    } else if (sample_type == 'protein') {
                        has_protein_samples_sync = true
                    }
                }
            }
        }
    }

    // Create separate RNA and protein channels using filter
    // Skip channel creation if in multisample-only mode (no samplesheet to parse)
    ch_rna_samples = !is_multisample_only ?
        ch_samplesheet.filter { meta, _data ->
            meta.containsKey('sample_type') && meta.sample_type.toLowerCase() == 'rna'
        } : channel.empty()

    ch_protein_samples = !is_multisample_only ?
        ch_samplesheet.filter { meta, _data ->
            meta.containsKey('sample_type') && meta.sample_type.toLowerCase() == 'protein'
        } : channel.empty()

    //
    // Count RNA samples to determine if RNA subworkflows should run
    // note: we branch the channel into three portions currently: one for counting, one for sample names, one for processing
    //
    ch_rna_samples
        .tap { ch_rna_for_count }
        .tap { ch_rna_for_names }
        .set { ch_rna_samples_filtered }

    // RNA sample names for matching with protein samples
    ch_rna_sample_names = ch_rna_for_names
        .map { meta, _data -> meta.sample_name }
        .toList()

    ch_rna_for_count
        .toList()
        .map { samples ->
            def count = samples.size()
            // Only log RNA sample detection in normal mode (not multisample-only mode)
            if (!is_multisample_only) {
                if (count > 0) {
                    if (skip_isocall_mode) {
                        log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.green} Detected ${count} RNA sample(s) - RNA analysis will start from S2 TRANSCRIPTOME (skipping S1 PACBIO ISOCALL)${colors.reset}-"
                    } else {
                        log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.green} Detected ${count} RNA sample(s) - RNA analysis subworkflows will run${colors.reset}-"
                    }
                } else {
                    log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.dim} No RNA samples detected - skipping RNA analysis subworkflows (S1 PACBIO ISOCALL, S2 TRANSCRIPTOME, S3 PREDICTED PROTEOME)${colors.reset}-"
                }
            }
            return count
        }
        .set { ch_rna_count }

    //
    // Handle GTF and FASTA decompression ONLY if RNA samples are present AND not in multisample-only mode
    //
    def gtf_file = !is_multisample_only && params.gtf ? file(params.gtf) : null
    def is_gtf_gzipped = gtf_file && gtf_file.name.endsWith('.gz')
    def fasta_file = !is_multisample_only && params.fasta ? file(params.fasta) : null
    def is_fasta_gzipped = fasta_file && fasta_file.name.endsWith('.gz')

    // Only process GTF and FASTA files if not in multisample-only mode and RNA samples are present
    if (!is_multisample_only) {
        ch_rna_count
            .map { count -> count > 0 }
            .set { ch_has_rna_samples }

        ch_gtf_input_conditional = ch_has_rna_samples
            .filter { has_rna -> has_rna && gtf_file != null }
            .map { _has_rna -> [[ id: 'gtf' ], gtf_file] }

        if (is_gtf_gzipped && gtf_file) {
            // GTF is already gzipped - decompress for TRANSCRIPTOME and keep original for PACBIO_ISOCALL
            GUNZIP_GTF(ch_gtf_input_conditional)

            // Conditionally sanitize PAR gene IDs: known _PAR_Y pattern exists for transcripts in GENCODE v25-43, but user could also provide their own GTF with this
            GUNZIP_GTF.out.gunzip
                .branch { meta, gtf ->
                    needs_sanitization: {
                        // Use exception-based early exit to avoid loading entire GTF into memory (versions after 43 get big)
                        try {
                            gtf.eachLine { line ->
                                if (line.contains('_PAR_Y')) {
                                    throw new RuntimeException('PAR_Y_FOUND')
                                }
                            }
                            false
                        } catch (RuntimeException e) {
                            if (e.message == 'PAR_Y_FOUND') {
                                true  
                            } else {
                                throw e  
                            }
                        }
                    }()
                    no_sanitization: true
                }
                .set { ch_gtf_branched }

            // Run sanitization only on files that need it
            SANITIZE_PAR_IDS(ch_gtf_branched.needs_sanitization)
            ch_gtf = ch_gtf_branched.no_sanitization
                .mix(SANITIZE_PAR_IDS.out.sanitized_gtf)
                .map { _meta, file -> file }

            // For gzipped GTF, create channel from filtered input (will be empty if no RNA samples)
            ch_gtf_gz = ch_gtf_input_conditional.map { _meta, file -> file }
        } else if (gtf_file) {
            // GTF is not gzipped - use as-is for TRANSCRIPTOME and compress for PACBIO_ISOCALL
            ch_gtf_for_check = ch_has_rna_samples
                .filter { has_rna -> has_rna }
                .map { _has_rna -> [[ id: 'gtf' ], gtf_file] }
                .ifEmpty([[ id: 'gtf' ], gtf_file])

            // Conditionally sanitize PAR gene IDs: known _PAR_Y pattern exists for transcripts in GENCODE v25-43, but user could also provide their own GTF with this
            ch_gtf_for_check
                .branch { meta, gtf ->
                    needs_sanitization: gtf.text.contains('_PAR_Y')
                    no_sanitization: true
                }
                .set { ch_gtf_branched }

            // Run sanitization only on files that need it
            SANITIZE_PAR_IDS(ch_gtf_branched.needs_sanitization)
            ch_gtf = ch_gtf_branched.no_sanitization
                .mix(SANITIZE_PAR_IDS.out.sanitized_gtf)
                .map { _meta, file -> file }

            GZIP_GTF(ch_gtf_input_conditional)
            ch_gtf_gz = GZIP_GTF.out.gzip.map { _meta, file -> file }
        } else {
            // No GTF file provided
            ch_gtf = channel.empty()
            ch_gtf_gz = channel.empty()
        }

        ch_fasta_input_conditional = ch_has_rna_samples
            .filter { has_rna -> has_rna && is_fasta_gzipped && fasta_file != null }
            .map { _has_rna -> [[ id: 'fasta' ], fasta_file] }

        if (is_fasta_gzipped) {
            GUNZIP_FASTA(ch_fasta_input_conditional)
            ch_fasta = GUNZIP_FASTA.out.gunzip.map { _meta, file -> file }
        } else {
            ch_fasta = ch_has_rna_samples
                .filter { has_rna -> has_rna && fasta_file != null }
                .map { _has_rna -> fasta_file }
                .ifEmpty(fasta_file ? channel.value(fasta_file) : channel.empty())
        }
    } else {
        // In multisample-only mode, skip all GTF/FASTA processing
        ch_gtf = channel.empty()
        ch_gtf_gz = channel.empty()
        ch_fasta = channel.empty()
    }

    //
    // Define sample metadata file (used by both RNA subworkflows and multisample analysis)
    //
    sample_metadata_file = params.sample_metadata ?: params.input

    //
    // RNA-specific subworkflows (only execute if RNA samples are present AND not in multisample-only mode)
    //
    if (has_rna_samples_sync && !is_multisample_only) {

        // Declare channels for GTF and counts that will be used by TRANSCRIPTOME
        // These will be populated either from PACBIO_ISOCALL or from external files
        def ch_called_gtf
        def ch_count_matrix

        if (!skip_isocall_mode) {
            //
            // SUBWORKFLOW: Run PacBio IsoCall analysis (normal mode)
            //
            // IsoCall requires config TOML and gzipped GTF reference for ISOCALL_PREP
            ch_isocall_config = channel.value(file("${projectDir}/bin/isocall_config.toml"))

            PACBIO_ISOCALL (
                ch_rna_samples_filtered,
                ch_fasta,
                ch_gtf_gz,
                ch_isocall_config
            )
            ch_versions = ch_versions.mix(PACBIO_ISOCALL.out.versions.ifEmpty([]))
            ch_called_gtf = PACBIO_ISOCALL.out.called_gtf
            ch_count_matrix = PACBIO_ISOCALL.out.count_matrix

        } else {
            //
            // Skip-IsoCall mode: Use external GTF and counts from user
            //
            log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Using external GTF and counts - skipping PACBIO_ISOCALL subworkflow${colors.reset}-"

            // Create channels from external files
            ch_called_gtf = channel.of([
                [id: params.dataset_name],
                file(params.external_gtf)
            ])

            ch_count_matrix = channel.of([
                [id: params.dataset_name],
                file(params.external_counts)
            ])
        }

        //
        // SUBWORKFLOW: Run SQANTI3 QC and filtering
        //
        TRANSCRIPTOME (
            ch_called_gtf
                .join(ch_count_matrix, by: 0)
                .map { meta, gtf, count ->
                    [meta, gtf, count] },
            ch_gtf,
            ch_fasta,
            file(sample_metadata_file),
            file(params.filter_script),
            file(params.hashlib_script),
            file(params.generate_hashids_script)
        )
        ch_versions = ch_versions.mix(TRANSCRIPTOME.out.versions.ifEmpty([]))

        //
        // SUBWORKFLOW: Run predicted proteome analysis
        //
        // Log species information
        if (params.genome && params.gencode_refs?.containsKey(params.genome)) {
            log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Auto-detected species from genome ${params.genome}: ${params.species}${colors.reset}-"
        } else if (params.species) {
            log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Species: ${params.species}${colors.reset}-"
        }

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
                .join(TRANSCRIPTOME.out.hashids_mapping, by: 0)
                .map { meta, fasta, gtf, classification, hashids, mapping ->
                    [meta, fasta, gtf, classification, hashids, mapping] },
            ch_gtf,
            hexamer_file,
            logit_model,
            file(params.filter_cpat_script),
            file(params.sqanti_protein_script),
            file(params.protein_class_script)
        )
        ch_versions = ch_versions.mix(PREDICTED_PROTEOME.out.versions.ifEmpty([]))
    }

    //
    // SUBWORKFLOW: Run proteomics analysis (only if protein samples are present)
    //
    // Count protein samples to determine if PROTEOMICS should run
    // Branch the channel into two: one for counting, one for processing
    ch_protein_samples
        .tap { ch_protein_for_count }
        .set { ch_protein_samples_filtered }

    // Count the samples
    ch_protein_for_count
        .toList()
        .map { samples ->
            def count = samples.size()
            if (count > 0) {
                log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.green} Detected ${count} protein sample(s) - PROTEOMICS subworkflow will run${colors.reset}-"
            } else {
                log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.dim} No protein samples detected - skipping PROTEOMICS subworkflow${colors.reset}-"
            }
            return count
        }
        .set { ch_protein_count }

    //
    // Determine protein FASTA for proteomics analysis
    // Priority: 1) User-provided --protein_fasta, 2) Auto-detect from GENCODE genome, 3) Skip if neither available
    //
    //def protein_fasta_path = params.protein_fasta

    // If not provided by user, get from GENCODE genome default
    //if (!protein_fasta_path && params.gencode_refs?.containsKey(params.genome)) {
    //    protein_fasta_path = params.gencode_refs[params.genome].protein_fasta
    //    log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Auto-detected protein FASTA from GENCODE genome ${params.genome}: ${protein_fasta_path}${colors.reset}-"
    //}
    
    // Resolve GENCODE protein FASTA from --genome (e.g., GRCh38.p14.v49)
    def gencode_protein_fasta_path = null
    if (params.gencode_refs?.containsKey(params.genome)) {
      gencode_protein_fasta_path = params.gencode_refs[params.genome].protein_fasta
      log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} GENCODE protein FASTA resolved from --genome ${params.genome}: ${gencode_protein_fasta_path}${colors.reset}-"
    }
    
    // Resolve protein FASTA for the proteomics search database (optional, user-provided)
    // Accepts LRP2 output or any custom FASTA; source is auto-detected in build_mass_spec_reference.R
    def s5_custom_protein_fasta_path = params.S5_custom_protein_fasta ?: null
    if (s5_custom_protein_fasta_path) {
      log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Protein FASTA provided for search database: ${s5_custom_protein_fasta_path}${colors.reset}-"
    }
    
    // Resolve CDS GTF paired with the protein FASTA (for novel peptide BED mapping)
    def s5_custom_cds_gtf_file = params.S5_custom_cds_gtf ? file(params.S5_custom_cds_gtf) : null
    if (s5_custom_cds_gtf_file) {
        log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Custom CDS GTF provided: ${s5_custom_cds_gtf_file}${colors.reset}-"
    }
    
    def s5_custom_counts_path = params.S5_custom_counts ?: null
    if (s5_custom_counts_path) {
      log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Count matrix provided for search database filtering: ${s5_custom_counts_path}${colors.reset}-"
    }
    //def protein_fasta_file = protein_fasta_path ? file(protein_fasta_path) : null
    //def is_protein_fasta_gzipped = protein_fasta_file && protein_fasta_file.name.endsWith('.gz')

    //if (is_protein_fasta_gzipped) {
    //    ch_protein_fasta_input = channel.of([[ id: 'protein_fasta' ], protein_fasta_file])
    //    GUNZIP_PROTEIN_FASTA(ch_protein_fasta_input)
    //    ch_protein_fasta = GUNZIP_PROTEIN_FASTA.out.gunzip.map { _meta, file -> file }
    //} else {
    //    ch_protein_fasta = protein_fasta_file ? channel.value(protein_fasta_file) : channel.empty()
    //}

    ch_has_protein_samples = ch_protein_count.map { count -> count > 0 }

    //
    // SUBWORKFLOW: Run proteomics analysis
    // Runs if protein_fasta is available (either user-provided or auto-detected) AND protein samples are present
    //
    //ch_protein_count
    //    .subscribe { count ->
    //        if (protein_fasta_path && count == 0) {
    //            log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} Protein FASTA available but no protein samples detected - skipping PROTEOMICS subworkflow${colors.reset}-"
    //        } else if (!protein_fasta_path && count > 0) {
    //            log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} Protein samples detected but no protein FASTA available (neither --protein_fasta provided nor auto-detected from GENCODE genome) - skipping PROTEOMICS subworkflow${colors.reset}-"
    //        } else if (!protein_fasta_path && count == 0) {
    //            log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.dim} No protein samples detected - skipping PROTEOMICS subworkflow${colors.reset}-"
    //        }
    //    }

    // note: pipeline will only execute PROTEOMICS if protein_fasta is available (user-provided or auto-detected) AND protein samples exist
    //if (protein_fasta_path) {
    if ((gencode_protein_fasta_path || s5_custom_protein_fasta_path) && has_protein_samples_sync) {
        ch_metamorpheus_config = channel.value(
            params.metamorpheus_config ?
                file(params.metamorpheus_config) :
                file("${projectDir}/sample_data/SearchTask.toml")
        )
        ch_mm_writable = channel.value(file("${projectDir}/assets/mm_writable_placeholder"))
        
        //
        // Decompress GENCODE protein FASTA if provided and gzipped (only when protein samples exist)
        //
        def gencode_protein_fasta_file = gencode_protein_fasta_path ? file(gencode_protein_fasta_path) : null
        def is_gencode_protein_fasta_gzipped = gencode_protein_fasta_file && gencode_protein_fasta_file.name.endsWith('.gz')

        if (is_gencode_protein_fasta_gzipped && gencode_protein_fasta_file) {
            ch_gencode_protein_fasta_input = ch_has_protein_samples
                .filter { it }
                .map { _has_protein -> [[ id: 'gencode_protein_fasta' ], gencode_protein_fasta_file] }
            GUNZIP_PROTEIN_FASTA(ch_gencode_protein_fasta_input)
            ch_gencode_protein_fasta = GUNZIP_PROTEIN_FASTA.out.gunzip.map { _meta, file -> file }
        } else if (gencode_protein_fasta_file) {
            ch_gencode_protein_fasta = channel.value(gencode_protein_fasta_file)
        } else {
            ch_gencode_protein_fasta = channel.value(file('NO_FILE'))
        }

        def s5_custom_protein_fasta_file = s5_custom_protein_fasta_path ? file(s5_custom_protein_fasta_path) : null
        def s5_custom_counts_file = s5_custom_counts_path ? file(s5_custom_counts_path) : null
        
        // BUILD_PROTEOME_REFERENCE search db creation logic:
        // - If RNA samples were processed, we build sample-specific references with LRP proteome + GENCODE concatenated
        // - If no RNA samples then we build GENCODE-only references per sample group

        // Extract outputs from RNA subworkflows if available (keep the RNA sample meta.id for use in filtering by CPM column name)
        if (has_rna_samples_sync) {
            ch_predicted_proteome_fasta = PREDICTED_PROTEOME.out.protein_all_orfs_fasta
                .map { _meta, fasta -> fasta }
                .first()
                .ifEmpty(s5_custom_protein_fasta_file ?: file('NO_FILE'))

            ch_transcript_counts_with_id = TRANSCRIPTOME.out.hashids_all
                .map { meta, counts -> [meta.id, counts] }
                .first()
                .ifEmpty(['NO_RNA_SAMPLE', file('NO_FILE')])

            ch_transcript_counts = ch_transcript_counts_with_id
                .map { rna_id, counts -> counts }

            // NOVEL_PEPTIDES Extract CDS GTF and ORF FASTA
            ch_lr_cds_gtf = PREDICTED_PROTEOME.out.cds_gtf
                .map { _meta, gtf -> gtf }
                .first()
                .ifEmpty(s5_custom_cds_gtf_file ?: file('NO_FILE'))

            ch_lr_orf_fasta = PREDICTED_PROTEOME.out.protein_all_orfs_fasta
                .map { _meta, fasta -> fasta }
                .first()
                .ifEmpty(s5_custom_protein_fasta_file ?: file('NO_FILE'))
        } else {
            // No RNA samples - use placeholder/user-provided files for proteomics-only mode
            ch_predicted_proteome_fasta = channel.value(s5_custom_protein_fasta_file ?: file('NO_FILE'))
            ch_transcript_counts_with_id = channel.value(['NO_RNA_SAMPLE', file('NO_FILE')])
            ch_transcript_counts = channel.value(s5_custom_counts_file ?: file('NO_FILE'))
            ch_lr_cds_gtf = channel.value(s5_custom_cds_gtf_file ?: file('NO_FILE'))
            ch_lr_orf_fasta = channel.value(s5_custom_protein_fasta_file ?: file('NO_FILE'))
        }

        // Resolve GENCODE annotation GTF for novel peptides BED mapping
        // Use the already-decompressed ch_gtf if available, otherwise resolve from params
        def gtf_for_novel = params.gtf ? file(params.gtf) : null
        ch_gtf_for_novel = gtf_for_novel
            ? channel.value(gtf_for_novel)
            : channel.value(file('NO_FILE'))
            
        // Create a channel that maps each protein sample to its sample_name for grouping
        // Group protein samples by sample_name (the biosample group)
        ch_protein_samples_grouped = ch_protein_samples_filtered
            .map { meta, file ->
                def sample_name = meta.sample_name
                return [sample_name, meta, file]
            }
            .groupTuple(by: 0)
            .map { sample_name, metas, files ->
                // Use the first meta but update id to be the sample_name for the grouped reference
                def grouped_meta = metas[0] + [id: sample_name]
                return [grouped_meta, files]
            }

        ch_build_ref_input = ch_protein_samples_grouped
            .combine(ch_predicted_proteome_fasta)
            .combine(ch_transcript_counts)
            .combine(ch_gencode_protein_fasta)
            .map { meta, _files, custom_fasta, counts, gencode_protein_fasta ->
                // Create unique placeholder names per sample to avoid Nextflow staging collisions
                def unique_custom = custom_fasta.name == 'NO_FILE' ? file("${meta.id}_NO_CUSTOM_FASTA") : custom_fasta
                def unique_counts = counts.name == 'NO_FILE' ? file("${meta.id}_NO_COUNTS") : counts
                def unique_gencode = gencode_protein_fasta.name == 'NO_FILE' ? file("${meta.id}_NO_GENCODE_PROTEIN_FASTA") : gencode_protein_fasta
                return [meta, unique_counts, unique_custom, unique_gencode]
            }
            
        // Script path for build_mass_spec_reference.R
        ch_build_proteome_script = channel.value(file("${projectDir}/bin/build_mass_spec_reference.R"))

        // BUILD_PROTEOME_REFERENCE runs once for each sample group to create per-group sample-specific references
        BUILD_PROTEOME_REFERENCE(
            ch_build_ref_input,
            ch_build_proteome_script,
            params.no_gencode ?: false
        )
        ch_versions = ch_versions.mix(BUILD_PROTEOME_REFERENCE.out.versions.ifEmpty([]))

        // Join the built references back to the protein samples
        // Map the reference output by sample_name and join with protein samples
        ch_protein_db_by_sample = BUILD_PROTEOME_REFERENCE.out.reference_fasta
            .map { meta, fasta -> [meta.id, fasta] }

        ch_protein_for_proteomics = ch_protein_samples_filtered
            .map { meta, file -> [meta.sample_name, meta, file] }
            .combine(ch_protein_db_by_sample, by: 0)
            .map { _sample_name, meta, file, protein_db -> [meta, file, protein_db] }

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
            params.fragpipe_license_accept,
            // Novel peptides inputs
            ch_rna_sample_names,    // List of RNA sample names to check for matches
            ch_lr_cds_gtf,          // Custom/LRP CDS GTF (for samples with matched RNA)
            ch_lr_orf_fasta,        // Custom/LRP ORF FASTA (for samples with matched RNA)
            ch_gtf_for_novel,   // Reference annotation GTF (for BED mapping)
            ch_gencode_protein_fasta   // GENCODE protein FASTA (for BED mapping)
        )
        ch_versions = ch_versions.mix(PROTEOMICS.out.versions.ifEmpty([]))
    }

    //
    // SUBWORKFLOW: Run differential analysis (automatic - based on RNA samples and conditions)
    //
    // Conditions we check to determine if MULTISAMPLE_ANALYSIS should run:
    // 1. RNA samples are present
    // 2. At least 2 unique conditions exist
    // 3. Each condition that runs has at least 2 samples (warn re. statistical robustness but still run if <3 per condition)
    //
    // Initialize should_run_multisample
    def should_run_multisample = false
    def lines
    def header_parts
    def sample_name_idx
    def condition_idx
    def multisample_sample_type_idx
    def samples_per_condition

    // Only parse metadata for validation if not in multisample-only mode
    if (!is_multisample_only) {
        // Parse metadata synchronously to count samples per condition
        def metadata_content = file(sample_metadata_file).text
        lines = metadata_content.split('\n')
        header_parts = lines[0].split(',')
        sample_name_idx = header_parts.findIndexOf { it.trim() == 'sample_name' }
        condition_idx = header_parts.findIndexOf { it.trim() == 'condition' }
        multisample_sample_type_idx = header_parts.findIndexOf { it.trim() == 'sample_type' }
        samples_per_condition = [:].withDefault { [] }
        lines.drop(1).each { line ->
            if (line.trim()) {
                def parts = line.split(',')
                if (parts.size() > multisample_sample_type_idx && parts[multisample_sample_type_idx].trim().toLowerCase() == 'rna') {
                    def condition = parts[condition_idx].trim()
                    def sample_name = parts[sample_name_idx].trim()
                    samples_per_condition[condition] << sample_name
                }
            }
        }
        def unique_conditions = samples_per_condition.keySet().size()
        def conditions_with_min_samples = samples_per_condition.findAll { k, v -> v.size() >= 2 }
        def conditions_with_robust_samples = samples_per_condition.findAll { k, v -> v.size() >= 3 }

        // Determine if we should run multisample analysis based on samples present and log final call
        should_run_multisample = unique_conditions >= 2 && conditions_with_min_samples.size() >= 2

        ch_rna_count.subscribe { count ->
            if (count == 0) {
                log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.dim} No RNA samples detected - skipping MULTISAMPLE_ANALYSIS${colors.reset}-"
            } else if (samples_per_condition.size() > 0) {
                if (unique_conditions < 2) {
                    log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} Only ${unique_conditions} condition detected - need at least 2 conditions for MULTISAMPLE_ANALYSIS${colors.reset}-"
                } else if (conditions_with_min_samples.size() < 2) {
                    log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} Only ${conditions_with_min_samples.size()} condition(s) with ≥2 samples - need at least 2 conditions with ≥2 samples each for MULTISAMPLE_ANALYSIS${colors.reset}-"
                } else {
                    // log pairwise comparisons
                    def valid_conditions_for_log = samples_per_condition.findAll { condition, samples -> samples.size() >= 2 }.keySet()
                    def conditions_list_for_log = valid_conditions_for_log.toList().sort()
                    def comparison_count = 0
                    def comparisons_to_log = []
                    conditions_list_for_log.eachWithIndex { cond1, i ->
                        conditions_list_for_log.eachWithIndex { cond2, j ->
                            if (j > i) {
                                comparison_count += 1
                                comparisons_to_log << [cond1, cond2]
                            }
                        }
                    }

                    log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Running MULTISAMPLE_ANALYSIS for ${comparison_count} pairwise comparison(s):${colors.reset}-"
                    comparisons_to_log.each { pair ->
                        log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan}   ${pair[0]} (${samples_per_condition[pair[0]].size()} samples) vs ${pair[1]} (${samples_per_condition[pair[1]].size()} samples)${colors.reset}-"
                    }

                    // Show warning if some conditions have < 3 samples
                    if (conditions_with_robust_samples.size() < conditions_with_min_samples.size()) {
                        log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} WARNING: Some conditions have <3 samples per condition. Please interpret results with caution as they may not be statistically robust, and consider adding more replicates if possible."
                    }
                }
            }
        }
    }

    // Prepare channels for MULTISAMPLE_ANALYSIS
    def ch_transcripts
    def ch_orfs
    def rna_metadata_file

    if (is_multisample_only) {
        // Multisample-only mode: create channels from provided files
        log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Creating channels from pre-computed files for multisample-only analysis${colors.reset}-"

        ch_transcripts = channel.of([
            [id: params.dataset_name],
            file(params.transcripts_gtf),
            file(params.transcript_counts)
        ])

        ch_orfs = channel.of([
            [id: params.dataset_name],
            file(params.orf_counts)
        ])

        rna_metadata_file = file(params.multisample_metadata)

    } else if (should_run_multisample) {
         // Normal mode: prepare channels from RNA subworkflow outputs
         ch_transcripts = TRANSCRIPTOME.out.corrected_gtf_filtered
             .join(TRANSCRIPTOME.out.hashids_filtered, by: 0)

         // Prepare ORF counts channel
         ch_orfs = PREDICTED_PROTEOME.out.hashids_orf

         // Create RNA-only metadata file by filtering the original samplesheet
         rna_metadata_file = file("${workDir}/${params.dataset_name}_rna_samples_metadata.csv")
         def sample_path_idx = header_parts.findIndexOf { it.trim() == 'sample_path' }

         // Filter to RNA samples and create new CSV with 'name' and 'group' columns
         def rna_metadata_lines = ['name,sample_path,group,sample_type']
         lines.drop(1).each { line ->
             if (line.trim()) {
                 def parts = line.split(',')
                 if (parts.size() > multisample_sample_type_idx && parts[multisample_sample_type_idx].trim().toLowerCase() == 'rna') {
                     def name = parts[sample_name_idx].trim()
                     def path = parts[sample_path_idx].trim()
                     def group = parts[condition_idx].trim()
                     rna_metadata_lines << "${name},${path},${group},RNA"
                 }
             }
         }
         rna_metadata_file.text = rna_metadata_lines.join('\n') + '\n'
    }

    if (is_multisample_only || should_run_multisample) {

         MULTISAMPLE_ANALYSIS (
             ch_transcripts,
             ch_orfs,
             rna_metadata_file,
             file(params.lr_leafcutter_script),
             file(params.leafcutter_ds_script),
             file(params.multisample_script),
             params.min_samples_per_intron,
             params.min_samples_per_group,
             params.min_usage_ratio,
             params.drimseq_min_gene_expr,
             params.drimseq_min_isoform_prop,
             params.dataset_name,
             params.leafcutter_threads
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
