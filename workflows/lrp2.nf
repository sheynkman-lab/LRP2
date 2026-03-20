/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GUNZIP as GUNZIP_FASTA         } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF           } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENCODE_FASTA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PROTEIN_FASTA } from '../modules/nf-core/gunzip/main'
include { GZIP as GZIP_GTF               } from '../modules/local/gzip/main'
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
    // Decompress reference files if they are gzipped (e.g. GENCODE files gzipped by default)
    //
    def fasta_file = params.fasta ? file(params.fasta) : null
    def is_fasta_gzipped = fasta_file && fasta_file.name.endsWith('.gz')

    if (is_fasta_gzipped) {
        ch_fasta_input = channel.of([[ id: 'genome_fasta' ], fasta_file])
        GUNZIP_FASTA(ch_fasta_input)
        ch_fasta = GUNZIP_FASTA.out.gunzip.map { _meta, file -> file }.first()
    } else {
        ch_fasta = fasta_file ? channel.value(fasta_file) : channel.empty()
    }

    // Handle GTF decompression and create both uncompressed and gzipped versions
    // ch_gtf: uncompressed GTF for TRANSCRIPTOME (SQANTI3)
    // ch_gtf_gz: gzipped GTF for PACBIO_ISOCALL (ISOCALL_PREP)
    def gtf_file = params.gencode_gtf ? file(params.gencode_gtf) : null
    def is_gtf_gzipped = gtf_file && gtf_file.name.endsWith('.gz')

    if (is_gtf_gzipped) {
        // GTF is already gzipped - decompress for TRANSCRIPTOME and keep original for PACBIO_ISOCALL
        ch_gtf_input = channel.of([[ id: 'gencode_gtf' ], gtf_file])
        GUNZIP_GTF(ch_gtf_input)
        ch_gtf = GUNZIP_GTF.out.gunzip.map { _meta, file -> file }.first()
        ch_gtf_gz = channel.value(gtf_file)
    } else {
        // GTF is not gzipped - use as-is for TRANSCRIPTOME and compress for PACBIO_ISOCALL
        ch_gtf = gtf_file ? channel.value(gtf_file) : channel.empty()
        if (gtf_file) {
            ch_gtf_gzip_input = channel.of([[ id: 'gencode_gtf' ], gtf_file])
            GZIP_GTF(ch_gtf_gzip_input)
            ch_gtf_gz = GZIP_GTF.out.gzip.map { _meta, file -> file }.first()
        } else {
            ch_gtf_gz = channel.empty()
        }
    }

    // Handle gencode_fasta decompression (used by TRANSCRIPTOME)
    def gencode_fasta_file = params.gencode_fasta ? file(params.gencode_fasta) : null
    def is_gencode_fasta_gzipped = gencode_fasta_file && gencode_fasta_file.name.endsWith('.gz')

    if (is_gencode_fasta_gzipped) {
        ch_gencode_fasta_input = channel.of([[ id: 'gencode_fasta' ], gencode_fasta_file])
        GUNZIP_GENCODE_FASTA(ch_gencode_fasta_input)
        ch_gencode_fasta = GUNZIP_GENCODE_FASTA.out.gunzip.map { _meta, file -> file }.first()
    } else {
        ch_gencode_fasta = gencode_fasta_file ? channel.value(gencode_fasta_file) : channel.empty()
    }

    //
    // Separate RNA and protein samples based on sample_type metadata
    //
    // Create separate RNA and protein channels using filter
    ch_rna_samples = ch_samplesheet
        .filter { meta, _data ->
            meta.containsKey('sample_type') && meta.sample_type.toLowerCase() == 'rna'
        }

    ch_protein_samples = ch_samplesheet
        .filter { meta, _data ->
            meta.containsKey('sample_type') && meta.sample_type.toLowerCase() == 'protein'
        }

    //
    // Count RNA samples to determine if RNA subworkflows should run
    // Branch the channel into two: one for counting, one for processing
    //
    ch_rna_samples
        .tap { ch_rna_for_count }
        .set { ch_rna_samples_filtered }

    // Count the samples
    ch_rna_for_count
        .toList()
        .map { samples ->
            def count = samples.size()
            if (count > 0) {
                log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.green} Detected ${count} RNA sample(s) - RNA analysis subworkflows will run${colors.reset}-"
            } else {
                log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.dim} No RNA samples detected - skipping RNA analysis subworkflows (PACBIO_ISOCALL, TRANSCRIPTOME, PREDICTED_PROTEOME)${colors.reset}-"
            }
            return count
        }
        .set { ch_rna_count }

    //
    // SUBWORKFLOW: Run PacBio IsoCall analysis (only if RNA samples present)
    //
    // IsoCall requires config TOML and gzipped GTF reference
    ch_isocall_config = channel.value(file("${projectDir}/bin/isocall_config.toml"))

    PACBIO_ISOCALL (
        ch_rna_samples_filtered,
        ch_fasta,
        ch_gtf_gz,
        ch_isocall_config
    )
    ch_versions = ch_versions.mix(PACBIO_ISOCALL.out.versions.ifEmpty([]))

    //
    // SUBWORKFLOW: Run SQANTI3 QC and filtering (only if RNA samples present)
    //
    // Note: samplesheet is used as sample_metadata
    sample_metadata_file = params.sample_metadata ?: params.input

    TRANSCRIPTOME (
        PACBIO_ISOCALL.out.called_gtf
            .join(PACBIO_ISOCALL.out.count_matrix, by: 0)
            .map { meta, gtf, count ->
                [meta, gtf, count] },
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
    def protein_fasta_path = params.protein_fasta

    // If not provided by user, get from GENCODE genome default
    if (!protein_fasta_path && params.gencode_refs?.containsKey(params.genome)) {
        protein_fasta_path = params.gencode_refs[params.genome].protein_fasta
        log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.cyan} Auto-detected protein FASTA from GENCODE genome ${params.genome}: ${protein_fasta_path}${colors.reset}-"
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
            if (protein_fasta_path && count == 0) {
                log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} Protein FASTA available but no protein samples detected - skipping PROTEOMICS subworkflow${colors.reset}-"
            } else if (!protein_fasta_path && count > 0) {
                log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} Protein samples detected but no protein FASTA available (neither --protein_fasta provided nor auto-detected from GENCODE genome) - skipping PROTEOMICS subworkflow${colors.reset}-"
            } else if (!protein_fasta_path && count == 0) {
                log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.dim} No protein samples detected - skipping PROTEOMICS subworkflow${colors.reset}-"
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

        // BUILD_PROTEOME_REFERENCE search db creation logic:
        // - If RNA samples were processed, we build sample-specific references with LRP proteome + GENCODE concatenated
        // - If no RNA samples then we build GENCODE-only references per sample group

        // Extract outputs from RNA subworkflows if available (keep the RNA sample meta.id for use in filtering by CPM column name)
        ch_predicted_proteome_fasta = PREDICTED_PROTEOME.out.protein_all_orfs_fasta
            .map { _meta, fasta -> fasta }
            .first()
            .ifEmpty(file('NO_FILE'))

        ch_transcript_counts_with_id = TRANSCRIPTOME.out.hashids_all
            .map { meta, counts -> [meta.id, counts] }
            .first()
            .ifEmpty(['NO_RNA_SAMPLE', file('NO_FILE')])

        ch_transcript_counts = ch_transcript_counts_with_id
            .map { rna_id, counts -> counts }

        ch_rna_sample_id = ch_transcript_counts_with_id
            .map { rna_id, counts -> rna_id }

        // Extract CDS GTF and ORF FASTA for novel peptides classification
        ch_lr_cds_gtf = PREDICTED_PROTEOME.out.protein_cds_gtf_copy
            .map { _meta, gtf -> gtf }
            .first()
            .ifEmpty(file('NO_FILE'))

        ch_lr_orf_fasta = PREDICTED_PROTEOME.out.protein_all_orfs_fasta
            .map { _meta, fasta -> fasta }
            .first()
            .ifEmpty(file('NO_FILE'))

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

        // Build proteome references per sample group
        // Prepare input channel for BUILD_PROTEOME_REFERENCE
        ch_build_ref_input = ch_protein_samples_grouped
            .map { meta, _files ->
                [meta, meta]  // Just pass meta for now
            }
            .combine(ch_predicted_proteome_fasta)
            .combine(ch_transcript_counts)
            .map { meta, _meta2, lrp_fasta, counts ->
                return [meta, lrp_fasta, counts]
            }

        // Script path for build_mass_spec_reference.R
        ch_build_proteome_script = channel.value(file("${projectDir}/bin/build_mass_spec_reference.R"))

        // BUILD_PROTEOME_REFERENCE runs once for each sample group to create per-group sample-specific references
        // Pass params.genome directly for file naming (e.g., GRCh38.p14.v49)
        BUILD_PROTEOME_REFERENCE(
            ch_build_ref_input,
            ch_build_proteome_script,
            params.genome,
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
            ch_lr_cds_gtf,
            ch_lr_orf_fasta
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
                log.info "-${colors.purple}[sheynkmanlab/lrp2]${colors.green} Differential analysis enabled and RNA samples detected - will run MULTISAMPLE_ANALYSIS${colors.reset}-"
            } else if (params.run_differential_analysis && count == 0) {
                log.warn "-${colors.purple}[sheynkmanlab/lrp2]${colors.yellow} Differential analysis enabled but no RNA samples detected - skipping MULTISAMPLE_ANALYSIS${colors.reset}-"
            }
        }

    if (params.run_differential_analysis) {
         // Prepare transcript channel with GTF and counts
         ch_transcripts = TRANSCRIPTOME.out.corrected_gtf_filtered
             .join(TRANSCRIPTOME.out.hashids_filtered, by: 0)

         // Prepare ORF counts channel
         ch_orfs = PREDICTED_PROTEOME.out.hashids_orf

         // Create RNA-only metadata file by filtering the original samplesheet
         def metadata_content = file(sample_metadata_file).text
         def lines = metadata_content.split('\n')
         def header_parts = lines[0].split(',')
         def sample_name_idx = header_parts.findIndexOf { it.trim() == 'sample_name' }
         def sample_path_idx = header_parts.findIndexOf { it.trim() == 'sample_path' }
         def condition_idx = header_parts.findIndexOf { it.trim() == 'condition' }
         def sample_type_idx = header_parts.findIndexOf { it.trim() == 'sample_type' }

         // Filter to RNA samples and create new CSV with 'name' and 'group' columns
         def rna_metadata_lines = ['name,sample_path,group,sample_type']
         lines.drop(1).each { line ->
             if (line.trim()) {
                 def parts = line.split(',')
                 if (parts.size() > sample_type_idx && parts[sample_type_idx].trim().toLowerCase() == 'rna') {
                     def name = parts[sample_name_idx].trim()
                     def path = parts[sample_path_idx].trim()
                     def group = parts[condition_idx].trim()
                     rna_metadata_lines << "${name},${path},${group},RNA"
                 }
             }
         }
         def rna_metadata_file = file("${workDir}/rna_samples_metadata.csv")
         rna_metadata_file.text = rna_metadata_lines.join('\n') + '\n'

         MULTISAMPLE_ANALYSIS (
             ch_transcripts,
             ch_orfs,
             rna_metadata_file,
             file(params.lr_leafcutter_script),
             file(params.leafcutter_ds_script),
             file(params.multisample_script),
             params.control_group,
             params.experimental_group,
             params.min_samples_per_intron,
             params.min_samples_per_group,
             params.min_usage_ratio,
             params.drimseq_min_gene_expr,
             params.drimseq_min_isoform_prop
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
