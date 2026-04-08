/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MSCONVERT_MZML         } from '../../modules/local/msconvert_mzml/main'
include { METAMORPHEUS           } from '../../modules/local/metamorpheus/main'
include { FRAGPIPE               } from '../../modules/local/fragpipe/main'
include { FRAGPIPE_AUTHENTICATE  } from '../../modules/local/fragpipe_authenticate/main'
include { NOVEL_PEPTIDES         } from '../../modules/local/novel_peptides/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PROTEOMICS {
    take:
    ch_ms_files_with_db         // channel: [meta, raw_or_mzml_file, protein_fasta]
    metamorpheus_config         // path: MetaMorpheus TOML configuration file
    mm_writable                 // path: writable directory for MetaMorpheus
    protein_search              // val: protein search engine ('fragpipe' or 'metamorpheus')
    fragpipe_workflow           // path: FragPipe workflow file (optional, null triggers auto-download)

    // FragPipe authentication parameters (only needed if protein_search='fragpipe')
    fragpipe_first_name         // val: User's first name
    fragpipe_last_name          // val: User's last name
    fragpipe_email              // val: User's email
    fragpipe_institution        // val: User's institution
    fragpipe_token              // val: 6-digit verification token
    fragpipe_license_accept     // val: License acceptance confirmation

    // Novel peptides inputs
    rna_sample_ids              // channel: List of RNA sample IDs (to check for matched RNA data)
    lr_cds_gtf                  // path: Long-read CDS GTF (only for samples with matched RNA sample)
    lr_orf_fasta                // path: Long-read ORF FASTA (only for samples with matched RNA sample)
    genome                      // val: Genome reference name (e.g., 'GRCh38.p14.v49')

    main:
    ch_versions = channel.empty()

    // Extract [meta, file] and protein_db from input channel
    ch_ms_files = ch_ms_files_with_db.map { meta, file, db -> [meta, file] }
    // Map protein DBs by sample_name and deduplicate (multiple files per sample have same DB)
    ch_protein_dbs = ch_ms_files_with_db
        .map { meta, file, db -> [meta.sample_name, db] }
        .unique()

    //
    // MODULE: Convert any protein sample files to .mzML format if needed
    //
    // Separate files that are already mzML from those that need conversion
    ch_ms_files
        .branch { meta, file ->
            mzml: file.name.endsWith('.mzML') || file.name.endsWith('.mzml')
                return [meta, file]
            needs_conversion: true
                return [meta, file]
        }
        .set { ch_ms_branched }

    MSCONVERT_MZML(
        ch_ms_branched.needs_conversion
    )
    ch_versions = ch_versions.mix(MSCONVERT_MZML.out.versions)

    // Combine converted and already-mzML files
    ch_mzml_files = MSCONVERT_MZML.out.mzml
        .mix(ch_ms_branched.mzml)

    //
    // fix: group mzML files by sample_name, not meta (different per file since different filename)
    //
    ch_mzml_grouped = ch_mzml_files
        .map { meta, file -> [meta.sample_name, meta, file] }
        .groupTuple(by: 0)
        .map { sample_name, metas, files ->
            def grouped_meta = metas[0]
            return [grouped_meta, files]
        }

    //
    // Filter samples based on protein search engine and mass_spec_type compatibility
    //
    ch_mzml_grouped
        .branch { meta, files ->
            // For FragPipe: accept both DDA and DIA
            fragpipe: protein_search == 'fragpipe' && (meta.mass_spec_type == 'DDA' || meta.mass_spec_type == 'DIA')
                return [meta, files]
            // For MetaMorpheus: only accept DDA (cannot handle DIA)
            metamorpheus: protein_search == 'metamorpheus' && meta.mass_spec_type == 'DDA'
                return [meta, files]
            // logging: samples that cannot be processed
            incompatible: true
                if (protein_search == 'metamorpheus' && meta.mass_spec_type == 'DIA') {
                    log.warn "INCOMPATIBLE: Sample ${meta.id} has mass_spec_type='DIA' but protein_search='metamorpheus'. MetaMorpheus cannot process DIA data. Please use --protein_search fragpipe for DIA samples. Skipping sample."
                } else {
                    log.warn "INCOMPATIBLE: Sample ${meta.id} with mass_spec_type='${meta.mass_spec_type}' cannot be processed with protein_search='${protein_search}'. Skipping sample."
                }
                return null
        }
        .set { ch_mzml_branched }

    //
    // MODULE: Run protein database search based on selected engine
    // Join mzML files with their corresponding protein databases (by sample_name)
    //
    ch_mzml_with_db = ch_mzml_branched.fragpipe
        .mix(ch_mzml_branched.metamorpheus)
        .map { meta, files -> [meta.sample_name, meta, files] }
        .join(ch_protein_dbs, by: 0)
        .map { sample_name, meta, files, protein_db -> [meta, files, protein_db] }

    if (protein_search == 'fragpipe') {
        //
        // FRAGPIPE WORKFLOW
        //
        // Step 1: Authenticate and download tools (runs once, cached for future runs)
        // Convert null values to empty strings for Nextflow compatibility
        def fp_first_name = fragpipe_first_name ?: ''
        def fp_last_name = fragpipe_last_name ?: ''
        def fp_email = fragpipe_email ?: ''
        def fp_institution = fragpipe_institution ?: ''
        def fp_token = fragpipe_token ?: ''
        def fp_license = fragpipe_license_accept ?: false

        FRAGPIPE_AUTHENTICATE(
            fp_first_name,
            fp_last_name,
            fp_email,
            fp_institution,
            fp_token,
            fp_license
        )
        ch_versions = ch_versions.mix(FRAGPIPE_AUTHENTICATE.out.versions)

        // Step 2: Group samples by sample_name for FragPipe
        // All files with the same sample_name run together in one FragPipe process
        ch_mzml_with_db
            .filter { meta, files, protein_db -> meta.mass_spec_type == 'DDA' || meta.mass_spec_type == 'DIA' }
            .map { meta, files, protein_db ->
                def sample_name = meta.sample_name
                return [sample_name, meta, files, protein_db]
            }
            .groupTuple(by: 0)
            .flatMap { sample_name, metas, files_list, protein_dbs ->
                // Run FragPipe once for each unique sample_name
                // Use the first meta but update id to be the sample_name for the grouped run
                def grouped_meta = metas[0] + [id: sample_name]
                def protein_db = protein_dbs[0]  // All samples with same sample_name should have same protein_db
                def all_files = files_list.flatten()
                return [[grouped_meta, all_files, protein_db]]
            }
            .set { ch_fragpipe_grouped }

        // Step 3: Run FragPipe
        // Prepare custom workflow file path or use null for auto-download
        def custom_workflow_file = fragpipe_workflow ? file(fragpipe_workflow) : file('NO_FILE')

        FRAGPIPE(
            ch_fragpipe_grouped,
            FRAGPIPE_AUTHENTICATE.out.msfragger_jar,
            FRAGPIPE_AUTHENTICATE.out.ionquant_jar,
            FRAGPIPE_AUTHENTICATE.out.diatracer_jar,
            custom_workflow_file
        )
        ch_versions = ch_versions.mix(FRAGPIPE.out.versions)
        ch_psm_table = FRAGPIPE.out.psm_table
        ch_search_results = FRAGPIPE.out.results

        // Step 4: Run NOVEL_PEPTIDES on FragPipe results
        // Extract GENCODE version from genome parameter (e.g., 'GRCh38.p14.v49' -> '49')
        def gencode_version = genome ? genome.tokenize('.')[-1].replaceAll(/[^0-9]/, '') : '49'

        // Merge DIA and DDA peptide files: since at least one will be present we make sure we get one of the actual files from FragPipe work directory, as fix for novel peptides stalling for larger datasets
        ch_fragpipe_peptide_files = FRAGPIPE.out.peptide_tsv
            .mix(FRAGPIPE.out.combined_peptide_tsv)
            .map { meta, peptide_file -> [meta.sample_name, meta, peptide_file] }

        // Join per-sample protein FASTA with FragPipe peptide files, then add LR files conditionally based on whether sample has matched RNA
        ch_novel_peptides_input = ch_fragpipe_peptide_files
            .combine(ch_protein_dbs, by: 0)
            .map { sample_name, meta, peptide_file, protein_fasta ->
                def rna_ids = rna_sample_ids.val
                def lr_gtf_file = lr_cds_gtf.val
                def lr_fasta_file = lr_orf_fasta.val
                def has_rna = rna_ids.contains(sample_name)
                def lr_gtf = has_rna ? lr_gtf_file : file("${sample_name}_NO_GTF")
                def lr_fasta = has_rna ? lr_fasta_file : file("${sample_name}_NO_FASTA")
                return [meta, peptide_file, protein_fasta, lr_gtf, lr_fasta]
            }

        def novel_peptides_script = file("${projectDir}/bin/novel_peptides.R")

        NOVEL_PEPTIDES(
            ch_novel_peptides_input,
            novel_peptides_script,
            gencode_version
        )
        ch_versions = ch_versions.mix(NOVEL_PEPTIDES.out.versions)
        ch_novel_peptides = NOVEL_PEPTIDES.out.novel_peptides
        ch_peptides_bed = NOVEL_PEPTIDES.out.peptides_bed
    }
    else {
        //
        // METAMORPHEUS WORKFLOW
        //
        // Group samples by sample_name
        ch_mzml_with_db
            .filter { meta, files, protein_db -> meta.mass_spec_type == 'DDA' }
            .map { meta, files, protein_db ->
                def sample_name = meta.sample_name
                return [sample_name, meta, files, protein_db]
            }
            .groupTuple(by: 0)
            .flatMap { sample_name, metas, files_list, protein_dbs ->
                // Run MetaMorpheus once for each unique sample_name
                // Use the first meta but update id to be the sample_name for the grouped run
                def grouped_meta = metas[0] + [id: sample_name]
                def protein_db = protein_dbs[0]  // All samples with same sample_name should have same protein_db
                def all_files = files_list.flatten()
                return [[grouped_meta, all_files, protein_db]]
            }
            .set { ch_metamorpheus_grouped }

        METAMORPHEUS(
            ch_metamorpheus_grouped.map { meta, files, db -> [meta, files] },
            ch_metamorpheus_grouped.map { meta, files, db -> db },
            metamorpheus_config,
            mm_writable
        )
        ch_versions = ch_versions.mix(METAMORPHEUS.out.versions)
        ch_psm_table = METAMORPHEUS.out.psm_table
        ch_search_results = METAMORPHEUS.out.results

        // NOTE: NOVEL_PEPTIDES only runs for FragPipe currently, not MetaMorpheus so we create empty channels to maintain consistent output structure
        // TO DO: incorporate novel peptides classification for MetaMorpheus results in the future and update this workflow accordingly
        ch_novel_peptides = Channel.empty()
        ch_peptides_bed = Channel.empty()
    }

    emit:
    // MSCONVERT outputs
    mzml_files                  = ch_mzml_files                     // [meta, *.mzML]

    // Protein search outputs
    psm_table                   = ch_psm_table                      // [meta, *.psmtsv or *.tsv]
    search_results              = ch_search_results                 // [meta, results/*]

    // Novel peptides outputs
    novel_peptides              = ch_novel_peptides                 // [meta, *_novel_peptides.tsv]
    peptides_bed                = ch_peptides_bed                   // [meta, *_all_peptides.bed] (optional)

    versions                    = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
