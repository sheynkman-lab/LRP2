/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MSCONVERT_MZML } from '../../modules/local/msconvert_mzml/main'
include { METAMORPHEUS   } from '../../modules/local/metamorpheus/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PROTEOMICS {
    take:
    ch_ms_files                 // channel: [meta, raw_or_mzml_file]
    protein_fasta               // path: protein database FASTA file (shared across all samples)
    metamorpheus_config         // path: MetaMorpheus TOML configuration file
    mm_writable                 // path: writable directory for MetaMorpheus

    main:
    ch_versions = channel.empty()

    //
    // MODULE: Convert any protein sample .raw files to .mzML format if needed
    //
    // Separate raw files from mzML files
    ch_ms_files
        .branch { meta, file ->
            raw: file.name.endsWith('.raw')
                return [meta, file]
            mzml: file.name.endsWith('.mzML') || file.name.endsWith('.mzml')
                return [meta, file]
        }
        .set { ch_ms_branched }

    MSCONVERT_MZML(
        ch_ms_branched.raw
    )
    ch_versions = ch_versions.mix(MSCONVERT_MZML.out.versions)

    // Combine converted and already-mzML files
    ch_mzml_files = MSCONVERT_MZML.out.mzml
        .mix(ch_ms_branched.mzml)

    //
    // Group mzML files by sample (meta.id) for multi-file samples
    //
    ch_mzml_grouped = ch_mzml_files
        .groupTuple(by: 0)
        .map { meta, files ->
            // If single file, unwrap from list
            def file_input = files.size() == 1 ? files[0] : files
            [meta, file_input]
        }

    //
    // MODULE: Run MetaMorpheus database search, with all samples searched against same protein database
    //
    METAMORPHEUS(
        ch_mzml_grouped,
        protein_fasta,
        metamorpheus_config,
        mm_writable
    )
    ch_versions = ch_versions.mix(METAMORPHEUS.out.versions)

    emit:
    // MSCONVERT outputs
    mzml_files                  = ch_mzml_files                     // [meta, *.mzML]

    // METAMORPHEUS outputs
    psm_table                   = METAMORPHEUS.out.psm_table        // [meta, *.psmtsv]
    search_results              = METAMORPHEUS.out.results          // [meta, results/*]

    versions                    = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
