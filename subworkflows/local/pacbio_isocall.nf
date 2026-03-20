/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ISOSEQ_ALIGN as ISOCALL_ALIGN } from '../../modules/local/isoseq_align/main'
include { ISOCALL_PROFILE                } from '../../modules/local/isocall_profile/main'
include { ISOCALL_MERGE                  } from '../../modules/local/isocall_merge/main'
include { ISOCALL_PREP                   } from '../../modules/local/isocall_prep/main'
include { ISOCALL_CALL                   } from '../../modules/local/isocall_call/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PACBIO_ISOCALL {
    take:
    ch_samples              // channel: [meta, flnc_bam] - each sample is already individual
    reference_fasta         // path: reference genome FASTA
    reference_gtf_gz        // path: reference genome GTF (gzipped - for ISOCALL_PREP)
    config_toml             // path: isocall configuration TOML file

    main:
    ch_versions = channel.empty()

    //
    // Samples come in as individual channel items, one per RNA sample
    //
    ch_flnc_bams = ch_samples
        .map { meta, bam ->
            [meta, file(bam)]
        }

    //
    // MODULE: Align each FLNC BAM to reference genome
    //
    ISOCALL_ALIGN (
        ch_flnc_bams,
        reference_fasta
    )
    ch_versions = ch_versions.mix(ISOCALL_ALIGN.out.versions)

    //
    // MODULE: Create junction chain profiles for each aligned BAM
    //
    ISOCALL_PROFILE (
        ISOCALL_ALIGN.out.bam
    )
    ch_versions = ch_versions.mix(ISOCALL_PROFILE.out.versions)

    //
    // MODULE: Prepare known isoforms database from GTF (runs independently)
    //
    ISOCALL_PREP (
        reference_gtf_gz
    )
    ch_versions = ch_versions.mix(ISOCALL_PREP.out.versions)

    //
    // Collect all profiles and create a single meta for merged analysis
    //
    ch_all_profiles = ISOCALL_PROFILE.out.profile
        .map { meta, profile -> profile }
        .collect()
        .map { profiles ->
            def meta = [id: params.dataset_name]
            [meta, profiles]
        }

    //
    // MODULE: Merge all sample profiles into single merged profile for input into isocall call 
    //
    ISOCALL_MERGE (
        ch_all_profiles
    )
    ch_versions = ch_versions.mix(ISOCALL_MERGE.out.versions)

    //
    // MODULE: Call isoforms using merged profiles and known isoforms
    //
    ISOCALL_CALL (
        ISOCALL_MERGE.out.merged_profile,
        ISOCALL_PREP.out.isoforms,
        reference_fasta,
        config_toml
    )
    ch_versions = ch_versions.mix(ISOCALL_CALL.out.versions)

    emit:
    // ISOCALL_ALIGN outputs (per sample)
    aligned_bam      = ISOCALL_ALIGN.out.bam                // [meta, *.aligned.bam]
    aligned_bai      = ISOCALL_ALIGN.out.bai                // [meta, *.aligned.bam.bai]

    // ISOCALL_PROFILE outputs (per sample)
    profiles         = ISOCALL_PROFILE.out.profile          // [meta, *_profile.gz]

    // ISOCALL_MERGE outputs
    merged_profile   = ISOCALL_MERGE.out.merged_profile     // [meta, *_merged_profiles.gz]

    // ISOCALL_PREP outputs
    known_isoforms   = ISOCALL_PREP.out.isoforms            // [*.isoforms.gz]

    // ISOCALL_CALL outputs (main outputs for downstream processing)
    called_gtf       = ISOCALL_CALL.out.gtf                 // [meta, *.isoforms.gtf.gz]
    count_matrix     = ISOCALL_CALL.out.count_matrix        // [meta, *.count_matrix.txt]

    versions         = ch_versions.unique().collectFile(name: 'versions.yml')
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
