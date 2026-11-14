/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CHECK_BAM_ALIGNMENT            } from '../../modules/local/check_bam_alignment/main'
include { PBMM2_ALIGN as ISOCALL_ALIGN   } from '../../modules/local/pbmm2_align/main'
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
    // Note: this is split by file type: FASTQ files vs BAM files.
    // For BAM files, store the original path in meta for later .bai file lookup
    //
    ch_samples
        .map { meta, sample_path ->
            def sample_file = file(sample_path)
            def meta_with_original = meta + [original_bam_path: sample_path.toString()]
            [meta_with_original, sample_file]
        }
        .branch { _meta, sample_file ->
            fastq: sample_file.toString() =~ /\.(fastq|fastq\.gz|fq|fq\.gz)$/
            bam: sample_file.toString() =~ /\.bam$/
        }
        .set { ch_input_split }

    //
    // PROCESS: Check all BAM files to determine if they're already pbmm2-aligned by checking PG header for pbmm2
    //
    CHECK_BAM_ALIGNMENT (
        ch_input_split.bam
    )
    ch_versions = ch_versions.mix(CHECK_BAM_ALIGNMENT.out.versions)

    //
    // Split BAMs into two channels based on alignment status
    //
    CHECK_BAM_ALIGNMENT.out.check
        .branch { _meta, _bam, is_aligned ->
            aligned: is_aligned.trim() == "true"
            unaligned: is_aligned.trim() == "false"
        }
        .set { ch_bam_split }

    //
    // Combine all files that need alignment: FASTQ files + unaligned BAMs
    //
    ch_needs_alignment = ch_input_split.fastq
        .mix(ch_bam_split.unaligned.map { meta, bam, _is_aligned -> [meta, bam] })

    //
    // MODULE: Align FASTQ files and unaligned BAMs to reference genome
    //
    ISOCALL_ALIGN (
        ch_needs_alignment,
        reference_fasta
    )
    ch_versions = ch_versions.mix(ISOCALL_ALIGN.out.versions)

    //
    // Extract pre-aligned BAMs (alignment check already logged in CHECK_BAM_ALIGNMENT)
    //
    ch_prealigned = ch_bam_split.aligned
        .map { meta, bam, _is_aligned -> [meta, bam] }

    //
    // Combine pre-aligned BAMs with newly aligned files
    //
    ch_all_aligned_bams = ch_prealigned
        .mix(ISOCALL_ALIGN.out.bam)

    //
    // MODULE: Create junction chain profiles for each aligned BAM
    //
    ISOCALL_PROFILE (
        ch_all_aligned_bams
    )
    ch_versions = ch_versions.mix(ISOCALL_PROFILE.out.versions)

    //
    // MODULE: Prepare known isoforms database from GTF (only runs if samples are present)
    //
    ISOCALL_PREP (
        reference_gtf_gz
    )
    ch_versions = ch_versions.mix(ISOCALL_PREP.out.versions)

    // Sample IDs in samplesheet order (ch_samples is emitted in samplesheet order)
    ch_sample_order = ch_samples
        .map { meta, _sample -> meta.id }
        .collect()

    // Collect profiles as they finish, then reorder them to match the samplesheet
    ch_all_profiles = ISOCALL_PROFILE.out.profile
        .map { meta, profile -> [meta.id, profile] }
        .collect(flat: false)
        .map { pairs -> [pairs] }
        .combine(ch_sample_order.map { ids -> [ids] })
        .map { pairs, ids ->
            def profile_by_id = pairs.collectEntries { id, profile -> [(id): profile] }
            [[id: params.dataset_name], ids.collect { id -> profile_by_id[id] }]
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
    // Combined aligned BAMs (both pre-aligned and newly aligned)
    aligned_bam      = ch_all_aligned_bams                  // [meta, *.bam]
    aligned_bai      = ch_prealigned
        .map { meta, bam ->
            // For pre-aligned BAMs, use the original path to find .bai file
            def bai_path = meta.original_bam_path ? file("${meta.original_bam_path}.bai") : file("${bam}.bai")
            [meta, bai_path]
        }
        .mix(ISOCALL_ALIGN.out.bai)                         // [meta, *.bam.bai]

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
