/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ISOSEQ_CLUSTER         } from '../modules/nf-core/isoseq/cluster/main'
include { PACBIO_ISOSEQ          } from '../subworkflows/local/pacbio_isoseq'
include { SQANTI_TRANSCRIPT      } from '../subworkflows/local/sqanti_transcript'
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
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Run PacBio IsoSeq analysis (merge, cluster, align, collapse)
    //
    PACBIO_ISOSEQ (
        ch_samplesheet,
        params.fasta
    )
    ch_versions = ch_versions.mix(PACBIO_ISOSEQ.out.versions)

    //
    // SUBWORKFLOW: Run SQANTI3 QC and filtering
    //
    // Note: samplesheet is used as sample_metadata if not explicitly provided
    sample_metadata_file = params.sample_metadata ?: params.input

    SQANTI_TRANSCRIPT (
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
    ch_versions = ch_versions.mix(SQANTI_TRANSCRIPT.out.versions)

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

