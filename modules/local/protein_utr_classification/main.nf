process PROTEIN_UTR_CLASSIFICATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/jtllab/lrp2-lite:latest"

    input:
    tuple val(meta), path(protein_classification), path(cds_gtf), path(corrected_fasta), path(hashids_cpm), path(all_orfs_mapped)
    path reference_gtf
    path protein_class_script

    output:
    tuple val(meta), path("protein_sqanti/*_protein_classification_all_isoforms.txt"), emit: protein_all_isoforms
    tuple val(meta), path("protein_sqanti/*.predicted.proteome.all_best_orfs.fa"), emit: protein_all_orfs_fasta
    tuple val(meta), path("protein_sqanti/*.predicted.proteome.high_confidence_ORF_cpm.txt"), emit: hashids_orf
    tuple val(meta), path("protein_sqanti/*.predicted.proteome.high_confidence_ORF.gtf"), emit: protein_gtf
    tuple val(meta), path("protein_sqanti/*.predicted.proteome.high_confidence_ORF.bed"), emit: protein_bed
    tuple val(meta), path("protein_sqanti/*.sqanti_protein_classification.tsv"), emit: protein_classification_copy
    tuple val(meta), path("orf_calling/*_corrected_filtered_CDS.gtf"), emit: cds_gtf_copy
    tuple val(meta), path("orf_calling/*_all_orfs_mapped.tsv"), emit: all_orfs_mapped_copy
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_junctions_after_stop = task.ext.min_junctions_after_stop ?: params.min_junctions_after_stop_codon
    def protein_class_keep = task.ext.protein_class_keep ?: params.protein_class_keep

    """
    export OUTPUT_DIR=.
    export OUTPUT_BASE_NAME=$prefix
    export GENCODE_GTF_FILE=\$(pwd)/$reference_gtf
    export MIN_JUNCTIONS_AFTER_STOP_CODON=$min_junctions_after_stop
    export PROTEIN_CLASS_KEEP="$protein_class_keep"
    export R_LIBS_USER=""
    export R_LIBS="/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"

    mkdir -p orf_calling protein_sqanti sqanti_transcript

    # Place files in expected locations with expected names
    ln -sf \$(pwd)/$cds_gtf orf_calling/${prefix}_corrected_filtered_CDS.gtf
    ln -sf \$(pwd)/$corrected_fasta sqanti_transcript/${prefix}_corrected_filtered.fasta
    ln -sf \$(pwd)/$protein_classification protein_sqanti/${prefix}.sqanti_protein_classification.tsv
    ln -sf \$(pwd)/$hashids_cpm sqanti_transcript/${prefix}_hashids_with_cpm_filtered.txt
    ln -sf \$(pwd)/$all_orfs_mapped orf_calling/${prefix}_all_orfs_mapped.tsv

    Rscript \$(pwd)/$protein_class_script $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/.*R version //g' | sed 's/ .*//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p protein_sqanti orf_calling
    touch protein_sqanti/${prefix}_protein_classification_all_isoforms.txt
    touch protein_sqanti/${prefix}.predicted.proteome.all_best_orfs.fa
    touch protein_sqanti/${prefix}.predicted.proteome.high_confidence_ORF_cpm.txt
    touch protein_sqanti/${prefix}.predicted.proteome.high_confidence_ORF.gtf
    touch protein_sqanti/${prefix}.predicted.proteome.high_confidence_ORF.bed
    touch protein_sqanti/${prefix}.sqanti_protein_classification.tsv
    touch orf_calling/${prefix}_corrected_filtered_CDS.gtf
    touch orf_calling/${prefix}_all_orfs_mapped.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.3.0
    END_VERSIONS
    """
}
