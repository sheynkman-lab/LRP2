process SQANTI_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker://docker.io/anaconesalab/sqanti3:5.2.2'

    input:
    tuple val(meta), path(isoforms_gtf), path(flnc_count)
    path reference_gtf
    path reference_fasta

    output:
    tuple val(meta), path("*_classification.txt"), emit: classification
    tuple val(meta), path("*_corrected.gtf"), emit: corrected_gtf
    tuple val(meta), path("*_corrected.fasta"), emit: corrected_fasta
    tuple val(meta), path("*_corrected.genePred"), emit: corrected_genepred
    tuple val(meta), path("*_corrected.gtf.cds.gff"), emit: corrected_cds_gff
    tuple val(meta), path("*.isoforms.gtf"), emit: isoforms_gtf
    tuple val(meta), path("*_junctions.txt"), emit: junctions
    tuple val(meta), path("*.params.txt"), emit: params
    tuple val(meta), path("refAnnotation_*.genePred"), emit: refannotation_genepred
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    source /conda/miniconda3/etc/profile.d/conda.sh
    conda activate SQANTI3.env

    # Decompress GTF if it's gzipped
    ISOFORMS_INPUT="$isoforms_gtf"
    if [[ "$isoforms_gtf" == *.gtf.gz ]]; then
        gunzip -c "$isoforms_gtf" > "${prefix}.isoforms.gtf"
        ISOFORMS_INPUT="${prefix}.isoforms.gtf"
    # Convert GFF to GTF if necessary (isoseq collapse outputs GFF but SQANTI3 requires GTF)
    elif [[ "$isoforms_gtf" == *.gff ]]; then
        gffread "$isoforms_gtf" -T -o "${prefix}.isoforms.gtf"
        ISOFORMS_INPUT="${prefix}.isoforms.gtf"
    fi

    sqanti3_qc.py \\
        --force_id_ignore \\
        --skipORF \\
        --output $prefix \\
        --dir . \\
        --cpus $task.cpus \\
        --report skip \\
        --fl_count $flnc_count \\
        "\$ISOFORMS_INPUT" \\
        $reference_gtf \\
        $reference_fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.2.2
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_classification.txt
    touch ${prefix}_corrected.gtf
    touch ${prefix}_corrected.fasta
    touch ${prefix}_corrected.genePred
    touch ${prefix}_corrected.gtf.cds.gff
    touch ${prefix}.isoforms.gtf
    touch ${prefix}_junctions.txt
    touch ${prefix}.params.txt
    touch refAnnotation_${prefix}.genePred

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.2.2
    END_VERSIONS
    """
}
