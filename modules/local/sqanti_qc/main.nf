process SQANTI_QC {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container 'docker://docker.io/anaconesalab/sqanti3:5.2.2'

    input:
    tuple val(meta), path(isoforms_gtf), path(flnc_count)
    path reference_gtf
    path reference_fasta

    output:
    tuple val(meta), path("*.transcriptome.SQANTI_classification.txt"), emit: classification
    tuple val(meta), path("*.transcriptome.gtf"), emit: corrected_gtf
    tuple val(meta), path("*.transcriptome.fasta"), emit: corrected_fasta
//    tuple val(meta), path("*.transcriptome.corrected.genePred"), emit: corrected_genepred
//    tuple val(meta), path("*.transcriptome.corrected.gtf.cds.gff"), emit: corrected_cds_gff
//    tuple val(meta), path("*.transcriptome.isoforms.gtf"), emit: isoforms_gtf
    tuple val(meta), path("*.transcriptome.junctions.txt"), emit: junctions
//    tuple val(meta), path("*.transcriptome.params.txt"), emit: params
//    tuple val(meta), path("refAnnotation.*.genePred"), emit: refannotation_genepred
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
        gunzip -c "$isoforms_gtf" > "${prefix}.isocall.isoforms.gtf"
        ISOFORMS_INPUT="${prefix}.isocall.isoforms.gtf"
    # Convert GFF to GTF if necessary (isoseq collapse outputs GFF but SQANTI3 requires GTF)
    elif [[ "$isoforms_gtf" == *.gff ]]; then
        gffread "$isoforms_gtf" -T -o "${prefix}.isocall.isoforms.gtf"
        ISOFORMS_INPUT="${prefix}.isocall.isoforms.gtf"
    fi

    sqanti3_qc.py \\
        --force_id_ignore \\
        --skipORF \\
        --output ${prefix}.transcriptome \\
        --dir . \\
        --cpus $task.cpus \\
        --report skip \\
        --fl_count $flnc_count \\
        "\$ISOFORMS_INPUT" \\
        $reference_gtf \\
        $reference_fasta \\
        $args
    
    mv ${prefix}.transcriptome_classification.txt ${prefix}.transcriptome.SQANTI_classification.txt
    mv ${prefix}.transcriptome_corrected.gtf ${prefix}.transcriptome.gtf
    mv ${prefix}.transcriptome_corrected.fasta ${prefix}.transcriptome.fasta
    mv ${prefix}.transcriptome_junctions.txt ${prefix}.transcriptome.junctions.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.2.2
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcriptome.SQANTI_classification.txt
    touch ${prefix}.transcriptome.gtf
    touch ${prefix}.transcriptome.fasta
    touch ${prefix}.transcriptome.junctions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.2.2
    END_VERSIONS
    """
}
