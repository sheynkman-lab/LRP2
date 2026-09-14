process DOWNLOAD_REMOTE_FILE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://docker.io/biocontainers/biocontainers:latest"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("${output_filename}"), emit: file
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Extract actual filename from URL
    def url_str = url.toString()
    if (url_str.contains('?file=')) {
        // MASSIVE-style URL with filename in query parameter
        def query_file_param = url_str.split('\\?file=')[1]
        output_filename = query_file_param.tokenize('/')[-1]
    } else {
        // Standard URL with filename in path
        def url_path = url_str.replaceAll(/\?.*$/, '')
        output_filename = url_path.tokenize('/')[-1]
    }

    """
    echo "Downloading remote file from: $url"
    echo "Saving as: ${output_filename}"

    # Try wget first (more commonly available), fallback to curl
    if command -v wget &> /dev/null; then
        echo "Using wget for download"
        wget -O "${output_filename}" "$url"
        download_tool="wget"
    elif command -v curl &> /dev/null; then
        echo "Using curl for download"
        curl -L -o "${output_filename}" "$url"
        download_tool="curl"
    else
        echo "ERROR: Neither wget nor curl available"
        exit 1
    fi

    if [ ! -f "${output_filename}" ]; then
        echo "ERROR: Failed to download file"
        exit 1
    fi

    file_size=\$(stat -c%s "${output_filename}" 2>/dev/null || stat -f%z "${output_filename}" 2>/dev/null || echo "0")
    echo "Downloaded file size: \${file_size} bytes"

    if [ "\${file_size}" -eq 0 ]; then
        echo "ERROR: Downloaded file is empty"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        download_tool: \${download_tool}
    END_VERSIONS
    """

    stub:
    def url_str = url.toString()
    if (url_str.contains('?file=')) {
        def query_file_param = url_str.split('\\?file=')[1]
        output_filename = query_file_param.tokenize('/')[-1]
    } else {
        def url_path = url_str.replaceAll(/\?.*$/, '')
        output_filename = url_path.tokenize('/')[-1]
    }

    """
    touch ${output_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: 7.68.0
    END_VERSIONS
    """
}
