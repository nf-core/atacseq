process GET_AUTOSOMES {
    tag "$fai"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path fai

    output:
    path '*.txt'       , emit: txt
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/atacseq/bin/
    """
    get_autosomes.py \\
        $fai \\
        ${fai.baseName}.autosomes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
