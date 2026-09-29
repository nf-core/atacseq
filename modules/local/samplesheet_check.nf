process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions
    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/atacseq/bin/
    def args = task.ext.args ?: ''
    """
    check_samplesheet.py \\
        $samplesheet \\
        $args

    """
}
