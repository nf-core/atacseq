process FEATURECOUNTS_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path('featurecounts/*')

    output:
    tuple val(meta), path("*.featureCounts.tsv"), emit: counts
    tuple val("${task.process}"), val('sed'), eval("sed --version | sed '1!d;s/.*GNU sed) //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    featurecounts_merge.sh \\
        ${prefix}.featureCounts.tsv \\
        \$(ls featurecounts/*.featureCounts.tsv | sort)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.featureCounts.tsv
    """
}
