process TSS_EXTRACT {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path bed

    output:
    path "*.bed"       , emit: tss
    tuple val("${task.process}"), val('sed'), eval("sed --version 2>&1 | tr '\\n' ' ' | sed 's/^.*GNU sed) //; s/ .*\$//'"), topic: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $bed | awk -v FS='\t' -v OFS='\t' '{ if(\$6=="+") \$3=\$2+1; else \$2=\$3-1; print \$1, \$2, \$3, \$4, \$5, \$6;}' > ${bed.baseName}.tss.bed

    """
}
