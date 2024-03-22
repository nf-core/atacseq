process MULTIQC {
    label 'process_medium'

    conda 'bioconda::multiqc=1.13'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13--pyhdfd78af_0' :
        'biocontainers/multiqc:1.13--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary

    path ('fastqc/*')
    path ('trimgalore/fastqc/*')
    path ('trimgalore/*')

    path ('alignment/library/*')
    path ('alignment/library/*')
    path ('alignment/library/*')

    path ('alignment/merged_library/unfiltered/*')
    path ('alignment/merged_library/unfiltered/*')
    path ('alignment/merged_library/unfiltered/*')
    path ('alignment/merged_library/unfiltered/picard_metrics/*')

    path ('alignment/merged_library/filtered/*')
    path ('alignment/merged_library/filtered/*')
    path ('alignment/merged_library/filtered/*')
    path ('alignment/merged_library/filtered/picard_metrics/*')

    path ('preseq/*')

    path ('deeptools/*')
    path ('deeptools/*')

    path ('macs2/merged_library/peaks/*')
    path ('macs2/merged_library/peaks/*')
    path ('macs2/merged_library/annotation/*')
    path ('macs2/merged_library/featurecounts/*')

    path ('genrich/merged_library/sep/peaks/*')
    path ('genrich/merged_library/sep/peaks/*')
    path ('genrich/merged_library/sep/annotation/*')

    path ('alignment/merged_replicate/*')
    path ('alignment/merged_replicate/*')
    path ('alignment/merged_replicate/*')
    path ('alignment/merged_replicate/picard_metrics/*')

    path ('macs2/merged_replicate/peaks/*')
    path ('macs2/merged_replicate/peaks/*')
    path ('macs2/merged_replicate/annotation/*')
    path ('macs2/merged_replicate/featurecounts/*')

    path ('deseq2_library/*')
    path ('deseq2_library/*')
    path ('deseq2_replicate/*')
    path ('deseq2_replicate/*')

    path ('genrich/merged_library/joint/peaks/*')
    path ('genrich/merged_library/joint/peaks/*')
    path ('genrich/merged_library/joint/annotation/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
