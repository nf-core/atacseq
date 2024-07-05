process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.1 conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bedtools_coreutils:ba273c06a3909a15':
        'community.wave.seqera.io/library/bedtools_coreutils:a623c13f66d5262b' }"

    input:
    tuple val(meta), path(bam), path(flagstat)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.txt")     , emit: scale_factor
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pe     = meta.single_end ? '' : '-pc'
    """
    SCALE_FACTOR=\$(grep '[0-9] mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -scale \$SCALE_FACTOR \\
        $pe \\
        $args \\
    > tmp.bg

    sort -T '.' -k1,1 -k2,2n tmp.bg > ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
