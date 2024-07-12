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
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pe     = meta.single_end ? '' : '-pc'
    def buffer = task.memory.toGiga().intdiv(2)
    """
    SCALE_FACTOR=\$(grep '[0-9] mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt

    bedtools \\
        genomecov \\
        -ibam $bam \\
        -bg \\
        -scale \$SCALE_FACTOR \\
        $args \\
    > tmp.bg

    ## ref: https://www.biostars.org/p/66927/
    ## ref in nf-core: https://github.com/nf-core/hicar/blob/d2d17a924e42d6f88640b79d48d8b332f33a953f/modules/local/atacreads/bedsort.nf#L23-L29
    LC_ALL=C sort \\
        --parallel=$task.cpus \\
        --buffer-size=${buffer}G \\
        -k1,1 -k2,2n \\
        $args2 \\
        tmp.bg \\
    > ${prefix}.bedGraph

    rm tmp.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        sort: \$(sort --version | head -n 1 | awk '{print \$4;}')
    END_VERSIONS
    """
}
