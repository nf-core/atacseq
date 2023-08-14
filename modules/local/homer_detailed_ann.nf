process HOMER_DETAIL_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    // configureHomer.pl does not work with docker or singularity, so we use conda (see: https://github.com/bioconda/bioconda-recipes/blob/master/recipes/homer/README.txt)
    conda "bioconda::homer=4.11"

    input:
    tuple val(meta), path(peak)
    val genome

    output:
    tuple val(meta), path("*annotatePeaks.detailed.txt"), emit: txt
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '4.11' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    installed=\$(perl \$CONDA_PREFIX/share/homer*/configureHomer.pl -list)
    line=\$(echo "\$installed" | grep "$genome")
    yes_or_no=\${line:0:1}

    if [ \$yes_or_no != "+" ]; then
        perl \$CONDA_PREFIX/share/homer*/configureHomer.pl -install $genome
    fi

    annotatePeaks.pl \\
        $peak \\
        $genome \\
        $args \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.detailed.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
