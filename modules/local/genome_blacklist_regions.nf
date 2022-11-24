process GENOME_BLACKLIST_REGIONS {
    tag "$sizes"

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    path sizes
    path blacklist

    output:
    path '*.bed'       , emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def file_out = "${sizes.simpleName}.include_regions.bed"
    def name_filter = params.mito_name ? "| awk '\$1 !~ /${params.mito_name}/ {print \$0}'": ''
    def mito_filter = params.keep_mito ? '' : name_filter
    if (blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes $mito_filter > $file_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes $mito_filter > $file_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
}
