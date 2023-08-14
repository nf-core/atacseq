process GENOME_BLACKLIST_REGIONS {
    tag "$sizes"

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0':
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    path sizes
    path blacklist
    val mito_name
    val keep_mito

    output:
    path '*.whitelist_regions.bed',     emit: whitelist_bed
    path '*.blacklist_regions.bed',     emit: blacklist_bed
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def whitelist = "${sizes.simpleName}.whitelist_regions.bed"
    def blacklist_out = "${sizes.simpleName}.blacklist_regions.bed"
    def name_filter = mito_name ? "| awk '\$1 !~ /${params.mito_name}/ {print \$0}'": ''
    def mito_filter = keep_mito ? '' : name_filter
    if (blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes $mito_filter > $whitelist
        complementBed -i $whitelist -g $sizes > $blacklist_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes $mito_filter > $whitelist
        complementBed -i $whitelist -g $sizes > $blacklist_out

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    }
}
