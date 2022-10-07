/*
 * Create IGV session file
 */
process IGV {

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3':
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    val aligner_dir
    val peak_dir
    path fasta
    path ("${aligner_dir}/mergedLibrary/bigwig/*")
    path ("${aligner_dir}/mergedLibrary/macs2/${peak_dir}/*")
    path ("${aligner_dir}/mergedLibrary/macs2/${peak_dir}/consensus/*")
    path ("mappings_lib/*")
    path ("${aligner_dir}/mergedReplicate/bigwig/*")
    path ("${aligner_dir}/mergedReplicate/macs2/${peak_dir}/*")
    path ("${aligner_dir}/mergedReplicate/macs2/${peak_dir}/consensus/*")
    path ("mappings_rep/*")

    output:
    path "*files.txt"  , emit: txt
    path "*.xml"       , emit: xml
    path "versions.yml", emit: versions

    script: // scripts are bundled with the pipeline in nf-core/atacseq/bin/
    def consensus_dir_lib = "${aligner_dir}/mergedLibrary/macs2/${peak_dir}/consensus/*"
    def consensus_dir_rep = "${aligner_dir}/mergedReplicate/macs2/${peak_dir}/consensus/*"
    """
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; > mLb_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; > mLb_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,0" \\; | { grep "^$consensus_dir_lib" || test \$? = 1; } > mLb_bed.igv.txt

    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; > mRp_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; > mRp_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,0" \\; | { grep "^$consensus_dir_rep" || test \$? = 1; } > mRp_bed.igv.txt

    cat mappings_lib/* > replace_paths.txt
    cat mappings_rep/* >> replace_paths.txt

    cat *.txt > igv_files_orig.txt
    igv_files_to_session.py igv_session.xml igv_files_orig.txt replace_paths.txt ../../genome/${fasta.getName()} --path_prefix '../../'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
// find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
// find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedReplicate/macs/${PEAK_TYPE}/consensus/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
