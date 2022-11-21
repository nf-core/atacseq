/*
 * Create IGV session file
 */
process IGV {

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3':
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path fasta
    path ("${bigwig_library_publish_dir}/*")
    path ("${peak_library_publish_dir}/*")
    path ("${consensus_library_publish_dir}/*")
    path ("${bigwig_replicate_publish_dir}/*")
    path ("${peak_replicate_publish_dir}/*")
    path ("${consensus_replicate_publish_dir}/*")
    val bigwig_library_publish_dir
    val peak_library_publish_dir
    val consensus_library_publish_dir
    val bigwig_replicate_publish_dir
    val peak_replicate_publish_dir
    val consensus_replicate_publish_dir
    // TODO: path differential_peaks from ch_macs_consensus_deseq_comp_igv.collect().ifEmpty([])

    output:
    // Publish fasta file while copyTo fails when the source and destination buckets are in different regions
    path "*files.txt"  , emit: txt
    path "*.xml"       , emit: xml
    path fasta         , emit: fasta
    path "versions.yml", emit: versions

    script: // scripts are bundled with the pipeline in nf-core/atacseq/bin/
    """
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$bigwig_library_publish_dir" || test \$? = 1; } > mLb_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$peak_library_publish_dir" || test \$? = 1; } > mLb_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,0" \\; | { grep "^$consensus_library_publish_dir" || test \$? = 1; } > mLb_bed.igv.txt
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$bigwig_replicate_publish_dir" || test \$? = 1; } > mRp_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$peak_replicate_publish_dir" || test \$? = 1; } > mRp_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,0" \\; | { grep "^$consensus_replicate_publish_dir" || test \$? = 1; } > mRp_bed.igv.txt

    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/${fasta.getName()} --path_prefix '../../'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
