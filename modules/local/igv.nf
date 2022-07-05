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
    path ("${bigwig_lib_publish_dir}/*")
    path ("${peak_lib_publish_dir}/*")
    path ("${consensus_lib_publish_dir}/*")
    path ("${bigwig_rep_publish_dir}/*")
    path ("${peak_rep_publish_dir}/*")
    path ("${consensus_rep_publish_dir}/*")
    val bigwig_lib_publish_dir
    val peak_lib_publish_dir
    val consensus_lib_publish_dir
    val bigwig_rep_publish_dir
    val peak_rep_publish_dir
    val consensus_rep_publish_dir
    // path differential_peaks from ch_macs_consensus_deseq_comp_igv.collect().ifEmpty([])

    output:
    path "*files.txt"  , emit: txt
    path "*.xml"       , emit: xml
    path "versions.yml", emit: versions

    script: // scripts are bundled with the pipeline in nf-core/atacseq/bin/
    """
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$bigwig_lib_publish_dir" || test \$? = 1; } > mLb_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$peak_lib_publish_dir" || test \$? = 1; } > mLb_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,0" \\; | { grep "^$consensus_lib_publish_dir" || test \$? = 1; } > mLb_bed.igv.txt
    find * -type l -name "*.bigWig" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$bigwig_rep_publish_dir" || test \$? = 1; } > mRp_bigwig.igv.txt
    find * -type l -name "*Peak" -exec echo -e ""{}"\\t0,0,178" \\; | { grep "^$peak_rep_publish_dir" || test \$? = 1; } > mRp_peaks.igv.txt
    find * -type l -name "*.bed" -exec echo -e ""{}"\\t0,0,0" \\; | { grep "^$consensus_rep_publish_dir" || test \$? = 1; } > mRp_bed.igv.txt

    cat *.txt > igv_files.txt
    igv_files_to_session.py igv_session.xml igv_files.txt ../../genome/${fasta.getName()} --path_prefix '../../'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
// find * -type f -name "${prefix}.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/"{}"\\t0,0,0" \\; > ${prefix}.bed.igv.txt
// find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/${antibody}/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt

// //L1016
//     find * -type f -name "*.bigWig" -exec echo -e "bwa/mergedLibrary/bigwig/"{}"\\t0,0,178" \\; > ${prefix}.bigWig.igv.txt
// //L1162
//     find * -type f -name "*.${PEAK_TYPE}" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/"{}"\\t0,0,178" \\; > ${prefix}_peaks.igv.txt
// //L1285
//     find * -type f -name "${prefix}.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/"{}"\\t0,0,0" \\; > ${prefix}.bed.igv.txt
// //L1403
// //     find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedLibrary/macs/${PEAK_TYPE}/consensus/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
// //L1604
//     find * -type f -name "*.bigWig" -exec echo -e "bwa/mergedReplicate/bigwig/"{}"\\t0,0,178" \\; > ${prefix}.bigWig.igv.txt
// //L1658
//     find * -type f -name "*.${PEAK_TYPE}" -exec echo -e "bwa/mergedReplicate/macs/${PEAK_TYPE}/"{}"\\t0,0,178" \\; > ${prefix}_peaks.igv.txt
// //L1780
//     find * -type f -name "${prefix}.bed" -exec echo -e "bwa/mergedReplicate/macs/${PEAK_TYPE}/consensus/"{}"\\t0,0,0" \\; > ${prefix}.bed.igv.txt
// //L1898
// //     find * -type f -name "*.FDR0.05.results.bed" -exec echo -e "bwa/mergedReplicate/macs/${PEAK_TYPE}/consensus/deseq2/"{}"\\t255,0,0" \\; > ${prefix}.igv.txt
