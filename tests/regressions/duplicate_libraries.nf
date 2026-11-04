nextflow.enable.dsl = 2

include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES } from '../../modules/nf-core/picard/markduplicates/main'

process MERGE_READ_GROUPS {
    container params.samtools_image
    input:
    path bams
    output:
    path 'combined.bam'
    script:
    """
    samtools merge -o combined.bam ${bams.sort().join(' ')}
    """
}

workflow FASTQ_ALIGN_CHROMAP {
    take:
    reads
    main:
    PICARD_ADDORREPLACEREADGROUPS(reads, channel.value([[:], [], []]))
    emit:
    bam = PICARD_ADDORREPLACEREADGROUPS.out.bam
}

workflow MERGED_LIBRARY_MARKDUPLICATES_PICARD {
    take:
    reads
    main:
    PICARD_MARKDUPLICATES(reads, channel.value([[:], [], []]))
    emit:
    bam = PICARD_MARKDUPLICATES.out.bam
}

workflow ATACSEQ {
    take:
    reads
    main:
    FASTQ_ALIGN_CHROMAP(reads)
    MERGE_READ_GROUPS(FASTQ_ALIGN_CHROMAP.out.bam.map { meta, bam -> bam }.collect())
    MERGED_LIBRARY_MARKDUPLICATES_PICARD(MERGE_READ_GROUPS.out.map { bam -> [[id: 'LIB_A_REP1', single_end: false], bam] })
    emit:
    bam = MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.bam
}

workflow NFCORE_ATACSEQ {
    take:
    reads
    main:
    ATACSEQ(reads)
    emit:
    bam = ATACSEQ.out.bam
}

workflow {
    reads = channel.of('LIB_A_REP1_T1', 'LIB_A_REP1_T2', 'LIB_B_REP1_T1')
        .map { id -> [[id: id, single_end: false], file("${params.fixture_dir}/${id}.bam", checkIfExists: true)] }
    NFCORE_ATACSEQ(reads)
    NFCORE_ATACSEQ.out.bam.subscribe { meta, bam ->
        file(params.result_path).text = groovy.json.JsonOutput.toJson([
            bam: bam.toString(), nextflow_version: nextflow.version.toString()
        ]) + '\n'
    }
}
