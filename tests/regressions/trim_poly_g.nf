nextflow.enable.dsl = 2

include { TRIMGALORE } from '../../modules/nf-core/trimgalore/main'

workflow {
    TRIMGALORE(channel.of([
        [id: 'poly_g', single_end: false],
        [file(params.read1, checkIfExists: true), file(params.read2, checkIfExists: true)]
    ]))
    TRIMGALORE.out.reads.subscribe { row ->
        def (meta, reads) = row
        file(params.receipt).text = groovy.json.JsonOutput.toJson([
            meta: meta, reads: [reads].flatten().collect { read -> read.toString() }
        ]) + '\n'
    }
}
