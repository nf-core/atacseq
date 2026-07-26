//
// Uncompress and prepare reference genome files
//

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_TSS_BED
    GUNZIP as GUNZIP_BLACKLIST } from '../../modules/nf-core/gunzip/main'

include {
    UNTAR as UNTAR_BWA_INDEX
    UNTAR as UNTAR_BOWTIE2_INDEX
    UNTAR as UNTAR_CHROMAP_INDEX
    UNTAR as UNTAR_STAR_INDEX    } from '../../modules/nf-core/untar/main'

include { GFFREAD              } from '../../modules/nf-core/gffread/main'
include { SAMTOOLS_FAIDX       } from '../../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX            } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD        } from '../../modules/nf-core/bowtie2/build/main'
include { CHROMAP_INDEX        } from '../../modules/nf-core/chromap/index/main'
include { KHMER_UNIQUEKMERS    } from '../../modules/nf-core/khmer/uniquekmers/main'

include { STAR_GENOMEGENERATE      } from '../../modules/local/star_genomegenerate'
include { GTF2BED                  } from '../../modules/local/gtf2bed'
include { GENOME_BLACKLIST_REGIONS } from '../../modules/local/genome_blacklist_regions'
include { GET_AUTOSOMES            } from '../../modules/local/get_autosomes'
include { TSS_EXTRACT              } from '../../modules/local/tss_extract'

workflow PREPARE_GENOME {
    take:
    genome             //  string: genome name
    genomes            //     map: genome attributes
    prepare_tool_index //  string: tool to prepare index for
    fasta              //    path: path to genome fasta file
    gtf                //    file: /path/to/genome.gtf
    gff                //    file: /path/to/genome.gff
    blacklist          //    file: /path/to/blacklist.bed
    gene_bed           //    file: /path/to/gene.bed
    tss_bed            //    file: /path/to/tss.bed
    mito_name          //  string: name of mitochondrial chromosome
    keep_mito          // boolean: keep mitochondrial chromosome
    bwa_index          //    file: /path/to/bwa/index/
    bowtie2_index      //    file: /path/to/bowtie2/index/
    chromap_index      //    file: /path/to/chromap/index/
    star_index         //    file: /path/to/star/index/
    macs_gsize         // integer: MACS3 genome size
    read_length        // integer: read length

    main:

    //
    // Uncompress genome fasta file if required
    //
    ch_fasta = channel.empty()
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { tuple -> tuple[1] }.first()
    } else {
        ch_fasta = channel.value(file(fasta, checkIfExists: true))
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_gtf = channel.value(file(gtf, checkIfExists: true))
        }
    } else if (gff) {
        if (gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], gff ] ).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_gff = channel.value(file(gff, checkIfExists: true))
        }
        ch_gtf      = GFFREAD ( ch_gff.map { gff_file -> [ [:], gff_file ] }, [] ).gtf.map { _meta, gtf_file -> gtf_file }
    }

    //
    // Uncompress blacklist file if required
    //
    ch_blacklist = channel.empty()
    if (blacklist) {
        if (blacklist.endsWith('.gz')) {
            ch_blacklist = GUNZIP_BLACKLIST ( [ [:], blacklist ] ).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_blacklist = channel.value(file(blacklist, checkIfExists: true))
        }
    }

    // Uncompress gene BED annotation file or create from GTF if required
    //
    // If --gtf is supplied along with --genome
    // Make gene bed from supplied --gtf instead of using iGenomes one automatically
    def make_bed = false
    if (!gene_bed) {
        make_bed = true
    } else if (genome && gtf) {
        if (genomes[ genome ].gtf != gtf) {
            make_bed = true
        }
    }

    if (make_bed) {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
    } else {
        if (gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], params.gene_bed ] ).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_gene_bed = channel.value(file(gene_bed, checkIfExists: true))
        }
    }

    if (!tss_bed) {
        ch_tss_bed = TSS_EXTRACT ( ch_gene_bed ).tss
    } else {
        if (tss_bed.endsWith('.gz')) {
            ch_tss_bed = GUNZIP_TSS_BED ( [ [:], tss_bed ] ).gunzip.map { tuple -> tuple[1] }
        } else {
            ch_tss_bed = channel.value(file(tss_bed, checkIfExists: true))
        }
    }

    //
    // Create chromosome sizes file
    //
    SAMTOOLS_FAIDX ( ch_fasta.map { item -> [ [:], item, [] ] }, true )
    ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { tuple -> tuple[1] }.first()
    ch_fai         = SAMTOOLS_FAIDX.out.fai.map { tuple -> tuple[1] }.first()

    //
    // Create autosomal chromosome list for ataqv
    //
    ch_genome_autosomes = channel.empty()
    GET_AUTOSOMES (
        ch_fai
    )
    ch_genome_autosomes = GET_AUTOSOMES.out.txt


    //
    // Prepare genome intervals for filtering by removing regions in blacklist file
    //
    ch_genome_filtered_bed = channel.empty()
    GENOME_BLACKLIST_REGIONS (
        ch_chrom_sizes,
        ch_blacklist.ifEmpty([]),
        mito_name ?: '',
        keep_mito
    )
    ch_genome_filtered_bed = GENOME_BLACKLIST_REGIONS.out.bed

    //
    // Uncompress BWA index or generate from scratch if required
    //
    ch_bwa_index = channel.empty()
    if (prepare_tool_index == 'bwa') {
        if (bwa_index) {
            if (bwa_index.endsWith('.tar.gz')) {
                ch_bwa_index = UNTAR_BWA_INDEX ( [ [:], bwa_index ] ).untar
            } else {
                ch_bwa_index = [ [:], file(params.bwa_index, checkIfExists: true)]
            }
        } else {
            ch_bwa_index = BWA_INDEX ( ch_fasta.map { item -> [ [:], item ] } ).index
        }
    }

    //
    // Uncompress Bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = channel.empty()
    if (prepare_tool_index == 'bowtie2') {
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [:], bowtie2_index ] ).untar
            } else {
                ch_bowtie2_index = [ [:], file(bowtie2_index, checkIfExists: true) ]
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta.map { item -> [ [:], item ] } ).index
        }
    }

    //
    // Uncompress CHROMAP index or generate from scratch if required
    //
    ch_chromap_index = channel.empty()
    if (prepare_tool_index == 'chromap') {
        if (chromap_index) {
            if (chromap_index.endsWith('.tar.gz')) {
                ch_chromap_index = UNTAR_CHROMAP_INDEX ( [ [:], chromap_index ] ).untar
            } else {
                ch_chromap_index = [ [:], file(chromap_index, checkIfExists: true) ]
            }
        } else {
            ch_chromap_index = CHROMAP_INDEX ( ch_fasta.map { item -> [ [:], item ] } ).index
        }
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = channel.empty()
    if (prepare_tool_index == 'star') {
        if (star_index) {
            if (star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], star_index ] ).untar.map{ tuple -> tuple[1] }
            } else {
                ch_star_index = channel.value(file(star_index, checkIfExists: true))
            }
        } else {
            ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
        }
    }

    //
    // Estimate MACS3 genome size
    //
    ch_macs_gsize = macs_gsize
    if (!macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_fasta.map { item -> [ [:], item ] },
            read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { _meta, kmers -> kmers.text.trim() }
    }

    emit:
    fasta         = ch_fasta                      //    path: genome.fasta
    fai           = ch_fai                        //    path: genome.fai
    gtf           = ch_gtf                        //    path: genome.gtf
    gene_bed      = ch_gene_bed                   //    path: gene.bed
    tss_bed       = ch_tss_bed                    //    path: tss.bed
    chrom_sizes   = ch_chrom_sizes                //    path: genome.sizes
    filtered_bed  = ch_genome_filtered_bed        //    path: *.include_regions.bed
    bwa_index     = ch_bwa_index                  //    path: bwa/index/
    bowtie2_index = ch_bowtie2_index              //    path: bowtie2/index/
    chromap_index = ch_chromap_index              //    path: genome.index
    star_index    = ch_star_index                 //    path: star/index/
    autosomes     = ch_genome_autosomes           //    path: *.autosomes.txt
    macs_gsize    = ch_macs_gsize                 // integer: MACS3 genome size
}
