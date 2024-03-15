
nextflow.enable.dsl=2

// Include module
include { SNPSIFT_SPLIT } from '../modules/nf-core/snpsift/split/main.nf'
include { SNPSIFT_SPLIT as SNPSIFT_JOIN_FILTERED_VCFS } from '../modules/nf-core/snpsift/split/main.nf'
include { SNPSIFT_SPLIT as SNPSIFT_JOIN_SMURF_VCFS } from '../modules/nf-core/snpsift/split/main.nf'
include { SMURF } from '../modules/local/SMuRF/main.nf'
include { TABIX_BGZIPTABIX } from '../modules/nf-core/tabix/bgziptabix/main.nf'
include { TABIX_BGZIP } from '../modules/nf-core/tabix/bgzip/main.nf'

workflow SMuRF {
    def vcf = file( params.vcf, checkIfExists: true )

    ch_input = Channel.value( vcf )
        .map{ vcf_file ->
         [ [ id: vcf_file.getBaseName(), split: true ], vcf_file ]
        }
        .branch{
           gz: it[1].toString() =~ /.*gz$/
           vcf: it[1].toString() =~ /.*vcf$/
        }

    ch_vcf = Channel.empty()
    ch_vcf = ch_vcf.mix( ch_input.vcf )
    
    TABIX_BGZIP( ch_input.gz )

    ch_vcf = ch_vcf.mix( TABIX_BGZIP.out.output )

    SNPSIFT_SPLIT( ch_vcf )

    ch_bams = extractBamsFromDir( params.bams, params.bulk_names )
        .groupTuple()

    ch_bgziptabix = SNPSIFT_SPLIT.out.out_vcfs
        .map{ meta, vcf_file ->
            [ vcf_file ]
        }
        .flatten()
        .map{ vcf_file -> 
            [ [ id: vcf_file.getBaseName(), split: true ], vcf_file ]
        }

    
    TABIX_BGZIPTABIX( ch_bgziptabix )

    ch_smurf = TABIX_BGZIPTABIX.out.gz_tbi
        .combine( ch_bams )
        .map{ meta, vcf_file, tbi, meta2, bam, bai ->
          meta = meta + [ bulk_names: meta2.bulk_names ]
          [ meta, vcf_file, tbi, bam, bai ]
        }

    SMURF( ch_smurf )

    ch_filtered_vcfs = SMURF.out.smurf_filtered_vcf
        .map{ meta, vcf_file -> 
            [ [ id: file(file(file(vcf_file.getBaseName()).getBaseName()).getBaseName()).getBaseName() ], vcf_file ]
        }
        .groupTuple()
    
    ch_smurf_vcfs = SMURF.out.smurf_vcf
        .map{ meta, vcf_file -> 
            [ [ id: file(file(file(vcf_file.getBaseName()).getBaseName()).getBaseName()).getBaseName() ], vcf_file ]
        }
        .groupTuple()

    SNPSIFT_JOIN_SMURF_VCFS( ch_smurf_vcfs )
    
}    

def extractBamsFromDir( bams_dir, bulk_names ) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    bams_dir = bams_dir.tokenize().collect{"$it/*.bam"}
    Channel
      .fromPath(bams_dir, type:'file')
      .ifEmpty { error "No bam files found in ${bams_dir}." }
      .map { bam_path ->
          bam_file = bam_path
          bai_file = bam_path+'.bai'
          [ [ bulk_names: bulk_names ], bam_file, bai_file]          
      }
}