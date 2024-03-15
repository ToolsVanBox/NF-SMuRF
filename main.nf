
nextflow.enable.dsl=2

include { SMuRF } from './workflows/smurf.nf'

workflow {
    SMuRF()
}