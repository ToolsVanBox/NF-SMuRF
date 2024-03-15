process SMURF {
  tag "$meta.id"
  label 'process_single'
  container = 'docker://vanboxtelbioinformatics/smurf:3.0.2'

  input:
    tuple val(meta), path(vcf), path(tbi), path(bams), path(bai)    

  output:
    tuple val(meta), path("*SMuRF.vcf"), emit: smurf_vcf
    tuple val(meta), path("*SMuRF.filtered.vcf"), emit: smurf_filtered_vcf
    tuple val(meta), path("*.pdf"), emit: smurf_pdf

  script:
    b = bams ? ' -b ' + bams.join(' -b ') : ''
    n = meta.bulk_names ? ' -n ' + meta.bulk_names.join(' -n ') : ''
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config = task.ext.config 
    
    """
    host=\$(hostname)
    echo \${host}

    export projectDir=${projectDir}

    python /smurf/SMuRF.py \
    -i ${vcf} \
    ${b} \
    ${n} \
    -t ${task.cpus} \
    -c ${config}

    """
}
