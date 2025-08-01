process vcfvep_process {
    tag "VCF-VEP merge on ${sample_id}"
    publishDir "${params.outdir}/8.VEPVCF", mode: 'copy'
    
    //conda 'pandas numpy'  // Add any dependencies your script needs
    
    input:
    tuple val(sample_id), path(vep_file), path(vcf_file)
    tuple val(sample_id), path(acmg_file)
    val(params.gender)
    val(params.clingen_file)
    
    output:
    tuple val(sample_id), path("${sample_id}*vcfvep.txt"), emit: vcf_output
    tuple val(sample_id), path("${sample_id}*output.csv"), emit: csv_output
    
    script:
    """
    python3 /usr/src/app/modules/vepvcf.py \
      --vep-file       ${vep_file} \
      --vcf-file       ${vcf_file} \
      --acmg-file      ${acmg_file} \
      --clingen-file   ${params.clingen_file} \
      --merged-vcf-vep ${sample_id}_vcfvep.txt \
      --output-file    ${sample_id}_output.csv \
      --sample-sex     ${params.gender}
    """
}
