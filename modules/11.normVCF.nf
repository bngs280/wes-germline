process VCFnorm {
    tag "Normalize vcf form TMB on ${sample_id}"
    publishDir "${params.outdir}/9.TMBScore", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(vcf)
    val(params.ref)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered_norm.vcf.gz"), emit: tmb_score
   // tuple val(sample_id), path("${sample_id}.TMB_results.log"), emit: tmb_score

    script:
    """
    bcftools norm -m- -f ${params.ref} -o ${sample_id}.filtered_norm.vcf.gz ${vcf}
    """
}
//bcftools norm -f ${params.refhg38} -o ${sample_id}.filtered_norm.vcf.gz ${vcf} \