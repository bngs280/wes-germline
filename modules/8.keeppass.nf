process KeepPASS {
    tag "Extract PASS variants on ${sample_id}"
    publishDir "${params.outdir}/6.filterDeepvariant", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path("${sample_id}_raw.vcf.gz")
    tuple val(sample_id), path(vcf)//worked
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS.vcf.gz"), emit: final_vcfs

    script: 
    """
    bcftools view -f PASS ${vcf} > ${sample_id}_filtered_PASS.vcf
    bgzip -c ${sample_id}_filtered_PASS.vcf > ${sample_id}_filtered_PASS.vcf.gz
    """
}
