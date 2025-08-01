process mutect2Filter {
    tag "Somatic Mutect2 Filter on ${sample_id}"
    publishDir "${params.outdir}/6.filterMutect", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path("${sample_id}_raw.vcf.gz")
    tuple val(sample_id), path(vcf), path(tbi), path(tsv) //worked

    val(params.ref)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), emit: filtered_vcfs

    script: 
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar FilterMutectCalls \
        -R ${params.ref} \
        -V ${vcf} \
        -O ${sample_id}_filtered.vcf.gz
    """
}
