process somVarCall {
    tag "Somatic Variant on ${sample_id}"
    publishDir "${params.outdir}/5.variantM", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai")
    tuple val(sample_id), path(bam), path(bai) //worked
    val(params.ref)
    //val(params.bed)
    val(params.bed)
    val(params.gnomad)
    val(params.db1000g)

    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz"), path("${sample_id}_raw.vcf.gz.tbi"), path("${sample_id}_raw.vcf.gz.stats"), emit: raw_vcfs

    script:
    """
    java -jar /usr/src/app/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Mutect2 \
        -R ${params.ref} \
        --native-pair-hmm-threads ${task.cpus} \
        -L ${params.bed} \
        --germline-resource ${params.gnomad} \
        --panel-of-normals ${params.db1000g} \
        -I ${bam} \
        -O ${sample_id}_raw.vcf.gz
    """
}
   /*-L ${params.bed} \
   // Check if the genome reference is hg19, and if so, process the hg19 BED file
    if (params.genome == 'hg19') {
        """
        # Preprocess hg19 BED file to match the expected format for Mutect2
        grep -E '^([1-9]|1[0-9]|2[0-2]|X|Y|M)\\b' ${params.bed} > hg37.bed
        """
    } else {
        """
        # If the genome is hg38, use the BED file as-is
        ${params.bed} = ${params.bedhg38}
        """
    }
    */