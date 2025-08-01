process mito_variant {
    tag "Mitochondrial variant calling on ${sample_id}"
    publishDir "${params.outdir}/12.mito_variants", mode : 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.mito_reference)
    val(params.mito_annotation)

    output:
    path "${sample_id}", emit:mitovariants

    script:
    """
    mkdir -p ${sample_id}
    /opt/mutserve/./mutserve call ${bam} --baseQ 20 --deletions --insertions --mode mutserve --reference ${params.mito_reference} --output ${sample_id}/${sample_id}.txt
    /opt/mutserve/./mutserve annotate --input ${sample_id}/${sample_id}.txt --annotation ${params.mito_annotation} --output ${sample_id}/${sample_id}_annotated.txt
    /opt/mutserve/./mutserve report --input ${sample_id}/${sample_id}_annotated.txt --output ${sample_id}/${sample_id}_report.html
    """
}
   
