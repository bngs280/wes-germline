process germvar {
    //container 'rgenxtools:latest'
    tag "Deepvariant on ${sample_id}"
    publishDir "${params.outdir}/12.germline", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(read_files), path(bai_file)
    val(params.refhg38)
    val(params.bedhg38)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), path("${sample_id}.gvcf.gz"), path("${sample_id}.gvcf.gz.tbi"), emit: germvarRes

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type WES \
    --ref ${params.refhg38} \
    --reads ${read_files} \
    --regions ${params.bedhg38} \
    --output_vcf ${sample_id}.vcf.gz \
    --output_gvcf ${sample_id}.gvcf.gz \
    --num_shards 4 \
    --intermediate_results_dir .
    """
}