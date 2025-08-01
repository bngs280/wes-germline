process germVarCall {
    //container 'germline_pipeline:1.1'
    tag "Germline Variant on ${sample_id}"
    publishDir "${params.outdir}/5.variantM", mode: 'copy'
    cpus 20

    input:
    //tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai")
    tuple val(sample_id), path(bam), path(bai) //worked
    path(params.ref)
    path(params.bed)
    
   
    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf.gz"), path("${sample_id}_raw.vcf.gz.tbi"), emit: raw_vcfs

    script:
    """
        echo "=== Diagnostic Information ==="
        echo "Sample ID: ${sample_id}"
        echo "BAM File: ${bam}"
        echo "BAI File: ${bai}"
       
        echo "Reference FASTA (params.ref): ${params.ref}"
        echo "=== End of Diagnostic Information ==="    
    /opt/deepvariant/bin/run_deepvariant \
    --model_type WES \
    --ref ${params.ref} \
    --reads ${bam} \
    --regions ${params.bed} \
    --output_vcf ${sample_id}_raw.vcf.gz \
    --output_gvcf ${sample_id}_raw.gvcf.gz \
    --num_shards ${task.cpus} \
    --intermediate_results_dir ./
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
     --regions ${params.bed} \
      echo "BED File: ${params.bed}"
    */
