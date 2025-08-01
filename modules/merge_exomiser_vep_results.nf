process mergeExomiserVEPResults {
    tag "Final P2G merge on ${sample_id}"
    publishDir "${params.outdir}/final_results", mode: 'copy'
    
    conda 'pandas numpy'
    
    input:
    tuple val(sample_id), path(vcf_output), path(variants_tsv), path(variants_json), path(snv_variants)
    
    output:
    tuple val(sample_id), path("${sample_id}_P2G_merged_results.csv"), emit: final_p2g_results
    
    script:
    """
    mkdir -p annotation_dir
    mkdir -p patient_dir
    
    python3 /usr/src/app/modules/exomiser_postprocess_trios.py \
        --vep        ${vcf_output} \
        --json       ${variants_json} \
        --variants   ${variants_tsv} \
        --output     ${sample_id}_P2G 

    """
}
