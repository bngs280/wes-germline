process lofDetection {
    tag "LoF detection on ${sample_id}"
    publishDir "${params.outdir}/lof_detection", mode: 'copy'
    
    conda 'pandas numpy'
    
    input:
    tuple val(sample_id), path(vep_file)
    path clingen_hi_file
    path gene_data_file
    
    output:
    tuple val(sample_id), path("${sample_id}_lof_results.csv"), emit: lof_results
    
    script:
    """
    python3 /usr/src/app/modules/pvs1_pipeline.py \\
        --vep-file ${vep_file} \\
        --clingen-hi ${clingen_hi_file} \\
        --gene-data ${gene_data_file} \\
        --output-file ${sample_id}_lof_results.csv
    """
}
