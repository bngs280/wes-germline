process runG2P {
    tag "G2P for ${sample_id}"
    publishDir "${params.outdir}/final_results", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_output)

    output:
    tuple val(sample_id), path("${sample_id}_G2P_filtered.csv"), emit: final_g2p_score
    tuple val(sample_id), path("${sample_id}_all_filtered.csv"), emit: primary_g2p_score

    script:
    """
    python3 /usr/src/app/modules/filtG2P3.py \
        --primary_filtered  ${sample_id}_all_filtered.csv \
        --output_file       ${sample_id}_G2P_filtered.csv \
        --merged_vcf_vep    ${vcf_output} \
        --omim_file         ${params.omim} \
        --population        ${params.frequency_column} 

    """
}
// --population       ${params.frequency_column}
//--moi_file         ${params.moi_dk} \

//--primary_filtered "MGL02_primary_filtered_1107g2p.csv" --output_file "MGL02_output_filtg2p1107.csv" --merged_vcf_vep "MGL02_vcf_vep.txt" --omim_file "gene.tsv" --population "SAS"
