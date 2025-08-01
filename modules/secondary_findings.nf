process runSecondaryFindings {
    tag "Secondary findings for ${sample_id}"
    publishDir "${params.outdir}/secondary_findings", mode: 'copy'
    
    input:
    tuple val(sample_id), path(pvs1_lof_ch)
    tuple val(_), path(exomiser_tsv)
    tuple val(_), path(vcf_output)
    
    output:
    path "${sample_id}_secondary_findings.csv", emit: secondary_findings
    
    script:
    """
    python3 /usr/src/app/modules/script5.py \
      --vep         ${vcf_output} \
      --exomiser    ${exomiser_tsv} \
      --acmg        ${params.secondaryfindings_gene_file} \
      --lof         ${pvs1_lof_ch} \
      --report-file ${sample_id}_secondary_findings.csv
    """
}
