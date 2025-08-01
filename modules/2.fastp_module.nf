process trim_with_fastp_paired {
    tag "FASTp on ${sample_id} (Paired-End)"
    publishDir "${params.outdir}/2.fastpM", mode: 'copy'

    input:
    tuple val(sample_id), path(read_pairs_pe)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz"), emit: fastp_out_pe

    when:
    params.platform == 'illumina_pe' || params.platform == 'bgi_pe'
    
    script:
    """
    fastp --detect_adapter_for_pe -q 20 \
          -i ${read_pairs_pe[0]} \
          -I ${read_pairs_pe[1]} \
          -o ${sample_id}_R1_trimmed.fastq.gz \
          -O ${sample_id}_R2_trimmed.fastq.gz \
          --html ${sample_id}_fastp.html \
          --report_title "Quality Control for ${sample_id}"
    """
}

process trim_with_fastp_single {
    tag "FASTp on ${sample_id} (Single-End)"
    publishDir "${params.outdir}/2.fastpM", mode: 'copy'

    input:
    tuple val(sample_id), path(read_pairs_se)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: fastp_out_se

    when:
    params.platform == 'illumina_se' || params.platform == 'nanopore' || params.platform == 'bgi_se'

    script:
    """
    fastp -q 20 \
          -i ${read_pairs_se} \
          -o ${sample_id}_trimmed.fastq.gz \
          --html ${sample_id}_fastp.html \
          --report_title "Quality Control for ${sample_id}"
    """
}
