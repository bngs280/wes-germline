process align_bwa_paired {
    tag "BWA on ${sample_id} (Paired-End)"
    publishDir "${params.outdir}/3.alignmentM", mode: 'copy'
    cpus 20

    input:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"), path("${sample_id}_R2_trimmed.fastq.gz")
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bams_pe
    
    when:
    params.platform == 'illumina_pe' || params.platform == 'bgi_pe'
    
    script:
    """
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${sample_id}" \
            ${params.ref} ${sample_id}_R1_trimmed.fastq.gz ${sample_id}_R2_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam 
    """
}

process align_bwa_single {
    tag "BWA on ${sample_id} (Single-End)"
    publishDir "${params.outdir}/3.alignmentM", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz")
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bams_se

    when:
    params.platform == 'illumina_se' || params.platform == 'nanopore' || params.platform == 'bgi_se'

    script:
    """
    bwa mem -t ${task.cpus} -M \
            -R "@RG\\tID:${sample_id}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:${sample_id}" \
            ${params.ref} ${sample_id}_trimmed.fastq.gz | samtools view -@ ${task.cpus} -b - | samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam 
    """
}
