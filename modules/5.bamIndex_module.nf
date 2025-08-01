process samtool_index {
    tag "Index final md BAM on ${sample_id}"
    publishDir "${params.outdir}/4.markDupliM", mode: 'copy'
    cpus 16

    input:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam")

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai"), emit: markdup_bams

    script:
    """
    samtools index -@ ${task.cpus} ${sample_id}_sorted_md.bam

    """
}
