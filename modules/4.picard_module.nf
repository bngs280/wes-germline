process markDup {
    tag "Mark Duplicate on ${sample_id}"
    publishDir "${params.outdir}/4.markDupliM", mode: 'copy'
    cpus 16

    input:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    output:
    tuple val(sample_id), path("${sample_id}_sorted_md.bam"), emit: markdup_bams

    script:
    """
    sambamba markdup -r -t ${task.cpus} ${sample_id}_sorted.bam ${sample_id}_sorted_md.bam

    """
}
/*
    sambamba markdup -r -t ${task.cpus} ${sample_id}_sorted.bam ${sample_id}_sorted_md.bam 

        java -Xmx50g -jar /usr/local/bin/picard.jar MarkDuplicates \
    	I=${sample_id}_sorted.bam \
    	O=${sample_id}_sorted_md.bam \
    	M=${sample_id}_sorted_md_metrics.txt
*/


