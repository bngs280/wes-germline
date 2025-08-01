process msiScore {
    tag "MSI Score Calculation on ${sample_id}"
    publishDir "${params.outdir}/8.MSIScore", mode: 'copy'
    cpus 4

    input:
    //tuple val(sample_id), path("${sample_id}_sorted_md.bam"), path("${sample_id}_sorted_md.bam.bai")
    tuple val(sample_id), path(bam), path(bai) //worked
    val(params.msimodel)

    output:
    path("${sample_id}.MSIscore.log"), emit: msi_score

    script:
    """
    /usr/src/app/msisensor2/msisensor2 msi -b 10 -M ${params.msimodel} -t ${bam} -o ${sample_id}.MSIscore.log
    """
}
