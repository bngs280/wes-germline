process fastqc_reads {
    tag "FASTQC on ${sample_id}"
    publishDir "${params.outdir}/1.QCM", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    output:
    path "${sample_id}", emit: fastqc_out

    script:
    """
    mkdir ${sample_id}
    fastqc -o ${sample_id} -f fastq -q ${read_files}
    """
}

