process SV_som {
    tag "Structural Variants on ${sample_id}"
    publishDir "${params.outdir}/10.SV_somatic", mode: 'copy'
    cpus 4

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)

    output:
    path "${sample_id}", emit: sv_TUmsomatic

    script:
    """
    # Activate the environment for Manta
    source activate p2Manta

    # Create a temporary directory in the work directory for Manta to run
    #mkdir -p 10.SV_somatic
    mkdir ${sample_id}

    # Run Manta config and workflow
    configManta.py --tumorBam ${bam} --referenceFasta ${params.ref} --runDir ${sample_id}
    ${sample_id}/runWorkflow.py
    """
}


process SV_germline {
    tag "Structural Variants on ${sample_id}"
    publishDir "${params.outdir}/10.SV_germline", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)

    output:
    path "${sample_id}", emit: sv_germline

    script:
    """
    # Activate the environment for Manta
    source activate p2Manta

    # Create a temporary directory in the work directory for Manta to run
    #mkdir -p 10.SV_germline
    mkdir ${sample_id}

    # Run Manta config and workflow
    configManta.py --bam ${bam} --referenceFasta ${params.ref} --runDir ${sample_id}
    ${sample_id}/runWorkflow.py -m local
    """
}
