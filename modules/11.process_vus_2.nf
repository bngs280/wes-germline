// First Process: Process VEP Output
process processVEPOutput {
    tag "Process VEP output for ${sample_id}"
    publishDir "${params.outdir}/10.VUSpriritised", mode: 'copy'
    cpus 14

    input:
    tuple val(sample_id), path(vus_vep)

    output:
    tuple val(sample_id), path("${sample_id}_processed.txt"), emit: processed_vcf

    script:
    """
    # Step 1: Remove header lines
    grep -v '^##' ${vus_vep} > ${sample_id}_tmp.txt

    # Step 2: Replace incorrect column header
    sed '1s/Eigen-phred_coding/Eigen-pred_coding/' ${sample_id}_tmp.txt > ${sample_id}_processed.txt
    rm ${sample_id}_tmp.txt
    """
}


// Second Process: Run VusPrize
process runVusPrioritization {
    tag "Run VusPrize for ${sample_id}"
    publishDir "${params.outdir}/10.VUSpriritised", mode: 'copy'
    cpus 14

    input:
    tuple val(sample_id), path(processed_vus) 
    path(params.vusprize_script)
    path(params.dependent_script)

    output:
    tuple val(sample_id), path("${sample_id}_vusScore.txt"), emit: vusprize_output

    script:
    """
    # Save the current working directory
    wd=\$(pwd)
    # Step 3: Run VusPrize.py
    unset PYTHONPATH
    # Copy the processed output file into the VusPrize directory
    cp ${processed_vus} /usr/src/app/VusPrize/vusprize/
    cd /usr/src/app/VusPrize/vusprize
    micromamba run -p /opt/miniconda3/envs/vusprize bash -c "
      python VusPrize.py \"\$(basename ${processed_vus})\" \"${sample_id}_vusScore.txt\"
    "
    # Copy the generated output file back to the original working directory
    cp ${sample_id}_vusScore.txt \$wd
    """
}

