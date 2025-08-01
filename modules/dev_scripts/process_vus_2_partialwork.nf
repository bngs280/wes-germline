/*
process processVEPOutput {
    tag "Process VEP output for ${sample_id}"
    publishDir "${params.outdir}/10.vusprize", mode: 'copy'
    cpus 14

    input:
    tuple val(sample_id), path(vus_vep)

    output:
    tuple val(sample_id), path("$(sample_id)_processed.txt"), emit: processed_vusprize

    script:
    """
    grep -v '^##' ${sample_id}_filtered_VUS_VEP.txt > ${sample_id}_tmp.txt
    sed '1s/Eigen-phred_coding/Eigen-pred_coding/' ${sample_id}_tmp.txt > ${sample_id}_processed.txt
    rm ${sample_id}_tmp.txt

    micromamba run -n vusprize /usr/src/app/VusPrize/vusprize/Vusprize.py ${sample_id}_processed.txt > ${sample_id}_vusprize_output.txt
    rm ${sample_id}_processed.txt
    """
*/

process processVEPOutput {
    tag "Process VEP output for ${sample_id}"
    publishDir "${params.outdir}/10.VUSpriritised", mode: 'copy'
    cpus 14

    input:
    tuple val(sample_id), path(vus_vep)

    output:
    tuple val(sample_id), path("${sample_id}_vusScore.txt"), emit: processed_vusprize

    script:
    """
    # Step 1: Remove header lines
    grep -v '^##' ${vus_vep} > ${sample_id}_tmp.txt

    # Step 2: Replace incorrect column header
    sed '1s/Eigen-phred_coding/Eigen-pred_coding/' ${sample_id}_tmp.txt > ${sample_id}_processed.txt
    rm ${sample_id}_tmp.txt

    # Step 3: Run VusPrize.py with absolute paths
    unset PYTHONPATH
    INPUT_FILE="\$(realpath ${sample_id}_processed.txt)"
    OUTPUT_FILE="\$(realpath ${sample_id}_vusScore.txt)"
    cd /usr/src/app/VusPrize/vusprize
    micromamba run -n vusprize python VusPrize.py "\$INPUT_FILE" "\$OUTPUT_FILE"
    """
}
