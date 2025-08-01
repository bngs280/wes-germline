process acmgpre {
    tag { "${sample_id}" }
    publishDir "${params.outdir}/20.ACMG", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_file)
    path(params.cache_bias)
    val(params.sd_bias)
    val(params.ref_bias)

    output:
    tuple val(sample_id), path("${sample_id}.json.gz"), emit: acmg_pre

    script:
    """
    dotnet /usr/src/app/Nirvana-v3.18.1/Nirvana.dll \
        --cache ${params.cache_bias} \
        --sd ${params.sd_bias} \
        --ref ${params.ref_bias} \
        --in ${vcf_file} \
        --out ${sample_id}
    """
}

process acmgclass {
    tag { "${sample_id}" }
    publishDir "${params.outdir}/20.ACMG", mode: 'copy'

    input:
    tuple val(sample_id), path(json_file)
    path(params.parameter_bias)

    output:
    tuple val(sample_id), path("${sample_id}_ACMG.tsv"), emit: acmg_class

    script:
    """
    python3 /usr/src/app/BIAS-2015/bias_2015.py \
        ${json_file} \
        ${params.parameter_bias} \
        ${sample_id}_ACMG.tsv
    """
}
