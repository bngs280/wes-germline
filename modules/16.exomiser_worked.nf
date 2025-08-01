process MODIFY_CONFIG {
    publishDir "${params.outdir}/16.prioritised", mode: 'copy'

    input:
    path template_yml
    val phenotypes

    output:
    path 'modified-config.yml', emit: config

    script:
    //mkdir -p ${params.outdir}/16.prioritised
    //def phenotypes_list = phenotypes.collect { """    - "${it}"""" }.join('\n')
    def phenotypes_list = phenotypes.collect { "  - '${it}'" }.join('\\n')
    """
    mkdir -p "${params.outdir}/16.prioritised"
    awk '
    /^analysis:/ {
        print \$0
        print "  proband:"
        print "  hpoIds:"
        print "${phenotypes_list}"
        next
    }
    { print }
    ' $template_yml > modified-config.yml
    """
}

process run_prioritizer {
    publishDir "${params.outdir}/16.prioritised", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_file)
    path config_yml
    val genome
    path exomiser_properties

    output:
    path "${sample_id}/*"
    tuple val(sample_id), path("${sample_id}/*variants.tsv"), emit: variants_tsv
    tuple val(sample_id), path("${sample_id}/*exomiser.json"), emit: variants_json

    script:
    """
    mkdir -p "${sample_id}"
    java -jar ${params.exomiser_jar} --analysis ${config_yml} --output-directory ${sample_id} --assembly ${genome} --vcf ${vcf_file} --spring.config.location=${exomiser_properties} --output-filename ${sample_id}
    """
}

