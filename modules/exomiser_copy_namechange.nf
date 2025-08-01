process MODIFY_CONFIG {
    publishDir "${params.outdir}/16.prioritised", mode: 'copy'

    input:
    path template_yml
    val phenotypes
    val sample_id

    output:
    tuple val(sample_id), path 'modified-config.yml', emit: config

    script:
    def phenotypes_list = phenotypes.collect { "  - '${it}'" }.join('\n')

    """
    mkdir -p "${params.outdir}/16.prioritised"

    # Debug: Print variables to check their values
    echo "Sample ID: ${sample_id}"
    echo "Phenotypes List: ${phenotypes_list}"

    awk -v sample="${sample_id}" -v phenotypes_list="${phenotypes_list}" '
    /^analysis:/ {
        print \$0
        print "  proband:"
        print "  hpoIds:"
        print phenotypes_list
        print "outputOptions:"
        print "  outputPrefix: " sample
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

    script:
    """
    mkdir -p "${sample_id}"
    java -jar ${params.exomiser_jar} --analysis ${config_yml} --output-directory ${sample_id} --assembly ${genome} --vcf ${vcf_file} --spring.config.location=${exomiser_properties}
    """
}

