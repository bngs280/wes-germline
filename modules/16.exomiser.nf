process MODIFY_CONFIG {
    publishDir "${params.outdir}/16.prioritised", mode: 'copy'

    input:
    path template_yml
    val hpoIds

    output:
    path 'modified-config.yml', emit: config

    script:
    //mkdir -p ${params.outdir}/16.prioritised
    //def phenotypes_list = phenotypes.collect { """    - "${it}"""" }.join('\n')
    //def phenotypes_list = phenotypes.collect { "  - '${it}'" }.join('\\n')
    def phenotypes_list = hpoIds.collect { "      - '${it}'" }.join('\\n')
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
    /^  hpoIds:/ { next }  # Skip the original hpoIds line
    /^- / { next }  # Skip any existing HPO ID lines
    { print }
    ' $template_yml > modified-config.yml
    """
}
/*
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
    java -jar ${params.exomiser_jar} --analysis ${config_yml} --output-directory ${sample_id} --assembly ${genome} --vcf ${vcf_file} --spring.config.location=${exomiser_properties} --output-filename ${sample_id}
    """
}
*/

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
    7z x -aos /usr/src/app/exomiser/db/2406_phenotype.zip -o/usr/src/app/exomiser/db
    # Determine which database to unzip based on the selected genome
    if [[ "$genome" == "hg19" ]]; then
        echo "Unzipping hg19 database..."
        unzip -o /usr/src/app/exomiser/db/2406_hg19.zip -d /usr/src/app/exomiser/db
    elif [[ "$genome" == "hg38" ]]; then
        echo "Unzipping hg38 database..."
        unzip -o /usr/src/app/exomiser/db/2406_hg38.zip -d /usr/src/app/exomiser/db
    else
        echo "Invalid genome specified: $genome"
        exit 1
    fi

    # Run Exomiser with the unzipped database
    java -jar ${params.exomiser_jar} --analysis ${config_yml} --output-directory ${sample_id} \
        --assembly ${genome} --vcf ${vcf_file} --spring.config.location=${exomiser_properties} \
        --output-filename ${sample_id}
    """
}

