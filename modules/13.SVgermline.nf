/*process SV_som {
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
*/

// process SV_germline {
//     tag "Structural Variants on ${sample_id}"
//     publishDir "${params.outdir}/10.SV_germline", mode: 'copy'

//     input:
//     tuple val(sample_id), path(bam), path(bai)
//     val(params.ref)

//     output:
//     path "${sample_id}", emit: sv_germline

//     script:
//     """
//     # Activate the environment for Manta
//     source /opt/miniconda3/etc/profile.d/conda.sh
//     conda activate p2Manta

//     # Create a temporary directory in the work directory for Manta to run
//     mkdir -p 10.SV_germline
//     mkdir ${sample_id}

//     # Run Manta config and workflow
//     configManta.py --bam ${bam} --referenceFasta ${params.ref} --runDir ${sample_id}
//     ${sample_id}/runWorkflow.py -m local
//     """
// }

// process SV_germline {
//     tag "Structural Variants on ${sample_id}"
//     publishDir "${params.outdir}/10.SV_germlineD", mode: 'copy'

//     input:
//     tuple val(sample_id), path(bam), path(bai)
//     val(params.ref)

//     output:
//     //path "${sample_id}", emit: sv_germline
//     tuple val(sample_id), path("${sample_id}/results/variants/diploidSV.vcf.gz"), path("${annot_dir}/${sample_id}_SV_annotated.tsv"), path("${sample_id}_SV_Prioritized.xlsx"), emit: sv_germline
//     script:
//     def vcf_file = "${sample_id}/results/variants/diploidSV.vcf.gz"
//     def annot_date = new Date().format('yyyyMMdd')
//     def annot_dir = "${annot_date}_AnnotSV"
//     def annot_file = "${annot_dir}/${sample_id}_SV_annotated.tsv"
//     def prioritized_file = "${sample_id}_SV_Prioritized.xlsx"
//     """
//     # Activate the environment for Manta
//     source /opt/miniconda3/etc/profile.d/conda.sh
//     conda activate p2Manta

//     # Create a temporary directory in the work directory for Manta to run
//     mkdir ${sample_id}

//     # Run Manta config and workflow
//     configManta.py --bam ${bam} --referenceFasta ${params.ref} --runDir ${sample_id}
//     ${sample_id}/runWorkflow.py -m local

//     # Deactivate Python 2 environment
//     conda deactivate

//     # Verify Manta output
//     if [ ! -f "${vcf_file}" ]; then
//         echo "Error: Manta VCF file was not created"
//         exit 1
//     fi

//     # Annotation SV
//     /usr/src/app/AnnotSV/bin/./AnnotSV -SVinputFile ${vcf_file} -outputFile ${annot_file}

//     # Verify AnnotSV output
//     if [ ! -f "${annot_file}" ]; then
//         echo "Error: AnnotSV output file was not created"
//         exit 1
//     fi

//     ## prioritization
//     python3 /usr/src/app/ref38/Validation_script/Prior_SV_main.py --annotSR ${annot_file} --output ${prioritized_file}

//     # Verify prioritization output
//     if [ ! -f "${prioritized_file}" ]; then
//         echo "Error: Prioritized output file was not created"
//         exit 1
//     fi
//     """
// }
process SV_germline {
    tag "Structural Variants on ${sample_id}"
    publishDir "${params.outdir}/10.SV_germlineD", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.ref)

    output:
    tuple val(sample_id), path("${sample_id}/results/variants/diploidSV.vcf.gz"), path("*/${sample_id}_SV_annotated.tsv"), path("${sample_id}_SV_Prioritized.xlsx"), emit: sv_germline

    script:
    """
    # Activate the environment for Manta
    source /opt/miniconda3/etc/profile.d/conda.sh
    conda activate p2Manta

    # Create directory
    mkdir ${sample_id}

    # Run Manta
    configManta.py --bam ${bam} --referenceFasta ${params.ref} --runDir ${sample_id}
    ${sample_id}/runWorkflow.py -m local

    conda deactivate

    # Verify VCF
    if [ ! -f "${sample_id}/results/variants/diploidSV.vcf.gz" ]; then
        echo "Error: Manta VCF file was not created"
        exit 1
    fi

    # Define annotation output folder using current date
    annot_dir=\$(date +%Y%m%d)_AnnotSV
    annot_file=\${annot_dir}/${sample_id}_SV_annotated.tsv

    # Annotate SVs
    /usr/src/app/AnnotSV/bin/AnnotSV -SVinputFile ${sample_id}/results/variants/diploidSV.vcf.gz -outputFile \${annot_file}

    if [ ! -f "\${annot_file}" ]; then
        echo "Error: AnnotSV output file was not created"
        exit 1
    fi

    # Prioritization
    python3 /usr/src/app/ref38/Validation_script/Prior_SV_main.py --annotSR \${annot_file} --output ${sample_id}_SV_Prioritized.xlsx

    if [ ! -f "${sample_id}_SV_Prioritized.xlsx" ]; then
        echo "Error: Prioritized output file was not created"
        exit 1
    fi
    """
}

