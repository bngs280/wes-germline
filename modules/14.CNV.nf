// process CNV_germline {
//     tag "Germline Copy Number on ${sample_id}"
//     publishDir "${params.outdir}/11.CNV_germline", mode: 'copy'
//     cpus 20

//     input:
//     tuple val(sample_id), path(bam), path(bai)
//     val(params.cnv_reference)

//     output:
//     tuple val(sample_id), path("${sample_id}"), emit: cnv_Germline

//     script:
//     """
//     echo "Current PATH: \$PATH"
//     echo "CNVkit location: \$(which cnvkit.py)"
//     source /opt/miniconda3/etc/profile.d/conda.sh
//     conda activate cnvkit_env
//     conda install pandas=1.5.3
//     mkdir -p ${params.outdir}
//     cnvkit.py batch ${bam} -r ${params.cnv_reference} -p ${task.cpus} -d ${sample_id}
//     cnvkit.py export vcf ${sample_id}/${sample_id}_sorted_md.call.cns -o ${sample_id}/${sample_id}_cnv.vcf

//     echo "CNVkit processing completed for $sample_id"
//     """
// }

// process CNV_germline {
//     tag "Germline Copy Number on ${sample_id}"
//     publishDir "${params.outdir}/11.CNV_germlineD", mode: 'copy'
//     cpus 20

//     input:
//     tuple val(sample_id), path(bam), path(bai)
//     val(params.cnv_baseline)

//     output:
//     //tuple val(sample_id), path("${sample_id}_CNV_ranked.csv"), emit: cnv_Germline
//     tuple val(sample_id), path("${sample_id}_CNV_ranked.csv"), path("${sample_id}_CNV_ranked.vcf"), path("${annot_dir}/${sample_id}_CNV_ranked.annotated.tsv"), path("${sample_id}_CNV_Prioritized.xlsx"), emit: cnv_Germline

//     script:
//     def annot_date = new Date().format('yyyyMMdd')
//     def annot_dir = "${annot_date}_AnnotSV"
//     def csv_file = "${sample_id}_CNV_ranked.csv"
//     def vcf_file = "${sample_id}_CNV_ranked.vcf"
//     def annot_file = "${annot_dir}/${sample_id}_CNV_ranked.annotated.tsv"
//     def prioritized_file = "${sample_id}_CNV_Prioritized.xlsx"
//     """
//     Rscript /usr/src/app/ref38/Validation_script/Exome_depth_script.R --test ${bam} --reference default --refdir ${params.cnv_baseline} --output ${csv_file}
//     # Run the Perl script to convert CSV to VCF
//     perl /usr/src/app/ref38/Validation_script/ed_vcf.pl ${csv_file}
    
//     # Verify VCF was created
//     if [ ! -f "${vcf_file}" ]; then
//         echo "Error: VCF file was not created by ed_vcf.pl"
//         exit 1
//     fi
    
//     # Annotation CNV
//     mkdir -p ${annot_dir}
//     /usr/src/app/AnnotSV/bin/./AnnotSV -SVinputFile ${vcf_file} -outputFile ${annot_file}

//     # CNV Prioritization
//     python3 /usr/src/app/ref38/Validation_script/priori_CNV.py --exomedepth ${csv_file} --annotSV ${annot_file} --output ${prioritized_file}
//     """
// }
process CNV_germline {
    tag "Germline Copy Number on ${sample_id}"
    publishDir "${params.outdir}/11.CNV_germlineD", mode: 'copy'
    cpus 20

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(params.cnv_baseline)

    output:
    tuple val(sample_id),
          path("${sample_id}_CNV_ranked.csv"),
          path("${sample_id}_CNV_ranked.vcf"),
          path("*/${sample_id}_CNV_ranked.annotated.tsv"),
          path("${sample_id}_CNV_Prioritized.xlsx"),
          emit: cnv_Germline

    script:
    def annot_date = new Date().format('yyyyMMdd')
    def annot_dir = "${annot_date}_AnnotSV"
    def csv_file = "${sample_id}_CNV_ranked.csv"
    def vcf_file = "${sample_id}_CNV_ranked.vcf"
    def annot_file = "${annot_dir}/${sample_id}_CNV_ranked.annotated.tsv"
    def prioritized_file = "${sample_id}_CNV_Prioritized.xlsx"

    """
    # Step 1: CNV calling using ExomeDepth
    Rscript /usr/src/app/ref38/Validation_script/Exome_depth_script.R \\
        --test ${bam} \\
        --reference default \\
        --refdir ${params.cnv_baseline} \\
        --output ${csv_file}

    # Step 2: Convert CSV to VCF
    perl /usr/src/app/ref38/Validation_script/ed_vcf.pl ${csv_file}

    # Check if VCF file was created
    if [ ! -f "${vcf_file}" ]; then
        echo "Error: VCF file was not created by ed_vcf.pl"
        exit 1
    fi

    # Step 3: Annotate CNVs using AnnotSV
    mkdir -p ${annot_dir}
    /usr/src/app/AnnotSV/bin/AnnotSV \\
        -SVinputFile ${vcf_file} \\
        -outputFile ${annot_file}

    # Step 4: Prioritize CNVs using in-house Python script
    python3 /usr/src/app/ref38/Validation_script/priori_CNV.py \\
        --exomedepth ${csv_file} \\
        --annotSV ${annot_file} \\
        --output ${prioritized_file}
    """
}
