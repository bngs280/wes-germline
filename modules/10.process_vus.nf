process vusVEP {
    tag "VEP on ${sample_id}"
    publishDir "${params.outdir}/9.vus_vep", mode: 'copy'
    cpus 15

    input:
    tuple val(sample_id), path(vcf)
    val(params.ref)
    val(params.cachee)
    val(params.dirplugin)
    val(params.dbnsfp4)
    val(params.loftool)
    val(params.CADDsnv)
    val(params.CADDindel)
    val(params.dbscSNV)
    val(params.assembly)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_VUS_VEP.txt"), emit: vus_vep

    script:
    """
    /usr/src/app/ensembl-vep/vep \
    --af --appris --buffer_size 50000 --cache --check_existing --distance 5000 \
    --assembly ${params.assembly} \
    --dir ${params.cachee} \
    --dir_plugins ${params.dirplugin} \
    --mane --pick \
    --fasta ${params.ref} \
    --force --fork 15 \
    --input_file ${vcf} --format vcf \
    --output_file ${sample_id}_filtered_VUS_VEP.txt --tab \
    --plugin MaxEntScan,/usr/src/app/vepDB/DBs/fordownload \
    --plugin Blosum62 \
    --plugin dbNSFP,${params.dbnsfp4},codon_degeneracy,Eigen-phred_coding,integrated_fitCons_score,GERP++_RS,phyloP100way_vertebrate \
    --plugin dbscSNV,${params.dbscSNV} \
    --plugin LoFtool,${params.loftool} \
    --plugin CADD,snv=${params.CADDsnv} \
    --refseq --quiet --regulatory --sift s --species homo_sapiens --symbol --transcript_version --tsl --safe --show_ref_allele --stats_text --uploaded_allele
    """
}

