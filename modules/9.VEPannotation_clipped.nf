process annotVEP {
    tag "VEP on ${sample_id}"
    publishDir "${params.outdir}/7.annotVEPP", mode: 'copy'
    cpus 16

    input:
    tuple val(sample_id), path(vcf)
    val(params.ref)
    val(params.cachee)
    val(params.dirplugin)
    val(params.dbNSFP)
    val(params.loftool)
    val(params.CADDsnv)
    val(params.CADDindel)
    val(params.dbscSNV)
    //val(params.DisGeNET)
    val(params.assembly)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_VEP.txt"), emit: annotVEP_vcfs

    script:
    """
    /home/ubuntu/VEP/ensembl-vep/ensembl-vep/./vep --biotype --buffer_size 50000 --cache --check_existing --database \
    --assembly ${params.assembly} \
    --dir ${params.cachee} \
    --dir_plugins ${params.dirplugin} \
    --fasta_dir ${params.ref} --force --fork 40 \
    --input_file ${vcf} --mane \
    --output_file ${sample_id}_filtered_PASS_VEP.txt --tab \
    --plugin dbNSFP,${params.dbNSFP},aapos,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,HGVSc_ANNOVAR,HGVSp_ANNOVAR,VEP_canonical,cds_strand,refcodon,codonpos \
    --quiet --safe --stats_text --symbol
    """
}
//--plugin DisGeNET,file=${params.DisGeNET} \ as not availble for hg19
//--offline 
//  --assembly GRCh37 ${params.assembly}
