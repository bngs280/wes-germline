process annotVEP {
    tag "VEP on ${sample_id}"
    publishDir "${params.outdir}/7.annotVEPP", mode: 'copy'
    cpus 80

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
    val(params.assembly)
    val(dosage_collins)
    val(spliceaiindel)
    val(spliceaisnv)
    val(pli)
    val(phenotypesfile)
    val(loeuf)
    val(clinpred)
    val(alphamissense)
    val(primateai)
    val(revel)
    val(clinvar)
    

    output:
    tuple val(sample_id), path("${sample_id}_filtered_PASS_VEP.txt"), emit: annotVEP_vcfs

    script:
    """
    /usr/src/app/ensembl-vep/./vep --af --appris --biotype --buffer_size 40000 --offline --cache --check_existing \
        --assembly ${params.assembly} \
        --dir ${params.cachee} \
        --dir_plugins ${params.dirplugin} \
        --distance 5000 \
        --fasta ${params.ref} \
        --force --fork 75 --mane \
        --input_file ${vcf} --output_file ${sample_id}_filtered_PASS_VEP.txt --tab --force_overwrite \
        --polyphen b --pubmed \
        --everything --total_length --numbers --canonical --hgvs \
        --plugin LoFtool,${params.loftool} \
        --plugin CADD,snv=${params.CADDsnv},indels=${params.CADDindel} \
        --plugin dbscSNV,${params.dbscSNV} \
        --plugin LOVD \
        --plugin DosageSensitivity,file=${params.dosage_collins} \
        --plugin SpliceRegion,Extended \
        --plugin SpliceAI,snv=${params.spliceaisnv},indel=${params.spliceaiindel} \
        --quiet --regulatory --safe --show_ref_allele --sift b --stats_text --symbol  --tsl --uploaded_allele \
        --plugin pLI,${params.pli} \
        --plugin Phenotypes,file=${params.phenotypesfile},include_types=Gene \
        --plugin dbNSFP,${params.dbNSFP},1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,1000Gp3_EUR_AF,1000Gp3_SAS_AF,ALFA_African_AF,ALFA_African_American_AF,ALFA_African_Others_AF,ALFA_Asian_AF,ALFA_East_Asian_AF,ALFA_European_AF,ALFA_Latin_American_1_AF,ALFA_Latin_American_2_AF,ALFA_Other_AF,ALFA_Other_Asian_AF,ALFA_South_Asian_AF,ALFA_Total_AF,AllofUs_AFR_AF,AllofUs_ALL_AF,AllofUs_AMR_AF,AllofUs_EAS_AF,AllofUs_EUR_AF,AllofUs_MID_AF,AllofUs_OTH_AF,AllofUs_POPMAX_AF,AllofUs_SAS_AF,Aloft_pred,Aloft_prob_Dominant,Aloft_prob_Recessive,Aloft_prob_Tolerant,AlphaMissense_pred,AlphaMissense_rankscore,AlphaMissense_score,BayesDel_addAF_pred,BayesDel_addAF_rankscore,BayesDel_addAF_score,BayesDel_noAF_pred,BayesDel_noAF_rankscore,BayesDel_noAF_score,ClinPred_pred,ClinPred_rankscore,ClinPred_score,DANN_rankscore,DANN_score,DEOGEN2_pred,DEOGEN2_rankscore,DEOGEN2_score,ESM1b_converted_rankscore,ESM1b_pred,ESM1b_score,Eigen-PC-phred_coding,Eigen-PC-raw_coding,Eigen-PC-raw_coding_rankscore,Eigen-phred_coding,Eigen-raw_coding,Eigen-raw_coding_rankscore,Ensembl_geneid,Ensembl_proteinid,Ensembl_transcriptid,GENCODE_basic,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,GERP_91_mammals,GERP_91_mammals_rankscore,HGVSc_VEP,HGVSc_snpEff,HGVSp_VEP,HGVSp_snpEff,Interpro_domain,LIST-S2_pred,LIST-S2_rankscore,LIST-S2_score,M-CAP_pred,M-CAP_rankscore,M-CAP_score,MANE,MPC_rankscore,MPC_score,MVP_rankscore,MVP_score,MetaLR_pred,MetaLR_rankscore,MetaLR_score,MetaRNN_pred,MetaRNN_rankscore,MetaRNN_score,MetaSVM_pred,MetaSVM_rankscore,MetaSVM_score,MutFormer_rankscore,MutFormer_score,MutPred_AAchange,MutPred_Top5features,MutPred_protID,MutPred_rankscore,MutPred_score,MutScore_rankscore,MutScore_score,MutationAssessor_pred,MutationAssessor_rankscore,MutationAssessor_score,MutationTaster_model,MutationTaster_pred,MutationTaster_rankscore,MutationTaster_score,MutationTaster_trees_benign,MutationTaster_trees_deleterious,PHACTboost_rankscore,PHACTboost_score,PROVEAN_converted_rankscore,PROVEAN_pred,PROVEAN_score,Polyphen2_HDIV_pred,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_score,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_score,PrimateAI_pred,PrimateAI_rankscore,PrimateAI_score,REVEL_rankscore,REVEL_score,RegeneronME_AFR_AF,RegeneronME_ALL_AF,RegeneronME_AMI_AF,RegeneronME_ASH_AF,RegeneronME_BI_AF,RegeneronME_C_EUR_AF,RegeneronME_EAS_AF,RegeneronME_EUR_AF,RegeneronME_E_AFR_AF,RegeneronME_E_ASIA_AF,RegeneronME_FIN_AF,RegeneronME_GBR_AF,RegeneronME_GHA_AF,RegeneronME_IAM_AF,RegeneronME_IND_AF,RegeneronME_IRE_AF,RegeneronME_ITA_AF,RegeneronME_MEA_AF,RegeneronME_MEX_AF,RegeneronME_NIN_AF,RegeneronME_N_EUR_AF,RegeneronME_PAK_AF,RegeneronME_SAS_AF,RegeneronME_SE_ASIA_AF,RegeneronME_SPA_AF,RegeneronME_S_AFR_AF,RegeneronME_S_EUR_AF,RegeneronME_W_AFR_AF,Reliability_index,SIFT4G_converted_rankscore,SIFT4G_pred,SIFT4G_score,SIFT_converted_rankscore,SIFT_pred,SIFT_score,TOPMed_frz8_AF,TSL,Uniprot_acc,Uniprot_entry,VARITY_ER_LOO_rankscore,VARITY_ER_LOO_score,VARITY_ER_rankscore,VARITY_ER_score,VARITY_R_LOO_rankscore,VARITY_R_LOO_score,VARITY_R_rankscore,VARITY_R_score,VEP_canonical,VEST4_rankscore,VEST4_score,aapos,bStatistic,bStatistic_converted_rankscore,cds_strand,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,clinvar_clnsig,clinvar_hgvs,clinvar_id,clinvar_review,clinvar_trait,clinvar_var_source,codon_degeneracy,codonpos,dbNSFP_POPMAX_AF,fathmm-XF_coding_pred,fathmm-XF_coding_rankscore,fathmm-XF_coding_score,gMVP_rankscore,gMVP_score,genename,gnomAD4.1_joint_AF,gnomAD4.1_joint_AFR_AF,gnomAD4.1_joint_AMI_AF,gnomAD4.1_joint_AMR_AF,gnomAD4.1_joint_ASJ_AF,gnomAD4.1_joint_EAS_AF,gnomAD4.1_joint_FIN_AF,gnomAD4.1_joint_MID_AF,gnomAD4.1_joint_NFE_AF,gnomAD4.1_joint_POPMAX_AF,gnomAD4.1_joint_SAS_AF,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,phastCons17way_primate,phastCons17way_primate_rankscore,phastCons470way_mammalian,phastCons470way_mammalian_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phyloP17way_primate,phyloP17way_primate_rankscore,phyloP470way_mammalian,phyloP470way_mammalian_rankscore,refcodon,rs_dbSNP \
        --pick --per_gene --filter_common \
        --plugin LOEUF,file=${params.loeuf},match_by=gene,match_by=gene \
        --plugin AlphaMissense,file=${params.alphamissense},cols=all \
        --plugin PrimateAI,${params.primateai} \
        --plugin ClinPred,file=${params.clinpred} \
        --plugin REVEL,file=${params.revel} \
        --plugin NearestExonJB,max_range=50000 \
        --plugin NMD \
        --plugin SpliceAI,snv=${params.spliceaisnv},indel=${params.spliceaiindel} \
        --custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNSIGCONF%CLNREVSTAT%CLNDN%CLNVI%CLNDISDB%ONCDN%ONCDNINCL%ONCDISDB%SCIDN%SCIDNINCL%SCIDISDB
    """
}
//--plugin DisGeNET,file=${params.DisGeNET} \ as not availble for hg19
//--offline 
//  --assembly GRCh37 ${params.assembly}
// /usr/src/app/ensembl-vep/./vep --biotype --buffer_size 50000 --cache --check_existing --database --offline \
//     --assembly ${params.assembly} \
//     --dir ${params.cachee} \
//     --dir_plugins ${params.dirplugin} \
//     --fasta_dir ${params.ref} --force --fork 50 \
//     --input_file ${vcf} --mane \
//     --output_file ${sample_id}_filtered_PASS_VEP.txt --tab \
//     --everything --total_length --numbers --canonical --hgvs \
//     --plugin dbNSFP,${params.dbNSFP},aapos,genename,Ensembl_geneid,Ensembl_transcriptid,Ensembl_proteinid,HGVSc_VEP,HGVSp_VEP,APPRIS,GENCODE_basic,refcodon,Ancestral_allele,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,VEST4_score,VEST4_rankscore,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,MetaRNN_score,MetaRNN_rankscore,MetaRNN_pred,M-CAP_score,M-CAP_rankscore,M-CAP_pred,REVEL_score,REVEL_rankscore,MutPred_score,MutPred_rankscore,MutPred_protID,MutPred_AAchange,MutPred_Top5features,MVP_score,MVP_rankscore,MPC_score,MPC_rankscore,DEOGEN2_score,DEOGEN2_rankscore,DEOGEN2_pred,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred,ClinPred_score,ClinPred_rankscore,ClinPred_pred,LIST-S2_score,LIST-S2_rankscore,LIST-S2_pred,VARITY_R_score,VARITY_R_rankscore,VARITY_ER_score,VARITY_ER_rankscore,VARITY_R_LOO_score,VARITY_R_LOO_rankscore,VARITY_ER_LOO_score,VARITY_ER_LOO_rankscore,AlphaMissense_score,AlphaMissense_rankscore,AlphaMissense_pred,Aloft_Fraction_transcripts_affected,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,DANN_score,DANN_rankscore,fathmm-XF_coding_score,fathmm-XF_coding_rankscore,fathmm-XF_coding_pred,Eigen-raw_coding,Eigen-raw_coding_rankscore,Eigen-phred_coding,Eigen-PC-raw_coding,Eigen-PC-raw_coding_rankscore,Eigen-PC-phred_coding,GenoCanyon_score,GenoCanyon_rankscore,integrated_fitCons_score,integrated_fitCons_rankscore,integrated_confidence_value,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,1000Gp3_AC,1000Gp3_AF,1000Gp3_SAS_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,1000Gp3_EUR_AF,UK10K_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_SAS_AC,ExAC_SAS_AF,gnomAD_exomes_flag,gnomAD_exomes_AF,gnomAD_exomes_SAS_AF,gnomAD_exomes_AFR_AF,gnomAD_exomes_AMR_AF,gnomAD_exomes_EAS_AF,gnomAD_exomes_FIN_AF,gnomAD_exomes_NFE_AF,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_var_source,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain \
//     --pubmed \
//     --plugin DosageSensitivity,file=${params.dosage_collins} \
//     --plugin SpliceRegion,Extended \
//     --plugin Phenotypes,file=${params.phenotypesfile},include_types=Gene \
//     --plugin LoFtool,${params.loftool} \
//     --plugin CADD,snv=${params.CADDsnv},indels=${params.CADDindel} \
//     --plugin dbscSNV,${params.dbscSNV} \
//     --plugin LOEUF,file=${params.loeuf},match_by=gene,match_by=gene \
//     --plugin AlphaMissense,file=${params.alphamissense},cols=all \
//     --plugin PrimateAI,${params.primateai} \
//     --plugin ClinPred,file=${params.clinpred} \
//     --plugin REVEL,file=${params.revel} \
//     --plugin NearestExonJB,max_range=50000 \
//     --plugin NMD \
//     --plugin pLI,${params.pli} \
//     --plugin SpliceAI,snv=${params.spliceaisnv},indel=${params.spliceaiindel} \
//     --custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
//     --plugin LOVD --pick --per_gene --filter_common --quiet --safe --stats_text --symbol --transcript_version --sift b --polyphen b --show_ref_allele --transcript_version --symbol
