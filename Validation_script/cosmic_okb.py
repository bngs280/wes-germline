import pandas as pd
import argparse
import os
from xlsxwriter import Workbook
def merge_tsv(tsv1_path, tsv2_path, output_path):
    # Predefined columns for matching, extracting, and removing
    tsv1_match_cols = ['Chr', 'Position', 'Ref', 'Alt']
    tsv2_match_cols = ['Chromosome', 'Start_Position', 'Reference_Allele', 'Allele']
    tsv2_extract_cols = [
        'Hugo_Symbol', 'Variant_Classification', 'dbSNP_RS', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Exon_Number', 
        'all_effects', 'RefSeq', 'IS-A-HOTSPOT', 'ANNOTATED', 'GENE_IN_ONCOKB', 'VARIANT_IN_ONCOKB', 'MUTATION_EFFECT', 
        'ONCOGENIC', 'LEVEL_1', 'LEVEL_2', 'LEVEL_3A', 'LEVEL_3B', 'LEVEL_4', 'LEVEL_R1', 'LEVEL_R2', 'HIGHEST_LEVEL', 
        'HIGHEST_SENSITIVE_LEVEL', 'HIGHEST_RESISTANCE_LEVEL', 'TX_CITATIONS', 'LEVEL_Dx1', 'LEVEL_Dx2', 'LEVEL_Dx3', 
        'HIGHEST_DX_LEVEL', 'DX_CITATIONS', 'LEVEL_Px1', 'LEVEL_Px2', 'LEVEL_Px3', 'HIGHEST_PX_LEVEL', 'PX_CITATIONS', 
        'GENE_SUMMARY', 'VARIANT_SUMMARY', 'TUMOR_TYPE_SUMMARY', 'DIAGNOSTIC_SUMMARY', 'PROGNOSTIC_SUMMARY', 
        'MUTATION_EFFECT_DESCRIPTION'
    ]
    tsv1_remove_cols = [
        '#Uploaded_variation', 'Location', 'FLAGS', 'MANE_PLUS_CLINICAL', 'SIFT', 'PolyPhen', 'CLIN_SIG', '1000Gp3_AC',
        '1000Gp3_AFR_AC', '1000Gp3_AMR_AC', '1000Gp3_EAS_AC', '1000Gp3_EUR_AC', '1000Gp3_SAS_AC', 'ALFA_African_AC',
        'ALFA_African_AN', 'ALFA_African_American_AC', 'ALFA_African_American_AN', 'ALFA_African_Others_AC',
        'ALFA_African_Others_AN', 'ALFA_Asian_AC', 'ALFA_Asian_AN', 'ALFA_East_Asian_AC', 'ALFA_East_Asian_AN',
        'ALFA_European_AC', 'ALFA_European_AN', 'ALFA_Latin_American_1_AC', 'ALFA_Latin_American_1_AF',
        'ALFA_Latin_American_1_AN', 'ALFA_Latin_American_2_AC', 'ALFA_Latin_American_2_AF', 'ALFA_Latin_American_2_AN',
        'ALFA_Other_AC', 'ALFA_Other_AN', 'ALFA_Other_Asian_AC', 'ALFA_Other_Asian_AN', 'ALFA_South_Asian_AC',
        'ALFA_South_Asian_AN', 'ALFA_Total_AC', 'ALFA_Total_AN', 'Aloft_Confidence', 'Aloft_Fraction_transcripts_affected',
        'Aloft_prob_Dominant', 'Aloft_prob_Recessive', 'Aloft_prob_Tolerant', 'AlphaMissense_rankscore',
        'AlphaMissense_score', 'BayesDel_addAF_rankscore', 'BayesDel_addAF_score', 'CADD_phred', 'CADD_phred_hg19',
        'CADD_raw', 'CADD_raw_hg19', 'CADD_raw_rankscore', 'CADD_raw_rankscore_hg19', 'ClinPred_rankscore', 'ClinPred_score',
        'DANN_rankscore', 'DEOGEN2_rankscore', 'DEOGEN2_score', 'ESM1b_rankscore', 'ESM1b_score', 'ESP6500_AA_AC',
        'ESP6500_AA_AF', 'ESP6500_EA_AC', 'ESP6500_EA_AF', 'EVE_Class10_pred', 'EVE_Class20_pred', 'EVE_Class25_pred',
        'EVE_Class30_pred', 'EVE_Class40_pred', 'EVE_Class50_pred', 'EVE_Class60_pred', 'EVE_Class70_pred', 'EVE_Class75_pred',
        'EVE_Class80_pred', 'EVE_Class90_pred', 'EVE_rankscore', 'EVE_score', 'Eigen-PC-phred_coding', 'Eigen-PC-raw_coding',
        'Eigen-PC-raw_coding_rankscore', 'Eigen-phred_coding', 'Eigen-raw_coding', 'Eigen-raw_coding_rankscore', 'Ensembl_geneid',
        'Ensembl_proteinid', 'Ensembl_transcriptid', 'ExAC_AC', 'ExAC_AFR_AC', 'ExAC_AMR_AC', 'ExAC_Adj_AC', 'ExAC_Adj_AF', 
        'ExAC_EAS_AC', 'ExAC_FIN_AC', 'ExAC_NFE_AC', 'ExAC_SAS_AC', 'ExAC_nonTCGA_AC', 'ExAC_nonTCGA_AF', 'ExAC_nonTCGA_AFR_AC', 
        'ExAC_nonTCGA_AFR_AF', 'ExAC_nonTCGA_AMR_AC', 'ExAC_nonTCGA_AMR_AF', 'ExAC_nonTCGA_Adj_AC', 'ExAC_nonTCGA_Adj_AF', 
        'ExAC_nonTCGA_EAS_AC', 'ExAC_nonTCGA_EAS_AF', 'ExAC_nonTCGA_FIN_AC', 'ExAC_nonTCGA_FIN_AF', 'ExAC_nonTCGA_NFE_AC', 
        'ExAC_nonTCGA_NFE_AF', 'ExAC_nonTCGA_SAS_AC', 'ExAC_nonTCGA_SAS_AF', 'ExAC_nonpsych_AC', 'ExAC_nonpsych_AF', 
        'ExAC_nonpsych_AFR_AC', 'ExAC_nonpsych_AFR_AF', 'ExAC_nonpsych_AMR_AC', 'ExAC_nonpsych_AMR_AF', 'ExAC_nonpsych_Adj_AC', 
        'ExAC_nonpsych_Adj_AF', 'ExAC_nonpsych_EAS_AC', 'ExAC_nonpsych_EAS_AF', 'ExAC_nonpsych_FIN_AC', 'ExAC_nonpsych_FIN_AF', 
        'ExAC_nonpsych_NFE_AC', 'ExAC_nonpsych_NFE_AF', 'ExAC_nonpsych_SAS_AC', 'ExAC_nonpsych_SAS_AF', 'FATHMM_converted_rankscore', 
        'FATHMM_score', 'GERP++_NR', 'GERP++_RS', 'GERP++_RS_rankscore', 'GM12878_confidence_value', 'GM12878_fitCons_rankscore', 
        'GM12878_fitCons_score', 'GTEx_V8_eQTL_gene', 'GTEx_V8_eQTL_tissue', 'GTEx_V8_sQTL_gene', 'GTEx_V8_sQTL_tissue', 
        'GenoCanyon_rankscore', 'GenoCanyon_score', 'H1-hESC_confidence_value', 'H1-hESC_fitCons_rankscore', 'H1-hESC_fitCons_score', 
        'HGVSc_VEP', 'HGVSc_snpEff', 'HGVSp_VEP', 'HGVSp_snpEff', 'HUVEC_confidence_value', 'HUVEC_fitCons_rankscore', 'HUVEC_fitCons_score', 
        'Interpro_domain', 'LINSIGHT', 'LINSIGHT_rankscore', 'LIST-S2_rankscore', 'LIST-S2_score', 'LRT_Omega', 'LRT_converted_rankscore', 
        'LRT_score', 'M-CAP_rankscore', 'M-CAP_score', 'MPC_rankscore', 'MVP_rankscore', 'MetaLR_rankscore', 'MetaLR_score', 'MetaRNN_rankscore', 
        'MetaRNN_score', 'MetaSVM_rankscore', 'MetaSVM_score', 'MutPred_rankscore', 'MutationAssessor_rankscore', 'MutationAssessor_score', 
        'MutationTaster_converted_rankscore', 'MutationTaster_score', 'PROVEAN_converted_rankscore', 'PROVEAN_score', 
        'Polyphen2_HDIV_rankscore', 'Polyphen2_HDIV_score', 'Polyphen2_HVAR_rankscore', 'Polyphen2_HVAR_score', 'PrimateAI_rankscore', 
        'PrimateAI_score', 'REVEL_rankscore', 'SIFT4G_converted_rankscore', 'SIFT4G_score', 'SIFT_converted_rankscore', 'SIFT_score', 
        'SiPhy_29way_logOdds', 'SiPhy_29way_logOdds_rankscore', 'SiPhy_29way_pi', 'TSL', 'UK10K_AC', 'Uniprot_acc', 'Uniprot_entry', 
        'VARITY_ER_LOO_rankscore', 'VARITY_ER_LOO_score', 'VARITY_ER_rankscore', 'VARITY_ER_score', 'VARITY_R_LOO_rankscore', 'VARITY_R_LOO_score', 
        'VARITY_R_rankscore', 'VARITY_R_score', 'VEP_canonical', 'VEST4_rankscore', 'VEST4_score', 'aapos', 'bStatistic', 'bStatistic_converted_rankscore', 
        'cds_strand', 'clinvar_MedGen_id', 'clinvar_Orphanet_id', 'clinvar_id', 'clinvar_review', 'clinvar_var_source', 'fathmm-MKL_coding_group', 
        'fathmm-MKL_coding_rankscore', 'fathmm-MKL_coding_score', 'fathmm-XF_coding_rankscore', 'fathmm-XF_coding_score', 'gMVP_rankscore', 
        'gMVP_score', 'genename', 'gnomAD_exomes_AC', 'gnomAD_exomes_AFR_AC', 'gnomAD_exomes_AFR_AN', 'gnomAD_exomes_AFR_nhomalt', 
        'gnomAD_exomes_AMR_AC', 'gnomAD_exomes_AMR_AN', 	'gnomAD_exomes_AMR_nhomalt', 'gnomAD_exomes_AN', 'gnomAD_exomes_ASJ_AC', 
        'gnomAD_exomes_ASJ_AN', 'gnomAD_exomes_ASJ_nhomalt', 'gnomAD_exomes_EAS_AC', 'gnomAD_exomes_EAS_AN', 'gnomAD_exomes_EAS_nhomalt', 'gnomAD_exomes_FIN_AC', 
        'gnomAD_exomes_FIN_AN', 'gnomAD_exomes_FIN_nhomalt', 'gnomAD_exomes_MID_AC', 'gnomAD_exomes_MID_AN', 'gnomAD_exomes_MID_nhomalt', 
        'gnomAD_exomes_NFE_AC', 'gnomAD_exomes_NFE_AN', 'gnomAD_exomes_NFE_nhomalt', 'gnomAD_exomes_POPMAX_AC', 'gnomAD_exomes_POPMAX_AN', 
        'gnomAD_exomes_SAS_AC', 'gnomAD_exomes_SAS_AN', 'gnomAD_exomes_SAS_nhomalt', 'gnomAD_exomes_flag', 'gnomAD_exomes_non_ukb_AC', 'gnomAD_exomes_non_ukb_AF', 
        'gnomAD_exomes_non_ukb_AFR_AC', 'gnomAD_exomes_non_ukb_AFR_AF', 'gnomAD_exomes_non_ukb_AFR_AN', 'gnomAD_exomes_non_ukb_AFR_nhomalt', 
        'gnomAD_exomes_non_ukb_AMR_AC', 'gnomAD_exomes_non_ukb_AMR_AF', 'gnomAD_exomes_non_ukb_AMR_AN', 'gnomAD_exomes_non_ukb_AMR_nhomalt', 
        'gnomAD_exomes_non_ukb_AN', 'gnomAD_exomes_non_ukb_ASJ_AC', 'gnomAD_exomes_non_ukb_ASJ_AF', 'gnomAD_exomes_non_ukb_ASJ_AN', 
        'gnomAD_exomes_non_ukb_ASJ_nhomalt', 'gnomAD_exomes_non_ukb_EAS_AC', 'gnomAD_exomes_non_ukb_EAS_AF', 'gnomAD_exomes_non_ukb_EAS_AN', 
        'gnomAD_exomes_non_ukb_EAS_nhomalt', 'gnomAD_exomes_non_ukb_FIN_AC', 'gnomAD_exomes_non_ukb_FIN_AF', 'gnomAD_exomes_non_ukb_FIN_AN', 
        'gnomAD_exomes_non_ukb_FIN_nhomalt', 'gnomAD_exomes_non_ukb_MID_AC', 'gnomAD_exomes_non_ukb_MID_AF', 'gnomAD_exomes_non_ukb_MID_AN', 
        'gnomAD_exomes_non_ukb_MID_nhomalt', 'gnomAD_exomes_non_ukb_NFE_AC', 'gnomAD_exomes_non_ukb_NFE_AF', 'gnomAD_exomes_non_ukb_NFE_AN', 
        'gnomAD_exomes_non_ukb_NFE_nhomalt', 'gnomAD_exomes_non_ukb_SAS_AC', 'gnomAD_exomes_non_ukb_SAS_AF', 'gnomAD_exomes_non_ukb_SAS_AN', 
        'gnomAD_exomes_non_ukb_SAS_nhomalt', 'gnomAD_exomes_non_ukb_nhomalt', 'integrated_confidence_value', 'integrated_fitCons_rankscore', 
        'phastCons100way_vertebrate', 'phastCons100way_vertebrate_rankscore', 'phastCons17way_primate', 'phastCons17way_primate_rankscore', 
        'phastCons470way_mammalian_rankscore', 'phyloP100way_vertebrate', 'phyloP100way_vertebrate_rankscore', 'phyloP17way_primate', 
        'phyloP17way_primate_rankscore', 'phyloP470way_mammalian_rankscore', 'LOVD', 'Classification', 'Allele', 'Gene', 'Feature_type', 
        'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'DISTANCE', 'HGNC_ID', 'BIOTYPE', 'SOMATIC', 'PHENO'
    ]
    # Extract the initial Sample ID from filenames
    tsv1_sample_id = os.path.basename(tsv1_path).split("_")[0]
    tsv2_sample_id = os.path.basename(tsv2_path).split("_")[0]

    # Check if the Sample IDs match
    if tsv1_sample_id != tsv2_sample_id:
        print(f"Sample IDs do not match: {tsv1_sample_id} vs {tsv2_sample_id}")
        return
    # Read TSV files
    tsv1 = pd.read_csv(tsv1_path, sep='\t')
    tsv2 = pd.read_csv(tsv2_path, sep='\t')
    
    # Drop unwanted columns from tsv1
    tsv1_cleaned = tsv1.drop(columns=tsv1_remove_cols, errors='ignore')
    
    # Merge tsv1 and tsv2 based on matching columns
    merged_df = pd.merge(tsv1_cleaned, tsv2[tsv2_extract_cols + tsv2_match_cols], 
                         left_on=tsv1_match_cols, right_on=tsv2_match_cols, how='left')
    
   
    # Drop unwanted columns
    merged_df = merged_df.drop(columns=['Chromosome', 'Start_Position', 'Reference_Allele', 'Allele'])

    # Rename specific columns
    rename_dict = {
        'SYMBOL': 'Gene',
        'Variant_Classification': 'Function',
        'ada_score': 'dbscSNV_ada',
        'rf_score': 'dbscSNV_rf'
    }
    merged_df = merged_df.rename(columns=rename_dict)

    # Reorder columns
    columns_order = [
        'Chr', 'Position', 'Ref', 'Alt', 'Gene', 'Hugo_Symbol', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Exon_Number',
        'all_effects', 'RefSeq', 'Consequence', 'Function', 'IMPACT', 'dbSNP_RS', 'AMP_ACMG_Classification', 'COSMIC_ID',
        'COSMIC_STATUS', 'HGVSc_ANNOVAR', 'HGVSp_ANNOVAR', 'Existing_variation', 'STRAND', '1000Gp3_AF', '1000Gp3_AFR_AF',
        '1000Gp3_AMR_AF', '1000Gp3_EAS_AF', '1000Gp3_EUR_AF', '1000Gp3_SAS_AF', 'ALFA_African_AF', 'ALFA_African_American_AF',
        'ALFA_African_Others_AF', 'ALFA_Asian_AF', 'ALFA_East_Asian_AF', 'ALFA_European_AF', 'ALFA_Other_AF', 'ALFA_Other_Asian_AF',
        'ALFA_South_Asian_AF', 'ALFA_Total_AF', 'ExAC_AF', 'ExAC_AFR_AF', 'ExAC_AMR_AF', 'ExAC_EAS_AF', 'ExAC_FIN_AF',
        'ExAC_NFE_AF', 'ExAC_SAS_AF', 'UK10K_AF', 'gnomAD_exomes_AF', 'gnomAD_exomes_AFR_AF', 'gnomAD_exomes_AMR_AF',
        'gnomAD_exomes_ASJ_AF', 'gnomAD_exomes_EAS_AF', 'gnomAD_exomes_FIN_AF', 'gnomAD_exomes_MID_AF', 'gnomAD_exomes_NFE_AF',
        'gnomAD_exomes_POPMAX_AF', 'gnomAD_exomes_SAS_AF', 'LIST-S2_pred', 'LRT_pred', 'M-CAP_pred', 'MPC_score', 'MVP_score',
        'MetaLR_pred', 'MetaRNN_pred', 'MetaSVM_pred', 'MutPred_Top5features', 'MutPred_score', 'MutationAssessor_pred',
        'MutationTaster_pred', 'PROVEAN_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'PrimateAI_pred', 'REVEL_score',
        'SIFT4G_pred', 'SIFT_pred', 'Aloft_pred', 'AlphaMissense_pred', 'BayesDel_addAF_pred', 'ClinPred_pred', 'DANN_score',
        'DEOGEN2_pred', 'ESM1b_pred', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred', 'FATHMM_pred', 'integrated_fitCons_score',
        'phastCons470way_mammalian', 'phyloP470way_mammalian', 'LoFtool', 'CADD_PHRED', 'CADD_RAW', 'dbscSNV_ada', 'dbscSNV_rf',
        'Feature', 'MANE', 'MANE_SELECT', 'PUBMED', 'clinvar_OMIM_id', 'clinvar_clnsig', 'clinvar_hgvs', 'clinvar_trait', 
        'SYMBOL_SOURCE', 'GERMQ', 'MBQ', 'MMQ', 'TLOD', 'GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'FAD', 'SB', 'IS-A-HOTSPOT',
        'ANNOTATED', 'GENE_IN_ONCOKB', 'VARIANT_IN_ONCOKB', 'MUTATION_EFFECT', 'ONCOGENIC', 'LEVEL_1', 'LEVEL_2', 'LEVEL_3A',
        'LEVEL_3B', 'LEVEL_4', 'LEVEL_R1', 'LEVEL_R2', 'HIGHEST_LEVEL', 'HIGHEST_SENSITIVE_LEVEL', 'HIGHEST_RESISTANCE_LEVEL',
        'TX_CITATIONS', 'LEVEL_Dx1', 'LEVEL_Dx2', 'LEVEL_Dx3', 'HIGHEST_DX_LEVEL', 'DX_CITATIONS', 'LEVEL_Px1', 'LEVEL_Px2',
        'LEVEL_Px3', 'HIGHEST_PX_LEVEL', 'PX_CITATIONS', 'GENE_SUMMARY', 'VARIANT_SUMMARY', 'TUMOR_TYPE_SUMMARY', 'DIAGNOSTIC_SUMMARY',
        'PROGNOSTIC_SUMMARY', 'MUTATION_EFFECT_DESCRIPTION'
    ]

    # Reordering columns based on the provided order
    merged_df = merged_df[columns_order]

    # Step 1: Pre-filter based on DP >= 10 or COSMIC_ID column
    pre_filtered_df = merged_df[(merged_df['DP'] >= 10) | (~merged_df['COSMIC_ID'].isna())]
    #main_filtered_df = merged_df.drop(pre_filtered_df.index)

    # Step 2: Filter based on ONCOGENIC containing 'Likely Oncogenic' or 'Oncogenic'
    oncogenic_filtered_df = merged_df[merged_df['ONCOGENIC'].str.contains('Likely Oncogenic|Oncogenic|Resistance|Likely Neutral', na=False, regex=True)]
    #main_filtered_df = merged_df.drop(final_filtered_df.index)

    # Step 3: filter based on DP >= 10 or COSMIC_ID column
    actionable_filtered_df = merged_df[(~merged_df['HIGHEST_LEVEL'].isna()) | (~merged_df['HIGHEST_SENSITIVE_LEVEL'].isna()) | (~merged_df['HIGHEST_RESISTANCE_LEVEL'].isna()) | (~merged_df['HIGHEST_DX_LEVEL'].isna()) | (~merged_df['HIGHEST_PX_LEVEL'].isna())] 
    #main_filtered_df = merged_df.drop(pre_filtered_df.index)

    # Step 4: Save the results to a new Excel file with multiple sheets
    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
        merged_df.to_excel(writer, sheet_name='All', index=False)
        pre_filtered_df.to_excel(writer, sheet_name='Pre_filtered', index=False)
        oncogenic_filtered_df.to_excel(writer, sheet_name='Oncogenic', index=False)
        actionable_filtered_df.to_excel(writer, sheet_name='Actionable', index=False)
    
        # Access the underlying XlsxWriter workbook and worksheet objects
        workbook  = writer.book
        final_sheet = writer.sheets['Actionable']

        # Set the tab color for 'Final_filtered' sheet to yellow
        final_sheet.set_tab_color('yellow')
    # Save the result to the output path
    #merged_df.to_csv(output_path, sep='\t', index=False)
    #print(f"Merged TSV saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge TSV files based on matching columns.")
    
    # Define arguments
    parser.add_argument('--cosmic', required=True, help="Path to the first TSV file (to be modified).")
    parser.add_argument('--okb', required=True, help="Path to the second TSV file (source for additional data).")
    parser.add_argument('--output', required=True, help="Output path for the merged TSV file.")
    
    args = parser.parse_args()
    
    # Call the merge function with parsed arguments
    merge_tsv(args.cosmic, args.okb, args.output)
