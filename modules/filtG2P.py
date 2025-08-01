# import argparse
# import os
# import pandas as pd
# import numpy as np
# import gzip
# import re
# from pandas.api.types import is_string_dtype
# from typing import Dict, Tuple
# import csv
# from typing import Optional

# # def load_gene_moi(moi_file: str) -> Dict[str, str]:

# #     """Loads gene-to-MOI mapping from TSV with header 'Gene symbol\tInheritance'"""

# #     gene_moi = {}

# #     with open(moi_file) as f:

# #         reader = csv.DictReader(f, delimiter=',')

# #         for row in reader:

# #             try:

# #                 gene = str(row['Gene']).strip().upper()  # Force string conversion

# #                 moi = str(row['Inheritance']).strip()

# #                 if gene and moi:  # Only add non-empty entries

# #                     gene_moi[gene] = moi

# #             except (KeyError, AttributeError) as e:

# #                 print(f"Warning: Skipping malformed row - {str(e)}")

# #                 continue

# #     return gene_moi



# # def infer_moi(chrom: str, gene: str, zygosity: str, gene_moi_map: Dict[str, str]) -> Tuple[str, str]:

# #     """Determines inheritance pattern with priority: known MOI > chromosome/zygosity rules"""

# #     # Convert all inputs to string and clean

# #     try:

# #         gene = str(gene).strip().upper()

# #         chrom = str(chrom).strip().upper().replace('CHR', '')

# #         zygosity = str(zygosity).strip().lower()

        

# #         # Handle NaN/None cases

# #         if gene == 'NAN' or not gene:

# #             return "UNK", "Missing gene symbol"

            

# #         # Check if we have a known MOI for this gene

# #         if gene in gene_moi_map:

# #             return gene_moi_map[gene], "Known MOI from reference"

        

# #         # Skip if chromosome information is missing

# #         if not chrom or chrom == 'NAN':

# #             return "UNK", "Missing chromosome information"

            

# #         # Infer from chromosome and zygosity if no known MOI

# #         if chrom == 'X':

# #             if zygosity == 'hemizygous':

# #                 return "XL", "X-linked recessive (hemizygous male)"

# #             elif zygosity == 'heterozygous':

# #                 return "XL", "X-linked dominant (heterozygous female)"

# #         elif chrom == 'MT':

# #             return "MT", "Mitochondrial inheritance"

# #         elif chrom in [str(i) for i in range(1, 23)] + ['AUTO']:

# #             if zygosity == 'homozygous':

# #                 return "AR", "Autosomal recessive (homozygous)"

# #             elif zygosity == 'heterozygous':

# #                 return "AD", "Autosomal dominant (heterozygous)"

                

# #     except Exception as e:

# #         print(f"Warning: MOI inference failed for {gene} - {str(e)}")

    

# #     return "UNK", "Unknown inheritance pattern"



# # def add_moi_annotations(df: pd.DataFrame, moi_file: str) -> pd.DataFrame:

# #     """Adds MOI annotations to DataFrame"""

# #     gene_moi_map = load_gene_moi(moi_file)

    

# #     # Ensure required columns exist and handle case sensitivity

# #     required_cols = {'location', 'symbol', 'zygosity'}

# #     available_cols = {col.lower() for col in df.columns}

    

# #     if not required_cols.issubset(available_cols):

# #         missing = required_cols - available_cols

# #         raise ValueError(f"Missing required columns: {missing}")

    

# #     # Standardize column names

# #     df = df.rename(columns={

# #         next(col for col in df.columns if col.lower() == 'location'): 'LOCATION',

# #         next(col for col in df.columns if col.lower() == 'symbol'): 'SYMBOL',

# #         next(col for col in df.columns if col.lower() == 'zygosity'): 'Zygosity'

# #     })

    

# #     # Clean data before processing

# #     df['SYMBOL'] = df['SYMBOL'].astype(str).str.strip().str.upper()

# #     df['Zygosity'] = df['Zygosity'].astype(str).str.strip().str.lower()

    

# #     # Apply MOI inference

# #     df[['MOI', 'MOI_Explanation']] = df.apply(

# #         lambda row: infer_moi(

# #             chrom=row['LOCATION'].split(':')[0] if ':' in str(row['LOCATION']) else str(row['LOCATION']),

# #             gene=row['SYMBOL'],

# #             zygosity=row['Zygosity'],

# #             gene_moi_map=gene_moi_map

# #         ),

# #         axis=1,

# #         result_type='expand'

# #     )

# #     return df

# def load_gene_moi(moi_file: str) -> Dict[str, str]:
#     """Loads gene-to-MOI mapping from TSV with header 'Gene symbol\tInheritance'"""
#     gene_moi = {}
#     with open(moi_file) as f:
#         reader = csv.DictReader(f, delimiter=',')
#         for row in reader:
#             try:
#                 gene = str(row['Gene']).strip().upper()
#                 moi = str(row['Inheritance']).strip()
#                 if gene and moi:
#                     gene_moi[gene] = moi
#             except (KeyError, AttributeError) as e:
#                 print(f"Warning: Skipping malformed row - {str(e)}")
#                 continue
#     return gene_moi

# def infer_moi(chrom: str, gene: str, zygosity: str, gene_moi_map: Dict[str, str]) -> Tuple[str, str]:
#     """Determines inheritance pattern from chromosome/zygosity when no known MOI exists"""
#     try:
#         gene = str(gene).strip().upper()
#         chrom = str(chrom).strip().upper().replace('CHR', '')
#         zygosity = str(zygosity).strip().lower()
        
#         if gene == 'NAN' or not gene:
#             return "UNK", "Missing gene symbol"
            
#         if not chrom or chrom == 'NAN':
#             return "UNK", "Missing chromosome information"
            
#         if chrom == 'X':
#             if zygosity == 'hemizygous':
#                 return "XL", "X-linked recessive (hemizygous male)"
#             return "XL", "X-linked dominant (heterozygous female)"
            
#         elif chrom == 'MT':
#             return "MT", "Mitochondrial inheritance"
            
#         elif chrom in [str(i) for i in range(1, 23)] + ['AUTO']:
#             if zygosity == 'homozygous':
#                 return "AR", "Autosomal recessive (homozygous)"
#             return "AD", "Autosomal dominant (heterozygous)"
            
#     except Exception as e:
#         print(f"Warning: MOI inference failed for {gene} - {str(e)}")
    
#     return "UNK", "Unknown inheritance pattern"

# def add_dual_moi_annotations(df: pd.DataFrame, clin_gen_file: str, omim_file: str) -> pd.DataFrame:
#     """Adds MOI annotations from both ClinGen and OMIM sources"""
#     # Load both MOI sources
#     clin_gen_moi = load_gene_moi(clin_gen_file)
#     omim_moi = load_gene_moi(omim_file)
    
#     # Standardize columns
#     required_cols = {'location', 'symbol', 'zygosity'}
#     available_cols = {col.lower() for col in df.columns}
    
#     if not required_cols.issubset(available_cols):
#         missing = required_cols - available_cols
#         raise ValueError(f"Missing required columns: {missing}")
    
#     df = df.rename(columns={
#         next(col for col in df.columns if col.lower() == 'location'): 'LOCATION',
#         next(col for col in df.columns if col.lower() == 'symbol'): 'SYMBOL',
#         next(col for col in df.columns if col.lower() == 'zygosity'): 'Zygosity'
#     })
    
#     # Clean data
#     df['SYMBOL'] = df['SYMBOL'].astype(str).str.strip().str.upper()
#     df['Zygosity'] = df['Zygosity'].astype(str).str.strip().str.lower()
    
#     # Add ClinGen MOI annotations
#     df[['MOI', 'MOI_Explanation']] = df.apply(
#         lambda row: (
#             (clin_gen_moi[row['SYMBOL']], "Known MOI from ClinGen")
#             if row['SYMBOL'] in clin_gen_moi
#             else infer_moi(
#                 chrom=row['LOCATION'].split(':')[0] if ':' in str(row['LOCATION']) else str(row['LOCATION']),
#                 gene=row['SYMBOL'],
#                 zygosity=row['Zygosity'],
#                 gene_moi_map=clin_gen_moi
#             )
#         ),
#         axis=1,
#         result_type='expand'
#     )
    
#     # Add OMIM MOI annotations
#     df[['MOI_omim', 'MOI_omim_explanation']] = df.apply(
#         lambda row: (
#             (omim_moi[row['SYMBOL']], "Known MOI from OMIM")
#             if row['SYMBOL'] in omim_moi
#             else infer_moi(
#                 chrom=row['LOCATION'].split(':')[0] if ':' in str(row['LOCATION']) else str(row['LOCATION']),
#                 gene=row['SYMBOL'],
#                 zygosity=row['Zygosity'],
#                 gene_moi_map=omim_moi
#             )
#         ),
#         axis=1,
#         result_type='expand'
#     )
    
#     return df 

# # Add this dictionary to map population names to gnomAD columns
# GNOMAD_POPULATIONS = {
#     'AFR': ('gnomAD_exomes_AFR_AF', 'African'),
#     'AMR': ('gnomAD_exomes_AMR_AF', 'Americas'),
#     'EAS': ('gnomAD_exomes_EAS_AF', 'East Asian'),
#     'FIN': ('gnomAD_exomes_FIN_AF', 'Finnish'),
#     'NFE': ('gnomAD_exomes_NFE_AF', 'Non-Finnish European'),
#     'SAS': ('gnomAD_exomes_SAS_AF', 'South Asian'),
#     'ALL': ('gnomAD_exomes_AF', 'Global')
# }

# def filter_vep_file(output_file, primary_filtered, merged_vcf_vep, moi_file, omim_file, population: Optional[str] = 'SAS'):
#     """
#     Filters the annotated VEP file and writes results to the output directory.

#     Parameters:
#         vep_file (str): Path to the annotated VEP file.
#         base_name (str): Base name for the output file.
#         output_dir (str): Directory to save filtered results.
#         vcf_local_path (str): path to raw vcf inside patient dir.
#     """
#     output_columns = ["LOCATION", "SYMBOL (Gene Name)", "GENE (GENE ID)", "REF", "ALT", "Zygosity", "Codons", "SNPs/Rsid", 'HGVSc_VEP', 'HGVSp_VEP', "LOVD", "Amino_acids", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "IMPACT", "DP", "GQ", "PL","QUAL", "FILTER",
#             "PHENOTYPES", "ClinVar_CLNSIG", "clinvar_review", "clinvar_Orphanet_id", 'clinvar_MedGen_id', "clinvar_OMIM_id", "ACMG DISEASE ID (OMIM/ORPHA ID)", "PUBMED", "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score", "Effect", "LRT_pred", "LRT_score", "MutationTaster_pred", "MutationTaster_score", "ada_score", "rf_score", "SIFT4G_score", "SIFT4G_pred", "REVEL_score", "GERP++_NR" ,"GERP++_RS",
#             "GERP++_RS_rankscore", "gnomAD_exomes_AF", "gnomAD_exomes_AFR_AF", "gnomAD_exomes_AMR_AF", "gnomAD_exomes_EAS_AF", "gnomAD_exomes_FIN_AF", "gnomAD_exomes_NFE_AF", "gnomAD_exomes_SAS_AF", "1000Gp3_AF", "1000Gp3_AFR_AF", "1000Gp3_AMR_AF", "1000Gp3_EAS_AF", "1000Gp3_EUR_AF", "1000Gp3_SAS_AF", "ESP6500_AA_AF", "ESP6500_EA_AF", "ExAC_SAS_AF", "final_score", "weighted_score", "SpliceAI_pred"]

#     #vcf_filename = os.path.basename(vcf_local_path)  # Extract VCF filename from S3 path
#     #base_name = os.path.splitext(os.path.splitext(vcf_filename)[0])[0]
#     #output_file = os.path.join(output_dir, f"{base_name}_G2P_filtered_results.csv")
#     #primary_filtered = os.path.join(output_dir, f"{base_name}_primary_results.csv")

#     try:
        
#         print(f"Filtering VEP file: {merged_vcf_vep}")
#         print(f"Saving filtered results to: {output_file}")

#         # Path to the original VCF file (assumes same base name as VEP file but with .vcf.gz extension)
#         #if not os.path.exists(vcf_local_path):
#         #    raise FileNotFoundError(f"VCF file not found: {vcf_local_path}")
#         if not os.path.exists(merged_vcf_vep):
#             raise FileNotFoundError(f"VEP file not found: {merged_vcf_vep}")

#         vcf_withqual = pd.read_csv(merged_vcf_vep, sep='\t')
#         #print(vcf_local_path)
#         # Add MOI annotations if MOI file provided

#         if moi_file or omim_file:
        
#             try:
#                 vcf_withqual = add_dual_moi_annotations(vcf_withqual, moi_file, omim_file)
#                 print("Successfully added MOI annotations")
#                 #print(vcf_withqual[['LOCATION', 'SYMBOL', 'Zygosity', 'MOI']].head(50))
#                 # Debug output
#                 moi_cols = [c for c in vcf_withqual.columns if 'MOI' in c]
#                 if moi_cols:
#                     print(vcf_withqual[['SYMBOL'] + moi_cols].head(50))
               
#             except Exception as e:
#                 print(f"Warning: MOI annotation failed - {str(e)}")
#         # Export file to check when start proccessing file, so first its adding MOI or not
#         print("DDDDDDDDDDDDDDDDDDDDDDDDDDDD")
#         #vcf_withqual.to_csv("MOI_add.csv", sep=',', index=False)
#         vcf_withqual.columns = [col.lstrip('#') for col in vcf_withqual.columns]
#         print("EEEEEEEEEEEEEEEEEEEEEEEEEEEE")
#         #filter_consequences = vcf_withqual[vcf_withqual['FILTER']=='PASS']
#         filter_consequences = vcf_withqual.loc[vcf_withqual['FILTER'] == 'PASS'].copy()
#         print("VVVVVVVVVVVVVVVVVVV")
#         # 🔑 force text columns to string dtype **before** any .str calls
#         text_cols = ['ClinVar_CLNSIG', 'Consequence', 'SIFT4G_pred', 'LRT_pred',
#                      'MutationTaster_pred', 'MutationAssessor_pred',
#                      'BayesDel_noAF_pred', 'fathmm-XF_coding_pred',
#                      'MetaSVM_pred', 'MetaLR_pred', 'M-CAP_pred', 'PROVEAN_pred']
#         for c in text_cols:
#             if c in filter_consequences.columns:
#                 filter_consequences[c] = (filter_consequences[c]
#                                           .astype('string', copy=False)
#                                           .fillna(''))

#         numeric_columns = ['REVEL_rankscore', 'MutPred_score']
       
#         filter_consequences[numeric_columns] = (
#             filter_consequences[numeric_columns]
#             .apply(pd.to_numeric, errors='coerce')
#         ) 
#         mask = (
#             (filter_consequences['SIFT4G_pred'] == 'D') |
#             (filter_consequences['LRT_pred']   == 'D') |
#             (filter_consequences['MutationTaster_pred'].isin(['D','A'])) |
#             (filter_consequences['MutationAssessor_pred'] == 'H') |
#             (filter_consequences['BayesDel_noAF_pred']    == 'D') |
#             (filter_consequences['fathmm-XF_coding_pred'] == 'D') |
#             (filter_consequences['REVEL_rankscore'] >= 0.6) |
#             (filter_consequences['PROVEAN_pred'] == 'D') |
#             (filter_consequences['MetaSVM_pred']  == 'D') |
#             (filter_consequences['MetaLR_pred']   == 'D') |
#             (filter_consequences['M-CAP_pred']    == 'D') |
#             (filter_consequences['MutPred_score'] > 0.5) |
#             (filter_consequences['ClinPred_pred'] == 'D')
#         )
#         filtered_bytool = filter_consequences.loc[mask].copy()
#         filtered_byclinvar = filtered_bytool.loc[
#             ~filtered_bytool['ClinVar_CLNSIG']
#                 .str.contains(
#                     "benign|likely_benign|uncertain_significance|conflicting",
#                     case=False, na=False
#                 )
#         ].copy()
 
#         print("FILTEREDDDDDDDD bY clINVARRRR , searching Zygosity :")
#         print(filtered_byclinvar[['Zygosity']])
#         filtered_byclinvar.to_csv("filtered_byclinvar.csv", sep=',', index=False) # check point 1st -- worked
#         # Convert list columns to strings
#         for col in filtered_bytool.columns:
#             if is_string_dtype(filtered_bytool[col]):
#                 continue
#             if filtered_bytool[col].apply(lambda x: isinstance(x, list)).any():
#                 filtered_bytool[col] = filtered_bytool[col].apply(lambda x: ", ".join(map(str, x)) if isinstance(x, list) else x)

#         for col in filtered_byclinvar.columns:
#             if is_string_dtype(filtered_bytool[col]):
#                 continue 
#             if filtered_byclinvar[col].apply(lambda x: isinstance(x, list)).any():
#                 filtered_byclinvar[col] = filtered_byclinvar[col].apply(lambda x: ", ".join(map(str, x)) if isinstance(x, list) else x)
#         # Get the gnomAD column to use based on population
#         try:
#             gnomad_col, pop_label = GNOMAD_POPULATIONS[population]
#             print(f"Using gnomAD population frequency column: {gnomad_col} ({pop_label})")
#         except KeyError:
#             raise ValueError(f"Invalid population code: {population}. Must be one of: {list(GNOMAD_POPULATIONS.keys())}")

#         filtered_bytool    = filtered_bytool.loc[:, ~filtered_bytool.columns.duplicated()]
#         filtered_byclinvar = filtered_byclinvar.loc[:, ~filtered_byclinvar.columns.duplicated()]
#         filtered_possible = pd.merge(filtered_bytool, filtered_byclinvar, how='outer')
#         columns_to_check = ['clinvar_MedGen_id', 'clinvar_OMIM_id', 'clinvar_Orphanet_id', 
#                             'ClinVar_CLNDN', 'PHENOTYPES', 'LOVD']
#         filtered_tool_withphenotype = filtered_possible[filtered_possible[columns_to_check].notna().any(axis=1)]
#         filtered_tool_withphenotype.to_csv("filtered_tool_withphenotype.csv", sep=',', index=False) # check point 2nd -- worked
#         rare_variants_threshold = 0.01  # Allele frequency cutoff for rare variants
#         filter_rare_variants = filtered_tool_withphenotype[(filtered_tool_withphenotype[gnomad_col].isna()) | \
#                                (filtered_tool_withphenotype[gnomad_col].fillna(0).astype(float) < rare_variants_threshold)]
#         filter_rare_variants.to_csv("filter_rare_variants.csv", sep=',', index=False) # check point 3rd -- worked
#         filter_rare_variants=filter_rare_variants[(filter_rare_variants['DP']>10) & (filter_rare_variants['GQ']>10)]
#         filter_rare_variants.to_csv("filter_rare_variants1.csv", sep=',', index=False) # check point 4th -- worked


#         # Filter variants based on VAF thresholds
#         #filtered_by_vaf = filter_rare_variants[
#          #   ((filter_rare_variants['Zygosity'] == 'Homozygous') & (filter_rare_variants['VAF'] >= 0.85)) |
#           #  ((filter_rare_variants['Zygosity'] == 'Heterozygous') & (0.35 <= filter_rare_variants['VAF']) & (filter_rare_variants['VAF'] <= 0.65)) |
#           #  ((filter_rare_variants['Zygosity'] == 'Low VAF') & (filter_rare_variants['VAF'] > 0.01))
#         #]
#        # filtered_by_vaf.to_csv("filtered_by_vaf0.5.csv", sep=',', index=False) # check point 5th -- 0 rows
#         # Optionally, further filter using existing criteria
#         filter_rare_variants = filter_rare_variants[(filter_rare_variants['DP'] > 10) & (filter_rare_variants['GQ'] > 10)]
#         #filter_rare_variants.to_csv("filtered_by_vaf.csv", sep=',', index=False) # check point 5th -- 0 rows
#         #print("checking zygosity, Filtered by VAF:")
#         #print(filtered_by_vaf[['Zygosity']])
#         excluded_consequences = ['synonymous_variant', 'intron_variant'] #, 'upstream_gene_variant', 'downstream_gene_variant']
#         filter_consequences = filter_rare_variants[~filter_rare_variants['Consequence'].isin(excluded_consequences)]

#         ## STEP 2nd Block
#         # Add scoring column to the DataFrame
#         filter_consequences['Score'] = 0

#         # Define thresholds for predictors with associated weights
#         predictor_thresholds = {
#             'CADD_PHRED': (20, 2),  # (threshold, weight)
#             'Polyphen2_HDIV_rankscore': (0.85, 1),
#             'Polyphen2_HVAR_rankscore': (0.85, 1),
#             'ClinPred_rankscore': (0.5, 2),
#             'BayesDel_noAF_rankscore': (0.5, 1),
#             'BayesDel_addAF_rankscore': (0.5, 1),
#             'AlphaMissense_rankscore': (0.5, 1),
#             'MetaRNN_rankscore': (0.8, 1),
#             'GERP++_RS_rankscore': (2.0, 1),
#             'integrated_fitCons_rankscore': (0.5, 1),
#             'ada_score': (0.5, 1),
#             'rf_score': (0.5, 1),
#             'REVEL_rankscore': (0.6, 2),
#             'MutPred_score': (0.5, 1),
#         }

#         # Increment score based on numeric thresholds
#         for predictor, (threshold, weight) in predictor_thresholds.items():
#             if predictor in filter_consequences.columns:
#                 filter_consequences['Score'] += filter_consequences[predictor].apply(
#                     lambda x: weight if pd.notna(x) and float(x) > threshold else 0
#                 )

#         # Define categorical predictors with weights
#         categorical_predictors = {
#             'SIFT4G_pred': (['D'], 1),  # (damaging values, weight)
#             'LRT_pred': (['D'], 1),
#             'MutationTaster_pred': (['D', 'A'], 1),
#             'MutationAssessor_pred': (['H','M'], 1),
#             'MetaSVM_pred': (['D'], 1),
#             'MetaLR_pred': (['D'], 1),
#             'BayesDel_noAF_pred': (['D'], 1),
#             'fathmm-XF_coding_pred': (['D'], 1),
#             'M-CAP_pred': (['D'], 1),
#             'PROVEAN_pred': (['D'], 1),
#         }

#         # Increment score for categorical predictors
#         for predictor, (damaging_values, weight) in categorical_predictors.items():
#             if predictor in filter_consequences.columns:
#                 filter_consequences['Score'] += filter_consequences[predictor].apply(
#                     lambda x: weight if str(x) in damaging_values else 0
#                 )

#         # Add scoring for clinical annotations
#         if 'ClinVar_CLNSIG' in filter_rare_variants.columns:
#             filter_consequences['Score'] += filter_consequences['ClinVar_CLNSIG'].str.contains(
#                 "pathogenic|likely_pathogenic", na=False).astype(int) * 2

#         # Add scoring for rarity
#         if gnomad_col in filter_consequences.columns:
#             filter_consequences['Score'] += filter_consequences[gnomad_col].apply(
#                 lambda x: 2 if pd.isna(x) or float(x) < rare_variants_threshold else 0
#             )

#         # Add scoring for variant impact
#         if 'Consequence' in filter_consequences.columns:
#             filter_consequences['Score'] += filter_consequences['Consequence'].str.contains(
#                 "frameshift|stop_gained|splice", na=False).astype(int) * 3

#         # Add scoring for quality metrics
#         if 'DP' in filter_consequences.columns and 'GQ' in filter_consequences.columns:
#             filter_consequences['Score'] += ((filter_consequences['DP'] > 10) & (filter_consequences['GQ'] > 10)).astype(int)

#         # Sort by the final score
#         prioritized_variants = filter_consequences.sort_values(by='Score', ascending=False)
#         prioritized_variants.to_csv("prioritized_variants.csv", sep=',', index=False) # check point 5th -- 0 rows

#         # Display top variants with relevant columns
#         columns_to_display = ['CHROM', 'POS', 'REF', 'ALT', 'Zygosity', 'Score', 'PHENOTYPES'] + list(predictor_thresholds.keys()) + list(categorical_predictors.keys())
#         print(prioritized_variants[columns_to_display])
#         print(prioritized_variants)

#         ## Block 3 #
#         tool_columns = {
#             'REVEL_score_normalized': 'REVEL_rankscore',
#             'LoFtool_normalized': 'LoFtool',  # Deleterious toward 0
#             'CADD_raw_rankscore_normalized': 'CADD_PHRED',
#             'Polyphen2_HDIV_score_normalized': 'Polyphen2_HDIV_rankscore',
#             'Polyphen2_HVAR_score_normalized': 'Polyphen2_HVAR_rankscore',
#             'AlphaMissense_score_normalized': 'AlphaMissense_rankscore',
#             'VARITY_ER_score_normalized': 'VARITY_ER_rankscore',
#             'SIFT4G_score_normalized': 'SIFT4G_converted_rankscore',  # Deleterious toward 0
#             'MutationTaster_score_normalized': 'MutationTaster_converted_rankscore',  # Deleterious toward 0
#             'fathmm-XF_coding_score_normalized': 'fathmm-XF_coding_rankscore',  # Deleterious toward 0
#             'MutationTaster_score_normalized': 'MutationTaster_converted_rankscore',
#             'BayesDel_addAF_score_normalized': 'BayesDel_addAF_rankscore',
#             'BayesDel_noAF_score_normalized': 'BayesDel_noAF_rankscore',
#             'DANN_score_normalized': 'DANN_rankscore',
#             'MetaLR_score_normalized': 'MetaLR_rankscore',
#             'MetaSVM_score_normalized': 'MetaSVM_rankscore',
#             'MutPred_score_normalized': 'MutPred_rankscore',
#             'VEST4_score_normalized': 'VEST4_rankscore'
#         }
#         weights = {
#             'REVEL_score_normalized': 0.7,
#             'LoFtool_normalized': 0.4,  # Higher weight since it's intolerance to changes
#             'CADD_raw_rankscore_normalized': 0.7,
#             'Polyphen2_HDIV_score_normalized': 0.6,
#             'Polyphen2_HVAR_score_normalized': 0.6,
#             'AlphaMissense_score_normalized': 0.6,
#             'VARITY_ER_score_normalized': 0.6,
#             'SIFT4G_score_normalized': 0.5,  # Toward 0
#             'MutationTaster_score_normalized': 0.5,  # Toward 0
#             'fathmm-XF_coding_score_normalized': 0.5,  # Toward 0
#             'BayesDel_addAF_score_normalized': 0.5,
#             'BayesDel_noAF_score_normalized': 0.5,
#             'DANN_score_normalized': 0.4,
#             'MetaLR_score_normalized': 0.4,
#             'MetaSVM_score_normalized': 0.4,
#             'MutPred_score_normalized': 0.3,
#             'VEST4_score_normalized': 0.3,
#             'MutationTaster_score_normalized' : 0.5
#         }
#         invert_tools = [
#             'LoFtool_normalized',
#             'SIFT4G_score_normalized'
#         ]

#         for tool, col in tool_columns.items():
#             if col in prioritized_variants.columns:
#                 # Attempt to convert to numeric
#                 prioritized_variants[col] = pd.to_numeric(prioritized_variants[col], errors='coerce')


#         def normalize_and_invert_scores(df, tool_columns, invert_tools):
#             """
#             Normalize the scores and optionally invert for tools where lower scores are deleterious.
#             """
#             for tool, col in tool_columns.items():
#                 if col in df.columns:
#                     max_val = df[col].max()
#                     min_val = df[col].min()
#                     # Avoid division by zero with a small offset
#                     if max_val > min_val:
#                         df[f'{tool}_normalized'] = (df[col] - min_val) / (max_val - min_val + 1e-9)
#                     else:
#                         df[f'{tool}_normalized'] = 0  # Assign 0 if all values are the same
            
#                     # Invert scores for tools deleterious toward 0
#                     if tool in invert_tools:
#                         df[f'{tool}_normalized'] = 1 - df[f'{tool}_normalized']
#             return df

#         # Normalize and invert scores as needed
#         normalized = normalize_and_invert_scores(prioritized_variants, tool_columns, invert_tools)

#         def compute_weighted_scores(df, tool_columns, weights):
#             df['weighted_score'] = 0  # Initialize the weighted score column
#             for tool, weight in weights.items():
#                 if tool in df.columns:
#                     df['weighted_score'] += df[tool] * weight
#             return df


#         # Compute weighted scores using updated normalized and inverted scores
#         weighted = compute_weighted_scores(normalized, tool_columns, weights)

#         # Sort by weighted score
#         sorted_weight = weighted.sort_values(by='weighted_score', ascending=False)

#         # Display results
#         columns_to_display = ['CHROM', 'POS', 'REF', 'ALT', 'weighted_score', 'ClinVar_CLNSIG', 'Zygosity', 'Score', 'PHENOTYPES']
#         topp_variants = sorted_weight.nlargest(10, 'weighted_score')  # Top 10 variants

#         # Define scoring criteria for clinical annotations
#         clinvar_scores = {
#             "pathogenic": 3,
#             "likely_pathogenic": 2,
#             "benign": -1,
#             "likely_benign": -1,
#             "uncertain_significance": 0
#         }

#         # Function to map CLIN_SIG to a score
#         def assign_clinvar_score(clinsig):
#             if pd.isna(clinsig):  # Handle missing values
#                 return 0  
#             terms = re.split(r'[,/]', clinsig.lower())  # Split by comma and convert to lowercase
#             scores = [clinvar_scores.get(term.strip(), 0) for term in terms]  # Get scores for each term
#             return max(scores)  # Return the highest score


#         # Apply scoring to the DataFrame
#         if 'ClinVar_CLNSIG' in sorted_weight.columns:
#             sorted_weight['ClinVar_Score'] = sorted_weight['ClinVar_CLNSIG'].apply(assign_clinvar_score)

#         # Ensure weighted_score exists
#         if 'weighted_score' not in sorted_weight.columns:
#             sorted_weight['weighted_score'] = 0

#         # Combine scores into a final composite score
#         sorted_weight['final_score'] = (
#             sorted_weight['weighted_score'] + sorted_weight['ClinVar_Score']
#         )
#         print("Checking zygosity in sorted_weight: ")
#         print(sorted_weight[['Zygosity']])
        
#         # Create a temporary boolean column: True if Codons is not null, False otherwise.
#         sorted_weight['has_Codons'] = sorted_weight['Codons'].notna()

#         # Sort first by the boolean (so that True comes first) and then by final_score in descending order.
#         sorted_weight = sorted_weight.sort_values(by=['has_Codons', 'final_score'], ascending=[False, False])

#         # Optionally, drop the temporary column if it's no longer needed.
#         sorted_weight.drop(columns=['has_Codons'], inplace=True)
#         # Fix the merge_prioritized function to handle empty DataFrames
#         def merge_prioritized(row):
#             # Check for OMIM (highest priority)
#             if pd.notnull(row.get('clinvar_OMIM_id', np.nan)):
#                 return f"OMIM:{row['clinvar_OMIM_id']}"
#             elif pd.notnull(row.get('clinvar_Orphanet_id', np.nan)):
#                 return f"Orpha:{row['clinvar_Orphanet_id']}"
#             elif pd.notnull(row.get('clinvar_MedGen_id', np.nan)):
#                 return f"MedGen:{row['clinvar_MedGen_id']}"
#             return ""

#         # Only try to add ACMG column if DataFrame is not empty
#         if not sorted_weight.empty:
#             sorted_weight['ACMG DISEASE ID (OMIM/ORPHA ID)'] = sorted_weight.apply(merge_prioritized, axis=1)
#         else:
#             sorted_weight['ACMG DISEASE ID (OMIM/ORPHA ID)'] = ""       
 
#         desired_columns = [
#             "Adjusted_Location", "SYMBOL", "Gene", "REF", "ALT", "Zygosity", "Codons", "Existing_variation", 'HGVSc_VEP', 'HGVSp_VEP', "LOVD", "Amino_acids", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "IMPACT", "DP", "GQ", "PL","QUAL", "FILTER",
#             "PHENOTYPES", "ClinVar_CLNSIG", "clinvar_review", "clinvar_Orphanet_id", 'clinvar_MedGen_id', "clinvar_OMIM_id", "ACMG DISEASE ID (OMIM/ORPHA ID)", "PUBMED", "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score", "Consequence", "LRT_pred", "LRT_score", "MutationTaster_pred", "MutationTaster_score", "ada_score", "rf_score", "SIFT4G_score", "SIFT4G_pred", "REVEL_score", "GERP++_NR" ,"GERP++_RS",
#             "GERP++_RS_rankscore", "gnomAD_exomes_AF", "gnomAD_exomes_AFR_AF", "gnomAD_exomes_AMR_AF", "gnomAD_exomes_EAS_AF", "gnomAD_exomes_FIN_AF", "gnomAD_exomes_NFE_AF", "gnomAD_exomes_SAS_AF", "1000Gp3_AF", "1000Gp3_AFR_AF", "1000Gp3_AMR_AF", "1000Gp3_EAS_AF", "1000Gp3_EUR_AF", "1000Gp3_SAS_AF", "ESP6500_AA_AF", "ESP6500_EA_AF", "ExAC_SAS_AF", "final_score", "weighted_score", "SpliceAI_pred"] 
#         sorted_weight = sorted_weight[desired_columns]

#         zygosity_mapping = {
#             'Heterozygous': 'HET',
#             'Homozygous': 'HOM'
#         }
#         sorted_weight.insert(sorted_weight.columns.get_loc('Zygosity') + 1,
#                              "Zygosity_label",
#                              sorted_weight["Zygosity"].map(zygosity_mapping))
#         clinvar_mapping = {
#             "likely benign": "LB",
#             "benign": "BEN",
#             "likely pathogenic": "LP",
#             "pathogenic": "PAT",
#             "uncertain significance": "VUS",
#             "conflicting interpretations of pathogenicity": "Conflicting" 
#         }
#         def map_clinsig(clinsig):
#             if pd.isna(clinsig):  # Handle missing values
#                 return "NA"
#             pattern = r'([^/,]+)'   # Matches any sequence of characters except `/` and `,`
#             def replace_match(match):
#                 term = match.group(0).strip().lower().replace('_', ' ')  #Extract term and normalize
#                 return clinvar_mapping.get(term, match.group(0).strip())
#             mapped_clinsig = re.sub(pattern, replace_match, clinsig)
#             return mapped_clinsig

#         sorted_weight.insert(sorted_weight.columns.get_loc('ClinVar_CLNSIG') + 1, "CLINSIG_label", sorted_weight['ClinVar_CLNSIG'].apply(map_clinsig))

#         sorted_weight['merged_id'] = sorted_weight[['clinvar_Orphanet_id', 'clinvar_MedGen_id', 'clinvar_OMIM_id']]\
#             .apply(lambda x: ';'.join(x.dropna().astype(str)), axis=1)
#         sorted_weight = sorted_weight.rename(columns={
#             'Gene': 'GENE (GENE ID)',
#             'SYMBOL': 'SYMBOL (Gene Name)',
#             'Existing_variation': 'SNPs/Rsid',
#             'CLINVAR_CLNSIG': 'ClinVar_CLNSIG',
#             'Adjusted_Location': 'LOCATION',
#             'Consequence': 'Effect'
#         })
        

#         # Convert the column to numeric, coercing any non-numeric values to NaN
#         sorted_weight['final_score'] = pd.to_numeric(sorted_weight['final_score'], errors='coerce')

       
#         # Display the top-ranked variants
#         top_variants = sorted_weight.nlargest(40, 'final_score')
#         print("TOPPPPPPP variants along with Zygosity:")
#         #print(top_variants[['CHROM', 'POS', 'REF', 'ALT', 'Zygosity', 'weighted_score', 'CLIN_SIG', 'Score','PHENOTYPES','final_score']])
#         #print(top_variants[['LOCATION', 'REF', 'ALT', 'Zygosity', 'Zygosity_label', 'CLINVAR CLNSIG', 'CLINSIG_label', 'PHENOTYPES','final_score', 'weighted_score']]) 
#         print(top_variants[output_columns])
#         # Check if DataFrame is empty before saving
#         if top_variants.empty:
#             print("No variants passed the filtering criteria. No results to save.")
#             empty_df = pd.DataFrame(columns=output_columns)
#             empty_df.to_csv(output_file, sep=',', index=False)
#             print(f"Empty result file created at: {output_file}")
#         else:
#             print(f"Saving filtered results to: {output_file}")
#             top_variants.to_csv(output_file, sep=',', index=False)
#             print(f"Filtering complete. Results saved in {output_file}.")
#         if filter_consequences.empty:
#             # primary filter had no rows → write an empty file with the same columns
#             pd.DataFrame(columns=filter_consequences.columns) \
#               .to_csv(primary_filtered, sep=',', index=False)
#             print(f"No primary hits – wrote empty file at {primary_filtered}")
#         else:
#             print(f"Saving primary filtered consequences to: {primary_filtered}")
#             filter_consequences.to_csv(primary_filtered, sep=',', index=False)
#             print(f"Filtering complete. Primary filtered saved in {primary_filtered}.")

#     except Exception as e:
#         print(f"Error filtering VEP file: {e}")
#         # Create empty DataFrame with the predefined columns
#         empty_df = pd.DataFrame(columns=output_columns)
#         empty_df.to_csv(output_file, sep=',', index=False)
#         print(f"Error occurred. Empty result file created at: {output_file}")
#         return  # Exit the function after handling the error
    


# if __name__ == "__main__":
#     # Parse command-line arguments
#     parser = argparse.ArgumentParser(description="Filter VEP annotated file.")
#     #parser.add_argument("--vep_file", required=True, help="Path to the VEP file.")
#     ##parser.add_argument("--base_name", required=True, help="Base name for the output files.")
#     parser.add_argument("--output_file", required=True, help="file to save filtered results.")
#     parser.add_argument("--primary_filtered", required=True, help="file ta save preliminaey result.")
#     parser.add_argument("--merged_vcf_vep", required=True, help="vcf and vep merged filr from annotation.")
#     parser.add_argument("--moi_file", help="Optional: Path to gene-to-MOI mapping file.")
#     parser.add_argument("--omim_file", help="Optional: Path to OMIM gene-to-MOI mapping file.")
#     parser.add_argument("--population", default="SAS", choices=list(GNOMAD_POPULATIONS.keys()), help="Population group for gnomAD filtering (default: SAS)")

#     args = parser.parse_args()


#     # Call the filtering function
#     filter_vep_file(args.output_file, args.primary_filtered, args.merged_vcf_vep, args.moi_file, args.omim_file, population=args.population)

# ## --population AFR

#### update on 08072025
import argparse
import os
import pandas as pd
import numpy as np
import gzip
import re
from pandas.api.types import is_string_dtype
from typing import Dict, Tuple
import csv
from typing import Optional


def load_gene_moi(moi_file: str) -> Dict[str, str]:
    """Loads gene-to-MOI mapping from TSV with header 'Gene symbol\tInheritance'"""
    gene_moi = {}
    with open(moi_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                gene = str(row['Gene']).strip().upper()
                moi = str(row['Inheritance']).strip()
                if gene and moi:
                    gene_moi[gene] = moi
            except (KeyError, AttributeError) as e:
                print(f"Warning: Skipping malformed row - {str(e)}")
                continue
    return gene_moi

def infer_moi(chrom: str, gene: str, zygosity: str, gene_moi_map: Dict[str, str]) -> Tuple[str, str]:
    """Determines inheritance pattern from chromosome/zygosity when no known MOI exists"""
    try:
        gene = str(gene).strip().upper()
        chrom = str(chrom).strip().upper().replace('CHR', '')
        zygosity = str(zygosity).strip().lower()
        
        if gene == 'NAN' or not gene:
            return "UNK", "Missing gene symbol"
            
        if not chrom or chrom == 'NAN':
            return "UNK", "Missing chromosome information"
            
        if chrom == 'X':
            if zygosity == 'hemizygous':
                return "XL", "X-linked recessive (hemizygous male)"
            return "XL", "X-linked dominant (heterozygous female)"
            
        elif chrom == 'MT':
            return "MT", "Mitochondrial inheritance"
            
        elif chrom in [str(i) for i in range(1, 23)] + ['AUTO']:
            if zygosity == 'homozygous':
                return "AR", "Autosomal recessive (homozygous)"
            return "AD", "Autosomal dominant (heterozygous)"
            
    except Exception as e:
        print(f"Warning: MOI inference failed for {gene} - {str(e)}")
    
    return "UNK", "Unknown inheritance pattern"

def add_dual_moi_annotations(df: pd.DataFrame, omim_file: str) -> pd.DataFrame:
    """Adds MOI annotations from both ClinGen and OMIM sources"""
    # Load both MOI sources
    #clin_gen_moi = load_gene_moi(clin_gen_file)
    omim_moi = load_gene_moi(omim_file)
    
    # Standardize columns
    required_cols = {'location', 'symbol', 'zygosity'}
    available_cols = {col.lower() for col in df.columns}
    
    if not required_cols.issubset(available_cols):
        missing = required_cols - available_cols
        raise ValueError(f"Missing required columns: {missing}")
    
    df = df.rename(columns={
        next(col for col in df.columns if col.lower() == 'location'): 'LOCATION',
        next(col for col in df.columns if col.lower() == 'symbol'): 'SYMBOL',
        next(col for col in df.columns if col.lower() == 'zygosity'): 'Zygosity'
    })
    
    # Clean data
    df['SYMBOL'] = df['SYMBOL'].astype(str).str.strip().str.upper()
    df['Zygosity'] = df['Zygosity'].astype(str).str.strip().str.lower()
    
    # Add OMIM MOI annotations
    df[['MOI', 'MOI_Explanation']] = df.apply(
        lambda row: (
            (omim_moi[row['SYMBOL']], "Known MOI from OMIM")
            if row['SYMBOL'] in omim_moi
            else infer_moi(
                chrom=row['LOCATION'].split(':')[0] if ':' in str(row['LOCATION']) else str(row['LOCATION']),
                gene=row['SYMBOL'],
                zygosity=row['Zygosity'],
                gene_moi_map=omim_moi
            )
        ),
        axis=1,
        result_type='expand'
    )
    
    return df

# Add this dictionary to map population names to gnomAD columns
GNOMAD_POPULATIONS = {
    'AFR': ('gnomAD_exomes_AFR_AF', 'African'),
    'AMR': ('gnomAD4.1_joint_AMR_AF', 'Americas'),
    'EAS': ('gnomAD4.1_joint_EAS_AF', 'East Asian'),
    'FIN': ('gnomAD4.1_joint_FIN_AF', 'Finnish'),
    'NFE': ('gnomAD4.1_joint_NFE_AF', 'Non-Finnish European'),
    'SAS': ('gnomAD4.1_joint_SAS_AF', 'South Asian'),
    'ALL': ('gnomAD4.1_joint_AF', 'Global')
}

def filter_vep_file(primary_filtered, output_file, merged_vcf_vep, omim_file, population: Optional[str] = 'SAS'):
    """
    Filters the annotated VEP file and writes results to the output directory.

    Parameters:
        vep_file (str): Path to the annotated VEP file.
        base_name (str): Base name for the output file.
        output_dir (str): Directory to save filtered results.
        vcf_local_path (str): path to raw vcf inside patient dir.
        desired_columns = [
            "Adjusted_Location", "SYMBOL", "Gene", "REF", "ALT", "Zygosity", "Codons", "Existing_variation", 'HGVSc', 'HGVSp', "LOVD", "Amino_acids", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "IMPACT", "DP", "GQ", "PL","QUAL", "FILTER",
            "PHENOTYPES", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN", "ClinVar_CLNVI", "MONDO","OMIM","Orpha","MedGen","HPO", "ACMG DISEASE ID (OMIM/ORPHA ID)", "PUBMED", "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score", "Consequence", "MutationTaster_pred", "MutationTaster_score",  "SIFT4G_score", "SIFT4G_pred", "REVEL_score", "GERP++_NR" ,"GERP++_RS","ada_score","rf_score",
            "GERP++_RS_rankscore", "gnomAD4.1_joint_AF", "gnomAD4.1_joint_AFR_AF", "gnomAD4.1_joint_AMR_AF", "gnomAD4.1_joint_EAS_AF", "gnomAD4.1_joint_FIN_AF", "gnomAD4.1_joint_NFE_AF", "gnomAD4.1_joint_SAS_AF", "1000Gp3_AF", "1000Gp3_AFR_AF", "1000Gp3_AMR_AF", "1000Gp3_EAS_AF", "1000Gp3_EUR_AF", "1000Gp3_SAS_AF","RegeneronME_ALL_AF","RegeneronME_E_ASIA_AF", "RegeneronME_SAS_AF","final_score", "weighted_score", "SpliceAI_pred", "MOI", "MOI_Explanation"]
    """
    try:
        
        print(f"Filtering VEP file: {merged_vcf_vep}")
        print(f"Saving filtered results to: {output_file}")

        # Path to the original VCF file (assumes same base name as VEP file but with .vcf.gz extension)
        #if not os.path.exists(vcf_local_path):
        #    raise FileNotFoundError(f"VCF file not found: {vcf_local_path}")
        if not os.path.exists(merged_vcf_vep):
            raise FileNotFoundError(f"VEP file not found: {merged_vcf_vep}")

        vcf_withqual = pd.read_csv(merged_vcf_vep, sep='\t',encoding='latin-1',low_memory=False)
        #print(vcf_local_path)
        # Add MOI annotations if MOI file provided

        if omim_file:
            try:
                vcf_withqual = add_dual_moi_annotations(vcf_withqual, omim_file)
                #vcf_withqual = add_moi_annotations(vcf_withqual, moi_file)
                print("Successfully added MOI annotations")
                #print(vcf_withqual[['LOCATION', 'SYMBOL', 'Zygosity', 'MOI']].head(50))
                # Debug output
                moi_cols = [c for c in vcf_withqual.columns if 'MOI' in c]
                if moi_cols:
                    print(vcf_withqual[['SYMBOL'] + moi_cols].head(50))
               
            except Exception as e:
                print(f"Warning: MOI annotation failed - {str(e)}")
        # Export file to check when start proccessing file, so first its adding MOI or not
        #vcf_withqual.to_csv("MOI_add.csv", sep=',', index=False)
        vcf_withqual.columns = [col.lstrip('#') for col in vcf_withqual.columns]
        #filter_consequences = vcf_withqual[vcf_withqual['FILTER']=='PASS']
        filter_consequences = vcf_withqual.loc[vcf_withqual['FILTER'] == 'PASS'].copy()
        filter_consequences = filter_consequences[(filter_consequences['DP'] > 10) & (filter_consequences['GQ'] > 10)]
        filter_consequences = filter_consequences[filter_consequences['IMPACT'] != 'LOW'].copy()
        # 🔑 force text columns to string dtype **before** any .str calls
        text_cols = ['ClinVar_CLNSIG', 'Consequence', 'SIFT4G_pred', #'LRT_pred',
                     'MutationTaster_pred', 'MutationAssessor_pred',
                     'BayesDel_noAF_pred', 'fathmm-XF_coding_pred',
                     'MetaSVM_pred', 'MetaLR_pred', 'M-CAP_pred', 'PROVEAN_pred','ClinPred_pred','BayesDel_addAF_pred']
        for c in text_cols:
            if c in filter_consequences.columns:
                filter_consequences[c] = (filter_consequences[c]
                                          .astype('string', copy=False)
                                          .fillna(''))

        numeric_columns = ['REVEL_rankscore', 'MutPred_score']
       
        filter_consequences[numeric_columns] = (
            filter_consequences[numeric_columns]
            .apply(pd.to_numeric, errors='coerce')
        ) 
        mask = (
            (filter_consequences['SIFT4G_pred'] == 'D') |
            #(filter_consequences['LRT_pred']   == 'D') |
            (filter_consequences['MutationTaster_pred'].isin(['D','A'])) |
            (filter_consequences['MutationAssessor_pred'] == 'H') |
            (filter_consequences['BayesDel_noAF_pred']    == 'D') |
            (filter_consequences['BayesDel_addAF_pred']    == 'D') |
            (filter_consequences['fathmm-XF_coding_pred'] == 'D') |
            (filter_consequences['REVEL_rankscore'] >= 0.6) |
            (filter_consequences['PROVEAN_pred'] == 'D') |
            (filter_consequences['MetaSVM_pred']  == 'D') |
            (filter_consequences['MetaLR_pred']   == 'D') |
            (filter_consequences['M-CAP_pred']    == 'D') |
            (filter_consequences['MutPred_score'] > 0.5) |
            (filter_consequences['ClinPred_pred'] == 'D')
        )
        filtered_bytool = filter_consequences.loc[mask].copy()
        tool_str = filtered_bytool['ClinVar_CLNSIG'].fillna('').astype(str)
        has_benign = tool_str.str.contains(r'\b(benign|likely_benign|uncertain_significance)\b',
                                          case=False, na=False)
        has_patho = tool_str.str.contains(r'\b(pathogenic|likely_pathogenic)\b',
                                         case=False, na=False)
        has_no_class = tool_str.str.contains(r'no_classification_for_the_single_variant',
                                            case=False, na=False)
        has_uncertain_risk = tool_str.str.contains(r'uncertain_risk_allele',
                                                  case=False, na=False)
        has_other = tool_str.str.contains(r'\b(other|not_provided)\b',
                                         case=False,na=False)

        # only drop rows that have benign AND _no_ pathogenic
        pure_benign = ((has_benign| has_no_class | has_uncertain_risk | has_other) & ~has_patho)
        filtered_bytool = filtered_bytool.loc[~pure_benign].copy()
        filtered_bytool.to_csv("filtered_bytool.csv", sep=',', index=False)

        clinvar = filter_consequences['ClinVar_CLNSIG'].fillna('').astype(str)
        has_benign = clinvar.str.contains(r'\b(benign|likely_benign|uncertain_significance)\b',
                                          case=False, na=False)
        has_patho = clinvar.str.contains(r'\b(pathogenic|likely_pathogenic)\b',
                                         case=False, na=False)
        has_no_class = clinvar.str.contains(r'no_classification_for_the_single_variant',
                                            case=False, na=False)
        has_uncertain_risk = clinvar.str.contains(r'\b(uncertain_risk_allele|risk_factor)\b',
                                                  case=False, na=False)
        has_other = clinvar.str.contains(r'\b(other|not_provided)\b',
                                         case=False,na=False)
        is_blank = clinvar == ''

        # only drop rows that have benign AND _no_ pathogenic
        pure_benign = ((has_benign | is_blank | has_no_class | has_uncertain_risk | has_other) & ~has_patho)

        filtered_byclinvar = filter_consequences.loc[~pure_benign].copy()
 
        print("FILTEREDDDDDDDD bY Tools and clINVARRRR , searching pathogenic :")
        print(filtered_byclinvar[['SYMBOL','ClinVar_CLNSIG','LOCATION']])
        filtered_byclinvar.to_csv("filtered_byclinvar.csv", sep=',', index=False) # check point 1st -- worked
        # Convert list columns to strings
        for col in filtered_bytool.columns:
            if is_string_dtype(filtered_bytool[col]):
                continue
            if filtered_bytool[col].apply(lambda x: isinstance(x, list)).any():
                filtered_bytool[col] = filtered_bytool[col].apply(lambda x: ", ".join(map(str, x)) if isinstance(x, list) else x)

        for col in filtered_byclinvar.columns:
            if is_string_dtype(filtered_byclinvar[col]):
                continue 
            if filtered_byclinvar[col].apply(lambda x: isinstance(x, list)).any():
                filtered_byclinvar[col] = filtered_byclinvar[col].apply(lambda x: ", ".join(map(str, x)) if isinstance(x, list) else x)
        # Get the gnomAD column to use based on population
        try:
            gnomad_col, pop_label = GNOMAD_POPULATIONS[population]
            print(f"Using gnomAD population frequency column: {gnomad_col} ({pop_label})")
        except KeyError:
            raise ValueError(f"Invalid population code: {population}. Must be one of: {list(GNOMAD_POPULATIONS.keys())}")

        filtered_bytool    = filtered_bytool.loc[:, ~filtered_bytool.columns.duplicated()]
        filtered_byclinvar = filtered_byclinvar.loc[:, ~filtered_byclinvar.columns.duplicated()]
        filtered_possible = pd.concat([filtered_bytool, filtered_byclinvar],ignore_index=True,sort=False)
        filtered_possible = filtered_possible.drop_duplicates()
        filtered_possible.to_csv("filtered_possible.csv",index=False)
        columns_to_check = ['MedGen', 'OMIM', 'Orpha', 
                            'ClinVar_CLNDN', 'PHENOTYPES', 'LOVD']
        filtered_tool_withphenotype = filtered_possible[filtered_possible[columns_to_check].notna().any(axis=1)]
        print("variants after filtering by absence of disease id columns:")
        print(filtered_tool_withphenotype)
        #filtered_tool_withphenotype.to_csv("filtered_tool_withphenotype.csv", sep=',', index=False) # check point 2nd -- worked
        rare_variants_threshold = 0.01  # Allele frequency cutoff for rare variants
        filter_rare_variants = filtered_possible[(filtered_possible[gnomad_col].isna()) | \
                               (filtered_possible[gnomad_col].fillna(0).astype(float) < rare_variants_threshold)].copy()
        #filter_rare_variants.to_csv("filter_rare_variants.csv", sep=',', index=False) # check point 3rd -- worked
        #filter_rare_variants=filter_rare_variants[(filter_rare_variants['DP']>10) & (filter_rare_variants['GQ']>10)]
        #filter_rare_variants.to_csv("filter_rare_variants1.csv", sep=',', index=False) # check point 4th -- worked


        # Filter variants based on VAF thresholds
        #filtered_by_vaf = filter_rare_variants[
         #   ((filter_rare_variants['Zygosity'] == 'Homozygous') & (filter_rare_variants['VAF'] >= 0.85)) |
          #  ((filter_rare_variants['Zygosity'] == 'Heterozygous') & (0.35 <= filter_rare_variants['VAF']) & (filter_rare_variants['VAF'] <= 0.65)) |
          #  ((filter_rare_variants['Zygosity'] == 'Low VAF') & (filter_rare_variants['VAF'] > 0.01))
        #]
       # filtered_by_vaf.to_csv("filtered_by_vaf0.5.csv", sep=',', index=False) # check point 5th -- 0 rows
        # Optionally, further filter using existing criteria
        #filter_rare_variants = filter_rare_variants[(filter_rare_variants['DP'] > 10) & (filter_rare_variants['GQ'] > 10)]
        #filter_rare_variants.to_csv("filtered_by_vaf.csv", sep=',', index=False) # check point 5th -- 0 rows
        print("checking zygosity, Filtered by VAF:")
        print(filter_rare_variants)

        parts = filter_rare_variants['SpliceAI_pred'] \
            .fillna('0|0|0|0|0') \
            .str.split('|', expand=True)

        # Extract DS_AG–DS_DL, replace empty/“.” with 0, convert to float
        filter_rare_variants[['DS_AG','DS_AL','DS_DG','DS_DL']] = (
            parts.iloc[:, 1:5]
                .replace({'': 0, '.': 0})
                .astype(float)
        )

        # Compute the max delta‐score
        filter_rare_variants['SpliceAI_DS'] = filter_rare_variants[['DS_AG','DS_AL','DS_DG','DS_DL']].max(axis=1)
        no_splice = filter_rare_variants['SpliceAI_DS'] < 0.3
        syn = filter_rare_variants['Consequence'] == 'synonymous_variant'
        intronic_far = (filter_rare_variants['Consequence'] == 'intron_variant') & (filter_rare_variants['DISTANCE'] > 7)
        non_conserved = (filter_rare_variants['phyloP470way_mammalian'] < 0) & (filter_rare_variants['phastCons100way_vertebrate'] < 0.3)
        bp7_mask = (syn | intronic_far) & non_conserved & no_splice
        filter_rare_variants = filter_rare_variants[~bp7_mask]


        ## STEP 2nd Block
        # Add scoring column to the DataFrame
        filter_rare_variants['Score'] = 0

        # Define thresholds for predictors with associated weights
        predictor_thresholds = {
            'CADD_PHRED': (20, 2),  # (threshold, weight)
            'Polyphen2_HDIV_rankscore': (0.85, 1),
            'Polyphen2_HVAR_rankscore': (0.85, 1),
            'ClinPred_rankscore': (0.5, 2),
            'BayesDel_noAF_rankscore': (0.5, 1),
            'BayesDel_addAF_rankscore': (0.5, 1),
            'AlphaMissense_rankscore': (0.5, 1),
            'MetaRNN_rankscore': (0.8, 1),
            'GERP++_RS_rankscore': (2.0, 1),
            'phastCons470way_mammalian':(0.5, 1),
            #'integrated_fitCons_rankscore': (0.5, 1),
            'ada_score': (0.5, 1),
            'rf_score': (0.5, 1),
            'REVEL_rankscore': (0.6, 2),
            'MutPred_score': (0.5, 1),
        }

        # Increment score based on numeric thresholds
        for predictor, (threshold, weight) in predictor_thresholds.items():
            if predictor in filter_rare_variants.columns:
                filter_rare_variants['Score'] += filter_rare_variants[predictor].apply(
                    lambda x: weight if pd.notna(x) and float(x) > threshold else 0
                )

        # Define categorical predictors with weights
        categorical_predictors = {
            'SIFT4G_pred': (['D'], 1),  # (damaging values, weight)
            #'LRT_pred': (['D'], 1),
            'MutationTaster_pred': (['D', 'A'], 1),
            'MutationAssessor_pred': (['H','M'], 1),
            'MetaSVM_pred': (['D'], 1),
            'MetaLR_pred': (['D'], 1),
            #'BayesDel_noAF_pred': (['D'], 1),
            'fathmm-XF_coding_pred': (['D'], 1),
            'M-CAP_pred': (['D'], 1),
            'PROVEAN_pred': (['D'], 1),
        }

        # Increment score for categorical predictors
        for predictor, (damaging_values, weight) in categorical_predictors.items():
            if predictor in filter_rare_variants.columns:
                filter_rare_variants['Score'] += filter_rare_variants[predictor].apply(
                    lambda x: weight if str(x) in damaging_values else 0
                )
        
        def assign_clinvar_score(clinsig):
            if pd.isna(clinsig):
                return 0
            terms = re.split(r'[\,\/\|]\s*', clinsig.lower())
            scores = {
                'pathogenic': 3,
                'likely_pathogenic': 2,
                'benign': -1,
                'likely_benign': -1,
                'uncertain_significance': 0,
            }
            return max(scores.get(t.strip(), 0) for t in terms)
        
        def assign_clinvar_score(clinsig: str) -> int:
             if not isinstance(clinsig, str) or not clinsig.strip():
                 return 0
             
             tokens = re.split(r'\W+', clinsig.lower())
             tokens = [t for t in tokens if t]
             if "pathogenic" in tokens:
                return 3
             if "likely_pathogenic" in tokens:
                return 2
             if "benign" in tokens or "likely_benign" in tokens:
                return -1
             if "uncertain_significance" in tokens:
                return 0

             # everything else (risk_factor, no_classification_for_the_single_variant, blank, etc.)
             return 0

        if 'ClinVar_CLNSIG' in filter_rare_variants.columns:
            filter_rare_variants['ClinVar_Score'] = (
                filter_rare_variants['ClinVar_CLNSIG']
                    .apply(assign_clinvar_score)
            )
        else:
            filter_rare_variants['ClinVar_Score'] = 0

        if 'ClinVar_CLNSIG' in filter_rare_variants.columns:
            filter_rare_variants['ClinVar_Score'] = filter_rare_variants['ClinVar_CLNSIG'].str.contains(r'(?<![A-Za-z0-9_])pathogenic(?![A-Za-z0-9_])|(?<![A-Za-z0-9_])likely_pathogenic(?![A-Za-z0-9_])',
                                                                                                        case=False, na=False).astype(int) * 2
            
        else:
            filter_rare_variants['ClinVar_Score'] = 0
        filter_rare_variants.to_csv("filtered_possible.csv",index=False)

        # Add scoring for rarity
        if gnomad_col in filter_rare_variants.columns:
            filter_rare_variants['Score'] += filter_rare_variants[gnomad_col].apply(
                lambda x: 2 if pd.isna(x) or float(x) < rare_variants_threshold else 0
            )

        # Add scoring for variant impact
        if 'Consequence' in filter_rare_variants.columns:
            filter_rare_variants['Score'] += filter_rare_variants['Consequence'].str.contains(
                "frameshift_variant|feature_truncation|feature_elongation|stop_lost|stop_gained|start_lost|splice_acceptor_variant|splice_donor_variant|transcript_ablation|transcript_amplification", na=False).astype(int) * 3


        min_s, max_s = filter_rare_variants['Score'].min(), filter_rare_variants['Score'].max()
        filter_rare_variants['Score_norm'] = (filter_rare_variants['Score'] - min_s) / (max_s - min_s + 1e-9)
        primary_results = filter_rare_variants

        ## Block 3 #
        tool_columns = {
            'REVEL_score_normalized': 'REVEL_rankscore',
            'LoFtool_normalized': 'LoFtool',  # Deleterious toward 0
            'CADD_raw_rankscore_normalized': 'CADD_PHRED',
            'Polyphen2_HDIV_score_normalized': 'Polyphen2_HDIV_rankscore',
            'Polyphen2_HVAR_score_normalized': 'Polyphen2_HVAR_rankscore',
            'AlphaMissense_score_normalized': 'AlphaMissense_rankscore',
            'VARITY_ER_score_normalized': 'VARITY_ER_rankscore',
            'SIFT4G_score_normalized': 'SIFT4G_converted_rankscore',  # Deleterious toward 0
            'MutationTaster_score_normalized': 'MutationTaster_converted_rankscore',  # Deleterious toward 0
            'fathmm-XF_coding_score_normalized': 'fathmm-XF_coding_rankscore',  # Deleterious toward 0
            'BayesDel_addAF_score_normalized': 'BayesDel_addAF_rankscore',
            'BayesDel_noAF_score_normalized': 'BayesDel_noAF_rankscore',
            'DANN_score_normalized': 'DANN_rankscore',
            'MetaLR_score_normalized': 'MetaLR_rankscore',
            'MetaSVM_score_normalized': 'MetaSVM_rankscore',
            'MutPred_score_normalized': 'MutPred_rankscore',
            'VEST4_score_normalized': 'VEST4_rankscore'
        }
        weights = {
            'REVEL_score_normalized': 0.7,
            'LoFtool_normalized': 0.4,  # Higher weight since it's intolerance to changes
            'CADD_raw_rankscore_normalized': 0.7,
            'Polyphen2_HDIV_score_normalized': 0.6,
            'Polyphen2_HVAR_score_normalized': 0.6,
            'AlphaMissense_score_normalized': 0.6,
            'VARITY_ER_score_normalized': 0.6,
            'SIFT4G_score_normalized': 0.5,  # Toward 0
            'MutationTaster_score_normalized': 0.5,  # Toward 0
            'fathmm-XF_coding_score_normalized': 0.5,  # Toward 0
            'BayesDel_addAF_score_normalized': 0.5,
            'BayesDel_noAF_score_normalized': 0.5,
            'DANN_score_normalized': 0.4,
            'MetaLR_score_normalized': 0.4,
            'MetaSVM_score_normalized': 0.4,
            'MutPred_score_normalized': 0.3,
            'VEST4_score_normalized': 0.3

        }
        invert_tools = [
            'LoFtool_normalized',
            'SIFT4G_score_normalized'
        ]

        for raw_col in tool_columns.values():
            if raw_col in filter_rare_variants.columns:
                # Attempt to convert to numeric
                filter_rare_variants[raw_col] = pd.to_numeric(filter_rare_variants[raw_col], errors='coerce')


        def normalize_and_invert_scores(df, tool_columns, invert_tools):
            """
            Normalize the scores and optionally invert for tools where lower scores are deleterious.
            """
            for norm_col, raw_col in tool_columns.items():
                if raw_col in df.columns:
                    max_val = df[raw_col].max()
                    min_val = df[raw_col].min()
                    # Avoid division by zero with a small offset
                    if max_val > min_val:
                        df[norm_col] = (df[raw_col] - min_val) / (max_val - min_val + 1e-9)
                    else:
                        df[norm_col] = 0  # Assign 0 if all values are the same
            
                    # Invert scores for tools deleterious toward 0
                    if norm_col in invert_tools:
                        df[norm_col] = 1 - df[norm_col]
            return df

        # Normalize and invert scores as needed
        filter_rare_variants = normalize_and_invert_scores(filter_rare_variants, tool_columns, invert_tools)

        def compute_weighted_scores(df, weights):
            df['weighted_score'] = 0  # Initialize the weighted score column
            for tool, weight in weights.items():
                if tool in df.columns:
                    df['weighted_score'] += df[tool].fillna(0) * weight
            return df


        # Compute weighted scores using updated normalized and inverted scores
        filter_rare_variants = compute_weighted_scores(filter_rare_variants, weights)

        
        w_score = 0.5
        filter_rare_variants['final_score'] = (
            filter_rare_variants['weighted_score']
          + w_score * filter_rare_variants['Score_norm']
          + filter_rare_variants.get('ClinVar_Score', 0)
        )
        
        # Sort by weighted score
        filter_rare_variants['has_Codons'] = (
            filter_rare_variants.get('Codons', pd.Series([], dtype=bool)).notna()
        )
        def clin_rank(sig):
            s = str(sig or '').lower()
            if re.search(r'\bpathogenic\b', s):
                return 2
            if re.search(r'\blikely_pathogenic\b', s):
                return 1
            return 0

        filter_rare_variants['clinvar_rank'] = filter_rare_variants['ClinVar_CLNSIG'].apply(clin_rank)
        filter_rare_variants.to_csv("filter_rare_variants_all2.csv",index=False)

        # 3) Sort by clinvar_rank (desc), then final_score (desc)
        sorted_weight = filter_rare_variants.sort_values(['clinvar_rank','final_score'], ascending=[False,False])
        sorted_weight.to_csv("sorted and filter_rare_variants_all.csv",index=False)

        # Display results
        columns_to_display = ['CHROM', 'POS', 'REF', 'ALT', 'weighted_score', 'ClinVar_CLNSIG', 'Zygosity', 'Score', 'PHENOTYPES']
        high_variants = sorted_weight[columns_to_display].head(10)  # Top 10 variants
        topp_variants = sorted_weight.head(40)
        
        print("variants just after adding weight")
        print(high_variants)

        def merge_prioritized(row):
            # Check for OMIM (highest priority)
            if pd.notnull(row.get('OMIM', np.nan)):
                return f"OMIM:{row['OMIM']}"
            elif pd.notnull(row.get('Orpha', np.nan)):
                return f"Orpha:{row['Orpha']}"
            elif pd.notnull(row.get('MedGen', np.nan)):
                return f"MedGen:{row['MedGen']}"
            return ""

        # Only try to add ACMG column if DataFrame is not empty
        if not sorted_weight.empty:
            sorted_weight['ACMG DISEASE ID (OMIM/ORPHA ID)'] = sorted_weight.apply(merge_prioritized, axis=1)
        else:
            sorted_weight['ACMG DISEASE ID (OMIM/ORPHA ID)'] = ""       
 
        desired_columns = [
            "Adjusted_Location", "SYMBOL", "Gene", "REF", "ALT", "Zygosity", "Codons", "Existing_variation", 'HGVSc', 'HGVSp', "LOVD", "Amino_acids", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "IMPACT", "DP", "GQ", "PL","QUAL", "FILTER",
            "PHENOTYPES", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT", "ClinVar_CLNDN", "ClinVar_CLNVI", "MONDO","OMIM","Orpha","MedGen","HPO", "ACMG DISEASE ID (OMIM/ORPHA ID)", "PUBMED", "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score", "Consequence", "MutationTaster_pred", "MutationTaster_score",  "SIFT4G_score", "SIFT4G_pred", "REVEL_score", "GERP++_NR" ,"GERP++_RS","ada_score","rf_score",
            "GERP++_RS_rankscore", "gnomAD4.1_joint_AF", "gnomAD4.1_joint_AFR_AF", "gnomAD4.1_joint_AMR_AF", "gnomAD4.1_joint_EAS_AF", "gnomAD4.1_joint_FIN_AF", "gnomAD4.1_joint_NFE_AF", "gnomAD4.1_joint_SAS_AF", "1000Gp3_AF", "1000Gp3_AFR_AF", "1000Gp3_AMR_AF", "1000Gp3_EAS_AF", "1000Gp3_EUR_AF", "1000Gp3_SAS_AF","RegeneronME_ALL_AF","RegeneronME_E_ASIA_AF", "RegeneronME_SAS_AF","final_score", "weighted_score", "SpliceAI_pred", "MOI", "MOI_Explanation"
        ]

        sorted_weight = sorted_weight[desired_columns]

        zygosity_mapping = {
            'Heterozygous': 'HET',
            'Homozygous': 'HOM'
        }
        sorted_weight.insert(sorted_weight.columns.get_loc('Zygosity') + 1,
                             "Zygosity_label",
                             sorted_weight["Zygosity"].map(zygosity_mapping))
        clinvar_mapping = {
            "likely benign": "LB",
            "benign": "BEN",
            "likely pathogenic": "LP",
            "pathogenic": "PAT",
            "uncertain significance": "VUS",
            "conflicting interpretations of pathogenicity": "Conflicting" 
        }
        def map_clinsig(clinsig):
            if pd.isna(clinsig):  # Handle missing values
                return "NA"
            pattern = r'([^/,]+)'   # Matches any sequence of characters except `/` and `,`
            def replace_match(match):
                term = match.group(0).strip().lower().replace('_', ' ')  #Extract term and normalize
                return clinvar_mapping.get(term, match.group(0).strip())
            mapped_clinsig = re.sub(pattern, replace_match, clinsig)
            return mapped_clinsig

        sorted_weight.insert(sorted_weight.columns.get_loc('ClinVar_CLNSIG') + 1, "CLINSIG_label", sorted_weight['ClinVar_CLNSIG'].apply(map_clinsig))

        sorted_weight['merged_id'] = sorted_weight[['Orpha', 'MedGen', 'OMIM']]\
            .apply(lambda x: ';'.join(x.dropna().astype(str)), axis=1)
        sorted_weight = sorted_weight.rename(columns={
            'Gene': 'GENE (GENE ID)',
            'SYMBOL': 'SYMBOL (Gene Name)',
            'Existing_variation': 'SNPs/Rsid',
            'ClinVar_CLNSIG': 'CLINVAR CLNSIG',
            'Adjusted_Location': 'LOCATION',
            'Consequence': 'Effect'
        })
        
        # Convert the column to numeric, coercing any non-numeric values to NaN
        sorted_weight['final_score'] = pd.to_numeric(sorted_weight['final_score'], errors='coerce')
      
        # Display the top-ranked variants
        #top_variants = sorted_weight.head(30)
        #print("TOPPPPPPP variants along with Zygosity:")
        #print(top_variants)

        if topp_variants.empty:
            print("No variants passed the filtering criteria. No results to save.")
            empty_df = pd.DataFrame(columns=desired_columns)
            empty_df.to_csv(output_file, sep=',', index=False)
            print(f"Empty result file created at: {output_file}")
        else:
            print(f"Saving filtered results to: {output_file}")
            topp_variants.to_csv(output_file, sep=',', index=False)
            print(f"Filtering complete. Results saved in {output_file}.")
        if sorted_weight.empty:
            # primary filter had no rows → write an empty file with the same columns
            pd.DataFrame(columns=sorted_weight.columns) \
              .to_csv(primary_filtered, sep=',', index=False)
            print(f"No primary hits – wrote empty file at {primary_filtered}")
        else:
            print(f"Saving primary filtered consequences to: {primary_filtered}")
            sorted_weight.to_csv(primary_filtered, sep=',', index=False)
            print(f"Filtering complete. Primary filtered saved in {primary_filtered}.")

    except Exception as e:
        print(f"Error filtering VEP file: {e}")
        # Create empty DataFrame with the predefined columns
        empty_df = pd.DataFrame(columns=desired_columns)
        empty_df.to_csv(output_file, sep=',', index=False)
        print(f"Error occurred. Empty result file created at: {output_file}")
        return  # Exit the function after handling the error



if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Filter VEP annotated file.")
    #parser.add_argument("--vep_file", required=True, help="Path to the VEP file.")
    parser.add_argument("--primary_filtered", required=True, help="Base name for the output files.")
    parser.add_argument("--output_file", required=True, help="Directory to save filtered results.")
    #parser.add_argument("--vcf_local_path", required=True, help="raw vcf path in local.")
    parser.add_argument("--merged_vcf_vep", required=True, help="vcf and vep merged filr from annotation.")
    parser.add_argument("--omim_file", help="Optional: Path to gene-to-MOI mapping file.")
    parser.add_argument("--population", default="SAS", choices=list(GNOMAD_POPULATIONS.keys()), help="Population group for gnomAD filtering (default: SAS)")

    args = parser.parse_args()


    # Call the filtering function
    filter_vep_file(args.primary_filtered, args.output_file, args.merged_vcf_vep, args.omim_file, population=args.population)

## --population AFR

