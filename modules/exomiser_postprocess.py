#!/usr/bin/env python3

import os
import re
import json
import pandas as pd
import argparse

def merge_outputs(variants_file, validated_vcf_path, merged_vcf_vep, variants_json_file, output_file):
    try:
        print("printing validated_vcf_path")
        print(validated_vcf_path)

        vcf_filename = os.path.basename(validated_vcf_path)
        base_name = os.path.splitext(os.path.splitext(vcf_filename)[0])[0]
        vep_df = pd.read_csv(merged_vcf_vep, sep='\t')
        vep_df.columns = [col.lstrip('#') for col in vep_df.columns]
        vep_df["Standardized_Variant_ID"] = (
            vep_df["CHROM"] + "-" +
            vep_df["POS"].astype(str) + "-" +
            vep_df["REF"] + "-" +
            vep_df["ALT"]
        )

        with open(variants_json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        records = []

        for gene in data:
            gene_symbol = gene.get('geneSymbol')
            # Pull out the HPO matches under HIPHIVE_PRIORITY → phenotypeEvidence
            pr = gene.get('priorityResults', {}).get('HIPHIVE_PRIORITY', {})
            phen_evidence = pr.get('phenotypeEvidence', [])
            matches = []
            for ev in phen_evidence:
                for m in ev.get('bestModelPhenotypeMatches', []):
                    matches.append({
                        'HPO_query_ID':   m['query']['id'],
                        'HPO_query_label': m['query']['label'],
                        'HPO_match_ID':   m['match']['id'],
                        'HPO_match_label': m['match']['label']
                    })
            # 3) For each variant that contributed to the gene score
            for v in gene.get('variantEvaluations', []):
                if not v.get('contributesToGeneScore'):
                    continue
                base = {
                    'Gene':           gene_symbol,
                    'Chr':            v.get('contigName'),
                    'Pos':            v.get('start'),
                    'Ref':            v.get('ref'),
                    'Alt':            v.get('alt'),
                    'VariantScore':   v.get('variantScore'),
                    'PhredScore':     v.get('phredScore'),
                    'VariantEffect':  v.get('variantEffect'),
                    'FilterStatus':   v.get('filterStatus')
                }
               # 4) Pair each variant with each HPO match
                for match in matches:
                    rec = {**base, **match}
                    records.append(rec)
        df_variant_pheno = pd.DataFrame(records)
        
        exomiser_df = pd.read_csv(variants_file, sep="\t")
        exomiser_df.columns = exomiser_df.columns.str.lstrip('#')
        exomiser_df['variant_key'] = exomiser_df['ID'].str.split('_').str[0]
        exomiser_df[['Chr', 'Pos', 'Ref', 'Alt']] = (
            exomiser_df['variant_key']
                  .str.split('-', expand=True)
                  .iloc[:, :4]
        )
        exomiser_df['Pos'] = exomiser_df['Pos'].astype(int)

        # 4) Merge on Chr, Pos, Ref, Alt
        json_tsv_merged = pd.merge(
            df_variant_pheno,
            exomiser_df,
            on=['Chr', 'Pos', 'Ref', 'Alt'],
            how='outer',
            suffixes=('_json','_tsv')
        )

        # List of columns to drop
        columns_to_drop = [
            "Gene", "VariantScore", "PhredScore", "FilterStatus", "ENTREZ_GENE_ID", "P-VALUE",
            "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_GENE_VARIANT_SCORE",
            "EXOMISER_VARIANT_SCORE", "CONTRIBUTING_VARIANT", "WHITELIST_VARIANT", "VCF_ID", "CONTIG",
            "CHANGE_LENGTH", "QUAL", "FILTER", "GENOTYPE", "HGVS", "CLINVAR_STAR_RATING",
            "GENE_CONSTRAINT_LOEUF", "GENE_CONSTRAINT_LOEUF_LOWER", "GENE_CONSTRAINT_LOEUF_UPPER",
            "MAX_FREQ_SOURCE", "MAX_FREQ", "ALL_FREQ", "MAX_PATH_SOURCE", "MAX_PATH", "ALL_PATH", "variant_key"
        ]

        # Drop columns
        json_tsv_merged.drop(columns=columns_to_drop, inplace=True, errors='ignore')
        hpo_cols      = [c for c in json_tsv_merged.columns if c.startswith("HPO_")]
        group_cols    = [c for c in json_tsv_merged.columns if c not in hpo_cols]
        agg_dict = {col: "first"           for col in group_cols}
        agg_dict.update({col: lambda x: list(dict.fromkeys(x))  for col in hpo_cols})
        collapsed = (
            json_tsv_merged
            .groupby(group_cols,  dropna=False, as_index=False)
            .agg(agg_dict)
        )
        for col in hpo_cols:
            collapsed[col] = collapsed[col].apply(
                lambda L: ", ".join(str(x) for x in L) if isinstance(L, list) else str(L)
            )

        collapsed["Standardized_Variant_ID"] = "chr" + collapsed["ID"].str.replace(r"_.*$", "", regex=True)
        print("json and variants file of exomiser merged columns:")
        print(list(collapsed.columns))

        merged_dfff = pd.merge(
            collapsed,
            vep_df,
            on="Standardized_Variant_ID",
            how="left"
        )
        column_rename_map = {
            "ID_y": "ID",
            "REF_y": "REF",
            "ALT_y": "ALT",
            "QUAL_y": "QUAL",
            "FILTER_y": "FILTER"
        }

        merged_dfff.rename(columns=column_rename_map, inplace=True)
        print("Merged dfff printing:")
        print(list(merged_dfff.columns))

        merged_dfff["Base"] = merged_dfff["ID_x"].str.rsplit("_", n=1).str[0]
        def combine_moi(moi_series):
            seen = []
            for moi in moi_series:
                if moi not in seen:
                    seen.append(moi)
            return ",".join(seen)
        # Group by 'Base' and combine the MOI values.
        combined = merged_dfff.groupby("Base")["MOI"].apply(combine_moi).reset_index().rename(columns={"MOI": "Combined_MOI"})
        combined
        print(f"print combined: {combined}")
        # Map the combined MOI back to every row using the 'Base' column.
        merged_dfff = merged_dfff.merge(combined, on="Base", how="left")
        print(f"print combined: {merged_dfff}")
        merged_dfff_unique = merged_dfff.drop_duplicates("Base")

        final_cols = [
            'EXOMISER_ACMG_CLASSIFICATION', 
            'EXOMISER_ACMG_EVIDENCE', 
            'EXOMISER_ACMG_DISEASE_ID', 
            'EXOMISER_ACMG_DISEASE_NAME',
            'HPO_query_ID',
            'HPO_match_ID', 
            'Combined_MOI',
            'MOI', 
            'Standardized_Variant_ID',
            'RANK'
        ] + [col for col in vep_df.columns if col not in merged_dfff.columns[:5]]

        final_df = merged_dfff[final_cols]
        final_df.columns = final_df.columns.str.replace("EXOMISER_", "", regex=False)
        final_df = final_df.drop_duplicates(subset=['Standardized_Variant_ID','MOI','ACMG_CLASSIFICATION','ACMG_EVIDENCE'])
        final_df = final_df.drop(columns=['Standardized_Variant_ID_vep', 'Standardized_Variant_ID'], errors='ignore')
       
        final_df = final_df.sort_values(by='RANK', ascending=True)
        print("Final_df columns: ", final_df.head()) 
        desired_columns = [
            "Adjusted_Location", "SYMBOL", "Gene", "REF", "ALT", "Zygosity", "Codons", "Existing_variation", 'HGVSc_VEP', 'HGVSp_VEP', "LOVD", "Amino_acids", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "ACMG_CLASSIFICATION", "IMPACT", "DP", "GQ", "PL","QUAL", "FILTER",
            "HPO_query_ID", "HPO_match_ID", "PHENOTYPES", "CLIN_SIG", "clinvar_review", "ACMG_DISEASE_ID", "Combined_MOI", "clinvar_Orphanet_id", 'clinvar_MedGen_id', "clinvar_OMIM_id", "PUBMED", "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score", "Consequence", "LRT_pred", "LRT_score", "MutationTaster_pred", "MutationTaster_score",  "ada_score", "rf_score", "SIFT4G_score", "SIFT4G_pred", "REVEL_score", "GERP++_NR" ,"GERP++_RS",
            "GERP++_RS_rankscore", "gnomAD_exomes_SAS_AF", "1000Gp3_SAS_AF", "ESP6500_EA_AF", "ExAC_SAS_AF", "RANK", "SpliceAI_pred"]
        final_df = final_df[desired_columns]

        final_df = final_df.rename(columns={
            'Gene': 'GENE (GENE ID)',
            'SYMBOL': 'SYMBOL (Gene Name)',
            'Existing_variation': 'SNPs/Rsid',
            'CLIN_SIG': 'CLINVAR CLNSIG',
            'ACMG_DISEASE_ID': 'ACMG DISEASE ID (OMIM/ORPHA ID)',
            'Combined_MOI': 'MOI',
            'Consequence': 'Effect',
            'Adjusted_Location': 'LOCATION'
        })

        zygosity_mapping = {
            'Heterozygous': 'HET',
            'Homozygous': 'HOM'
        }
        final_df.insert(final_df.columns.get_loc("Zygosity") + 1,
                "Zygosity_label",
                final_df["Zygosity"].map(zygosity_mapping))

        clinvar_mapping = {
            "LIKELY_BENIGN": "LB",
            "BENIGN": "Ben",
            "LIKELY_PATHOGENIC": "LP",
            "PATHOGENIC": "PAT",
            "UNCERTAIN_SIGNIFICANCE": "VUS",
            "CONFLICTING_INTERPRETATIONS_OF_PATHOGENICITY": "Conflicting"
        }
        
        def map_clinsig(clinsig):
            if pd.isna(clinsig):  # Handle missing values
                return "NA"
            pattern = r'([^/,]+)'   # Matches any sequence of characters except `/` and `,`
            def replace_match(match):
                term = match.group(0).strip().lower()  #Extract term and normalize
                return clinvar_mapping.get(term, term)
            mapped_clinsig = re.sub(pattern, replace_match, clinsig)
            return mapped_clinsig

        final_df.insert(final_df.columns.get_loc('CLINVAR CLNSIG') + 1, "CLINSIG_label", final_df['CLINVAR CLNSIG'].apply(map_clinsig))

        acmg_mapping = {
            'PATHOGENIC': 'PAT',
            'BENIGN': 'BEN',
            'UNCERTAIN_SIGNIFICANCE': 'VUS',
            'LIKELY_BENIGN': 'LB',
            'LIKELY_PATHOGENIC': 'LP'
        }
        final_df.insert(final_df.columns.get_loc("ACMG_CLASSIFICATION") + 1,
                "ACMG_CLASSIFICATION_label",
                final_df["ACMG_CLASSIFICATION"].map(acmg_mapping))

        final_df.to_csv(output_file, index=False)
        print(f"Merged results saved to {output_file}")

    except pd.errors.EmptyDataError:
        print(f"[ERROR] One or both files are empty")
        raise
    except FileNotFoundError as e:
        print(e)
        raise
    except Exception as e:
        print(f"Error merging outputs: {e}")
        raise 

def main():
    parser = argparse.ArgumentParser(description='Merge Exomiser annotation outputs')
    parser.add_argument('--variants_file', required=True, help='Exomiser variants TSV file')
    parser.add_argument('--validated_vcf_path', required=True, help='Validated VCF file path')
    parser.add_argument('--merged_vcf_vep', required=True, help='Merged VCF with VEP annotations resulted from the previous python script')
    parser.add_argument('--variants_json_file', required=True, help='Exomiser variants JSON file')
    parser.add_argument('--output_file', required=True, help='Output CSV file path')
    
    args = parser.parse_args()
    
    merge_outputs(
        variants_file=args.variants_file,
        validated_vcf_path=args.validated_vcf_path,
        merged_vcf_vep=args.merged_vcf_vep,
        variants_json_file=args.variants_json_file,
        output_file=args.output_file
    )

if __name__ == "__main__":
    main()
