#!/usr/bin/env python3
"""
Standalone script to merge Exomiser JSON + TSV outputs with VEP annotations,
and collapse HPO and MOI fields into final CSV.

Usage:
    python merge_exomiser_vep.py \
        --vep      path/to/vep_annotated.tsv \
        --json     path/to/exomiser.json \
        --variants path/to/exomiser.variants.tsv \
        --output   path/to/output_prefix

Produces:
    {output_prefix}_variant_pheno_HPO.csv
    {output_prefix}_json_tsv_merged.csv
    {output_prefix}_collapsed.csv
    {output_prefix}_merged_vep.csv
    {output_prefix}_final.csv
"""

import argparse
import json
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Merge Exomiser JSON & TSV with VEP and collapse HPO/MOI"
    )
    p.add_argument('--vep', required=True,
                   help="VEP-annotated TSV file (with #POS etc. removed)")
    p.add_argument('--json', required=True,
                   help="Exomiser JSON output file")
    p.add_argument('--variants', required=True,
                   help="Exomiser variants TSV output file")
    p.add_argument('--output', required=True,
                   help="Output file prefix (no extension)")
    return p.parse_args()


def load_vep(path: str) -> pd.DataFrame:
    vep = pd.read_csv(path, sep='\t', encoding='latin-1', dtype=str)
    # strip leading '#' from any column names
    vep.columns = [c.lstrip('#') for c in vep.columns]
    # drop rows without POS, cast to int
    vep['POS'] = pd.to_numeric(vep['POS'], errors='coerce')
    vep = vep.dropna(subset=['POS'])
    vep['POS'] = vep['POS'].astype(int)
    # build Standardized_Variant_ID
    vep['Standardized_Variant_ID'] = (
        vep['CHROM'] + '-' +
        vep['POS'].astype(str) + '-' +
        vep['REF'] + '-' +
        vep['ALT']
    )
    return vep


def extract_variant_pheno(json_path: str) -> pd.DataFrame:
    with open(json_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    records = []
    for gene in data:
        gene_symbol = gene.get('geneSymbol')
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
            for match in matches:
                rec = {**base, **match}
                records.append(rec)

    df = pd.DataFrame(records)
    out_csv = f"{args.output}_variant_pheno_HPO.csv"
    df.to_csv(out_csv, index=False)
    print(f"Wrote variant‐pheno matches to {out_csv}")
    return df


def merge_json_tsv(df_variant_pheno: pd.DataFrame, variants_tsv: str) -> pd.DataFrame:
    exo = pd.read_csv(variants_tsv, sep='\t', dtype=str)
    exo.columns = exo.columns.str.lstrip('#')
    exo['variant_key'] = exo['ID'].str.split('_').str[0]
    exo[['Chr', 'Pos', 'Ref', 'Alt']] = (
        exo['variant_key'].str.split('-', expand=True).iloc[:, :4]
    )
    exo['Pos'] = exo['Pos'].astype(int)

    merged = pd.merge(
        df_variant_pheno, exo,
        on=['Chr', 'Pos', 'Ref', 'Alt'],
        how='outer',
        suffixes=('_json','_tsv')
    )

    # drop unwanted cols
    drop_cols = [
        "Gene", "VariantScore", "PhredScore", "FilterStatus", "ENTREZ_GENE_ID", "P-VALUE",
        "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_PHENO_SCORE", "EXOMISER_GENE_VARIANT_SCORE",
        "EXOMISER_VARIANT_SCORE", "CONTRIBUTING_VARIANT", "WHITELIST_VARIANT", "VCF_ID", "CONTIG",
        "CHANGE_LENGTH", "QUAL", "FILTER", "GENOTYPE", "HGVS", "CLINVAR_STAR_RATING",
        "GENE_CONSTRAINT_LOEUF", "GENE_CONSTRAINT_LOEUF_LOWER", "GENE_CONSTRAINT_LOEUF_UPPER",
        "MAX_FREQ_SOURCE", "MAX_FREQ", "ALL_FREQ", "MAX_PATH_SOURCE", "MAX_PATH", "ALL_PATH", "variant_key"
    ]
    merged.drop(columns=drop_cols, errors='ignore', inplace=True)

    # collapse HPO columns
    hpo_cols = [c for c in merged.columns if c.startswith("HPO_")]
    group_cols = [c for c in merged.columns if c not in hpo_cols]

    agg_dict = {col: "first" for col in group_cols}
    agg_dict.update({col: lambda x: list(dict.fromkeys(x)) for col in hpo_cols})

    collapsed = (
        merged
        .groupby(group_cols, dropna=False, as_index=False)
        .agg(agg_dict)
    )
    # turn lists into comma‐separated strings
    for col in hpo_cols:
        collapsed[col] = collapsed[col].apply(
            lambda L: ", ".join(str(x) for x in L) if isinstance(L, list) else str(L)
        )

    collapsed['Standardized_Variant_ID'] = \
        "chr" + collapsed['ID'].str.replace(r"_.*$", "", regex=True)

    out_csv = f"{args.output}_json_tsv_merged.csv"
    collapsed.to_csv(out_csv, index=False)
    print(f"Wrote JSON‐TSV merged to {out_csv}")
    return collapsed


def merge_with_vep(collapsed: pd.DataFrame, vep_df: pd.DataFrame) -> pd.DataFrame:
    merged = pd.merge(
        collapsed, vep_df,
        on="Standardized_Variant_ID", how="left"
    )

    # rename and clean up
    merged.rename(columns={
        "ID_y": "ID", "REF_y": "REF", "ALT_y": "ALT",
        "QUAL_y": "QUAL", "FILTER_y": "FILTER"
    }, inplace=True)
    merged["Base"] = merged["ID_x"].str.rsplit("_", n=1).str[0]

    # combine MOI per Base
    def combine_moi(series):
        seen = []
        for v in series.dropna():
            if v not in seen:
                seen.append(v)
        return ",".join(seen)
    combined = merged.groupby("Base")["MOI"] \
                     .apply(combine_moi).reset_index() \
                     .rename(columns={"MOI": "Combined_MOI"})
    merged = merged.merge(combined, on="Base", how="left")

    out_csv = f"{args.output}_merged_vep.csv"
    merged.to_csv(out_csv, index=False)
    print(f"Wrote merged with VEP to {out_csv}")
    return merged


def build_final(merged: pd.DataFrame, vep_df: pd.DataFrame) -> pd.DataFrame:
    # select final cols
    final_cols = [
        'EXOMISER_ACMG_CLASSIFICATION', 'EXOMISER_ACMG_EVIDENCE',
        'EXOMISER_ACMG_DISEASE_ID', 'EXOMISER_ACMG_DISEASE_NAME',
        'HPO_query_ID', 'HPO_match_ID', 'Combined_MOI', 'MOI',
        'Standardized_Variant_ID', 'RANK'
    ] + [c for c in vep_df.columns if c not in merged.columns[:5]]

    df = merged[final_cols]
    df.columns = df.columns.str.replace("EXOMISER_", "", regex=False)
    df = df.drop_duplicates(subset=['Standardized_Variant_ID', 'MOI',
                                    'ACMG_CLASSIFICATION', 'ACMG_EVIDENCE'])
    df = df.sort_values(by='RANK', ascending=True)

    out_csv = f"{args.output}_P2G.csv"
    df.to_csv(out_csv, index=False)
    print(f"Wrote final data to {out_csv}")
    return df


if __name__ == "__main__":
    args = parse_args()

    # Step 1: load VEP
    vep_df = load_vep(args.vep)

    # Step 2: extract JSON→phenotype matches
    df_variant_pheno = extract_variant_pheno(args.json)

    # Step 3: merge JSON with TSV
    collapsed = merge_json_tsv(df_variant_pheno, args.variants)

    # Step 4: merge back with VEP
    merged_vep = merge_with_vep(collapsed, vep_df)

    # Step 5: build the final dataframe
    final_df = build_final(merged_vep, vep_df)

    print("All done.")
