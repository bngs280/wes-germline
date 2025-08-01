#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import re
import importlib

# your existing filter imports…
from lmna_filter    import filter_lmna
from ryr2_filter    import filter_ryr2
from ryr1_filter    import filter_ryr1
import hfe_filter; importlib.reload(hfe_filter)
from hfe_filter     import filter_hfe

def parse_arguments():
    p = argparse.ArgumentParser(
        description='Process VEP + Exomiser for secondary findings'
    )
    p.add_argument('--vep',         required=True,
                   help='VEP annotated TSV (with ## headers)')
    p.add_argument('--exomiser',    required=True,
                   help='Exomiser variants TSV')
    p.add_argument('--acmg',        required=True,
                   help='ACMG SF gene list (Excel or TSV)')
    p.add_argument('--lof',         required=True,
                   help='LoF variants CSV')
    p.add_argument('--report-file', required=True,
                   help='Where to write the final report TSV')
    return p.parse_args()

def main():
    args = parse_arguments()

    # ---------------------- load VEP ----------------------
    # grab and skip the header block
    vep_meta = []
    with open(args.vep, 'r') as f:
        for L in f:
            if L.startswith('##'):
                vep_meta.append(L)
            elif L.startswith('#'):
                vep_meta.append(L)
                break
    annotated = pd.read_csv(
        args.vep, sep='\t', skiprows=len(vep_meta), low_memory=False
    )
    annotated.columns = [c.lstrip('#') for c in annotated.columns]
    # replace '-' with NaN except key cols
    exclude = {"Location","#Uploaded_variation","Allele","UPLOADED_ALLELE"}
    to_replace = [c for c in annotated if c not in exclude]
    annotated[to_replace] = annotated[to_replace].replace('-', np.nan)

    # ---------------------- load Exomiser ----------------------
    exomiser = pd.read_csv(args.exomiser, sep='\t')
    # ---------------------- load ACMG SF list ----------------------
    sf_list = pd.read_excel(args.acmg) \
                 if args.acmg.lower().endswith(('.xls','.xlsx')) \
                 else pd.read_csv(args.acmg, sep='\t')
    # ---------------------- load LoF ----------------------
    plof_df = pd.read_csv(args.lof)

    # ---------------------- filter out special genes ----------------------
    exclude_genes = ["LMNA","TTN","RYR1","HFE","RYR2"]
    base = annotated[~annotated.SYMBOL.isin(exclude_genes)].copy()

    # ---------------------- merge to get inheritance ----------------------
    moi = sf_list.rename(columns={'Gene':'SYMBOL'})
    merged = pd.merge(
        base,
        moi[['SYMBOL','Inheritance ']],
        on='SYMBOL', how='inner'
    )

    # ---------------------- flag P/LP ----------------------
    merged['ClinVar_CLNSIG_norm'] = merged['ClinVar_CLNSIG'].fillna('').str.lower()
    def term_is_plp(cln_str):
        for t in re.split(r"[,/]", cln_str):
            t = t.strip()
            if t in ("pathogenic", "likely_pathogenic"):
                return True
        return False
    merged["is_plp"] = merged["ClinVar_CLNSIG_norm"].apply(term_is_plp)

    # AR: at least two P/LP in same gene
    ar = merged[(merged['Inheritance ']=='AR') & merged.is_plp]
    ar_counts = ar.groupby('SYMBOL').size().loc[lambda s: s>=2].index
    mask_ar = (merged['Inheritance ']=='AR') & merged.SYMBOL.isin(ar_counts) & merged.is_plp

    # AD/SD/XL: any P/LP
    mask_dom = merged['Inheritance '].isin(['AD','SD','XL']) & merged.is_plp

    report_df = merged[ mask_dom | mask_ar ].copy()

    # ---------------------- gene-specific rules ----------------------
    lmna_df = filter_lmna(plof_df, annotated)
    ryr2_df = filter_ryr2(annotated)
    ryr1_df = filter_ryr1(annotated)
    hfe_df  = filter_hfe(annotated)

    # ---------------------- final summary ----------------------
    # combine everything into one table if you like, or save report_df alone
    # here we just save the core secondary‐findings table:
    final = pd.concat([
        report_df,
        lmna_df,
        ryr1_df,
        ryr2_df,
        hfe_df
    ], ignore_index=True)

    final.to_csv(args.report_file, sep=',', index=False)
    print(f"Saved secondary findings report to {args.report_file}")

if __name__=='__main__':
    main()
