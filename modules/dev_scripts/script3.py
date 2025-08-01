#!/usr/bin/env python3
"""
Standalone script for processing VEP-annotated variants with secondary findings analysis.
Accepts input files as command-line arguments.
"""

import pandas as pd
import numpy as np
import re
import importlib
import argparse
from lmna_filter import filter_lmna

from ryr2_filter import filter_ryr2

from ryr1_filter import filter_ryr1

# (B) Load HFE filter; but in order to reload it at runtime, we must first
#     import the module itself, then reload, then re‐import the function:
import hfe_filter
import importlib
importlib.reload(hfe_filter)
from hfe_filter import filter_hfe

# (C) Similarly, import and reload ryr2_filter in case you edited it:
import ryr2_filter
importlib.reload(ryr2_filter)
from ryr2_filter import filter_ryr2

# (C) Similarly, import and reload ryr2_filter in case you edited it:
import ryr1_filter
importlib.reload(ryr1_filter)
from ryr1_filter import filter_ryr1


def parse_arguments():
    """Parse command-line arguments for input files."""
    parser = argparse.ArgumentParser(description='Process VEP-annotated variants for secondary findings analysis')
    
    parser.add_argument('--vep', required=True, help='VEP annotated output file (e.g., SRR22574939_raw_VEP_vcfvep.txt)')
    parser.add_argument('--vcf', required=True, help='Filtered VCF file (e.g., SRR22574939_filtered_PASS.vcf)')
    parser.add_argument('--clingen', required=True, help='ClinGen gene curation file (e.g., ClinGen_gene_curation_list_GRCh38_cleaned.tsv)')
    parser.add_argument('--exomiser', required=True, help='Exomiser variants file (e.g., croppedd.SRR13386345_S19_hg19-exomiser.variants.tsv)')
    parser.add_argument('--acmg', required=True, help='ACMG SF gene list (e.g., ACMG SF3.2 gene list.xlsx)')
    parser.add_argument('--lof', required=True, help='Loss-of-function variants file (e.g., full_lof_SRR22574939.csv)')
    parser.add_argument('--dbnsfp', required=True, help='dbNSFP annotated file (e.g., cleaned_SRR22574939_dbNSFP4.7_clinpredAMPIREVhg38_exonplugin4.txt)')
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    # ---------------------- Read Input Files ----------------------
    # Read VEP Output
    vep_file = args.vep
    vcf_file = args.vcf
    clingen_HI = args.clingen
    dummy_exomiser = args.exomiser

    # Read VEP headers
    vep_headers = []
    with open(vep_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                vep_headers.append(line)
            elif line.startswith('#'):
                vep_headers.append(line)
                break

    vep_metadata_headers = vep_headers[:-1]
    
    # Read annotated VEP file
    annotated_df = pd.read_csv(vep_file, sep='\t', skiprows=len(vep_metadata_headers), low_memory=False)
    
    # Replace '-' with np.nan except in specific columns
    cols_to_exclude = ["Location", "#Uploaded_variation", "Allele", "UPLOADED_ALLELE"]
    cols_to_replace = [col for col in annotated_df.columns if col not in cols_to_exclude]
    annotated_df.loc[:, cols_to_replace] = annotated_df.loc[:, cols_to_replace].replace('-', np.nan)

    # Adjust variant locations
    def adjust_location(row):
        try:
            allele = str(row.get("Allele", "complex")).strip()
            location = str(row["Location"]).strip()
            uploaded_var = str(row["Uploaded_variation"]).strip()

            chrom, pos = location.split(":")
            left_str = pos.split("-")[0] if "-" in pos else pos
            left = int(left_str)
            up = int(uploaded_var.split("_")[1])

            if len(allele) == 1 and allele.upper() in "ACGT":
                return f"{chrom}:{left_str}" if "-" in pos else location

            return f"{chrom}:{left - 1}" if left != up - 1 else f"{chrom}:{left}"

        except Exception as e:
            return row.get("Location", "NA")

    annotated_df["Adjusted_Location"] = annotated_df.apply(adjust_location, axis=1)

    # Read other data files
    exomiser = pd.read_csv(dummy_exomiser, sep="\t")

    sf_list = pd.read_excel(args.acmg)
    probable_sf = sf_list[~sf_list['Gene'].isin(exomiser['GENE_SYMBOL'])]
    sf_from_vep = annotated_df[annotated_df['SYMBOL'].isin(probable_sf['Gene'])]
    plof_df = pd.read_csv(args.lof)

    allvep = pd.read_csv(args.dbnsfp, sep="\t")
    
    # Merge data
    new_merge = pd.merge(sf_from_vep, allvep, how='left', left_on="Uploaded_variation", right_on="#Uploaded_variation")
    new_merge = new_merge.loc[:, ~new_merge.columns.str.endswith("_y")]
    new_merge.columns = new_merge.columns.str.replace(r'_x$', '', regex=True)
    
    sf_from_vep = pd.merge(allvep, sf_from_vep, how='right', left_on="#Uploaded_variation", right_on="Uploaded_variation")
    sf_from_vep = sf_from_vep.loc[:, ~sf_from_vep.columns.str.endswith("_y")]
    sf_from_vep.columns = sf_from_vep.columns.str.replace(r'_x$', '', regex=True)


    # ---------------------- Secondary Findings Analysis ----------------------
    vep_df = new_merge
    moi_df = probable_sf

    # 0) Exclude special genes (handled separately)
    exclude_list = ["LMNA", "TTN", "RYR1", "HFE", "RYR2"]
    vep_df = vep_df[~vep_df["SYMBOL"].isin(exclude_list)].copy()

    # 1) Merge VEP with MOI
    merged = pd.merge(
        vep_df,
        moi_df[["Gene", "Inheritance "]],
        left_on="SYMBOL",
        right_on="Gene",
        how="inner"
    )

    # 2) Create normalized ClinVar field and is_plp boolean
    merged["ClinVar_CLNSIG_norm"] = (
        merged["ClinVar_CLNSIG"]
          .fillna("")
          .str.lower()
    )

    def term_is_plp(cln_str):
        for t in cln_str.split(","):
            t = t.strip()
            if t in ("pathogenic", "likely_pathogenic"):
                return True
        return False

    merged["is_plp"] = merged["ClinVar_CLNSIG_norm"].apply(term_is_plp)

    # 3a) AD, SD, XL: report all P/LP variants
    mask_AD_SD_XL = (
        merged["Inheritance "].isin(["AD", "SD", "XL"])
        & (merged["is_plp"])
    )

    # 3b) AR: only genes with ≥ 2 P/LP variants
    ar_plp = merged[(merged["Inheritance "] == "AR") & (merged["is_plp"])]
    ar_counts = (
        ar_plp
          .groupby("SYMBOL")
          .size()
          .reset_index(name="plp_count")
    )
    genes_ar_to_report = set(
        ar_counts.loc[ar_counts["plp_count"] >= 2, "SYMBOL"]
    )
    mask_AR = (
        (merged["Inheritance "] == "AR")
        & (merged["SYMBOL"].isin(genes_ar_to_report))
        & (merged["is_plp"])
    )

    # 4) Combine and extract
    mask_report = mask_AD_SD_XL | mask_AR
    report_df = merged.loc[mask_report].copy()
    report_df.drop(columns=["is_plp", "ClinVar_CLNSIG_norm", "Gene_y"], inplace=True)

    # 5) Print summary
    print("\nUnique ClinVar statuses in report_df:")
    print(report_df["ClinVar_CLNSIG"].unique())

    desired_cols = ["SYMBOL", "Inheritance ", "ClinVar_CLNSIG", "Location", "PHENOTYPES"]
    try:
        minimal_report = report_df[desired_cols].copy()
    except KeyError:
        minimal_report = pd.DataFrame(columns=desired_cols)

    print("\nMinimal summary of secondary findings variants (original column names):")
    print(minimal_report)

    # ---------------------- Gene-Specific Filters ----------------------
    # LMNA rule
    lmna_report = filter_lmna(plof_df=plof_df, vep_df=new_merge)
    print("\nLMNA report rows:", len(lmna_report))
    print(lmna_report)

    # RYR2 rule
    ryr2_report = filter_ryr2(new_merge)
    print("\nRYR2 report rows:", len(ryr2_report))
    print(ryr2_report.head())

    # RYR1 rule
    ryr1_report = filter_ryr1(new_merge)
    print("\nRYR1 report rows:", len(ryr1_report))
    print(ryr1_report.head())

    # HFE rule
    hfe_report = filter_hfe(new_merge)
    print("\nHFE p.C282Y homozygote report rows:", len(hfe_report))
    print(hfe_report.head())

    # ---------------------- Final Summary ----------------------
    summary_counts = {
        "LMNA": len(lmna_report),
        "RYR2": len(ryr2_report),
        "RYR1": len(ryr1_report),
        "HFE":  len(hfe_report),
    }

    if not report_df.empty:
        generic_counts = report_df.groupby("SYMBOL").size().to_dict()
        for gene, cnt in generic_counts.items():
            summary_counts[gene] = cnt

    print("\n=== Final Variant Summary ===")
    for gene, cnt in summary_counts.items():
        print(f"{gene}: {cnt} variants")

if __name__ == "__main__":
    main()