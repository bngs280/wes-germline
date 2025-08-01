# lmna_filter.py

import pandas as pd

def filter_lmna(plof_df: pd.DataFrame, vep_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a DataFrame of all LMNA pLOF + missense variants
    that are (P/LP) AND have a cardiac‐related ClinVar trait.
    """
    gene = "LMNA"
    sig_col   = "ClinVar_CLNSIG"
    trait_col = "ClinVar_CLNDN"
    cardiac_terms = ["cardiomyopathy", "arrhythmia", "acm", "arvc", "dcm", "heart"]

    def is_plp(clnstr: str) -> bool:
        for term in re.split(r"[,/]", cln_str):
            if term.strip().lower() in ("pathogenic", "likely_pathogenic"):
                return True
        return False

    # 1) Pull out LMNA pLOF rows
    lmna_plof = plof_df[plof_df["SYMBOL"] == gene].copy()
    lmna_plof[sig_col]   = lmna_plof[sig_col].fillna("")
    lmna_plof[trait_col] = lmna_plof[trait_col].fillna("").astype(str)
    plof_mask = (
        lmna_plof[sig_col].apply(is_plp)
        &
        lmna_plof[trait_col].str.lower().apply(lambda s: any(ct in s for ct in cardiac_terms))
    )
    filtered_plof = lmna_plof[plof_mask]

    # 2) Pull out LMNA missense rows
    lmna_mis = vep_df[
        (vep_df["SYMBOL"] == gene) &
        (vep_df["Consequence"] == "missense_variant")
    ].copy()
    lmna_mis[sig_col]   = lmna_mis[sig_col].fillna("")
    lmna_mis[trait_col] = lmna_mis[trait_col].fillna("").astype(str)
    mis_mask = (
        lmna_mis[sig_col].apply(is_plp)
        &
        lmna_mis[trait_col].str.lower().apply(lambda s: any(ct in s for ct in cardiac_terms))
    )
    filtered_mis = lmna_mis[mis_mask]

    # 3) Concatenate and return
    return pd.concat([filtered_plof, filtered_mis], ignore_index=True)
