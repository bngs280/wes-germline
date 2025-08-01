# ryr1_filter.py

import pandas as pd

def filter_ryr1(vep_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns all RYR1 P/LP variants whose ClinVar disease name mentions malignant hyperthermia (MH).
    Expects vep_df to have at least these columns:
      - SYMBOL
      - ClinVar_CLNSIG         (ClinVar significance)
      - ClinVar_CLNDN    (ClinVar disease name)
      - Consequence      (e.g. "missense_variant")
      - plus any other columns you want to carry forward
    """
    gene = "RYR1"
    sig_col   = "ClinVar_CLNSIG"
    trait_col = "ClinVar_CLNDN"
    
    # Accept either a comma‐separated string or a pre‐split list
    def is_plp(clnstr):
        if isinstance(clnstr, list):
            terms = clnstr
        else:
            terms = re.split(r"[,/]", str(clnstr))
        for t in terms:
            if t.strip().lower() in ("pathogenic", "likely_pathogenic"):
                return True
        return False

    # 1) Subset to RYR1 only
    ryr1_df = vep_df.loc[vep_df["SYMBOL"] == gene, :].copy()
    print("▶ after SYMBOL filter, cols in ryr1 are:", ryr1_df.columns.tolist())

    # 2) Normalize ClinSig and clinvar_trait to lowercase
    #    (If CLIN_SIG is already a list, leave it; but converting to str won't hurt)
    ryr1_df[sig_col]     = ryr1_df[sig_col].fillna("").str.lower()
    ryr1_df[trait_col]   = ryr1_df[trait_col].fillna("").str.lower()

    # 3) Keep only those with ClinSig containing “pathogenic” or “likely_pathogenic”
    mask_clinsig = ryr1_df[sig_col].apply(is_plp)
    ryr1_clin = ryr1_df.loc[mask_clinsig, :].copy()
    print("▶ ryr1_clin.shape:", ryr1_clin.shape)
    print("▶ ryr1_clin columns:", ryr1_clin.columns.tolist())
    print("▶ columns containing 'CLNDN':", [c for c in ryr1_clin.columns if 'CLNDN' in c])

    # 4) Define malignant hyperthermia (MH) keywords
    mh_terms = [
        "malignant hyperthermia",
        "mh"
    ]

    # 5) Keep only those whose clinvar_trait mentions an MH term
    mask_mh = ryr1_clin[trait_col].apply(
        lambda s: any(term in s for term in mh_terms)
    )
    filtered_ryr1 = ryr1_clin.loc[mask_mh, :].copy()

    # 6) If no rows matched, return an empty DataFrame with the same columns as vep_df
    if filtered_ryr1.empty:
        return pd.DataFrame(columns=vep_df.columns)

    # 7) Otherwise, return the filtered rows (all columns)
    return filtered_ryr1

