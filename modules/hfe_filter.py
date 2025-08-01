# hfe_filter.py

import pandas as pd
from urllib.parse import unquote

def filter_hfe(vep_df: pd.DataFrame, mane_transcript: str = "NM_000410.3") -> pd.DataFrame:
    """
    Return all HFE p.C282Y homozygotes on a given MANE transcript.

    vep_df : pd.DataFrame
        A VEP‐annotated table. Must contain at least:
          - "SYMBOL"
          - "MANE_SELECT"
          - "HGVSp_VEP"
          - "Zygosity"

    mane_transcript : str, optional (default="NM_000410.3")
        Which MANE transcript to filter on.
    """

    # 1) Subset to HFE only
    hfe_df = vep_df.loc[vep_df["SYMBOL"] == "HFE", :].copy()

    # 2) Restrict to the specified MANE transcript
    hfe_df["MANE_SELECT"] = hfe_df["MANE_SELECT"].fillna("")  # avoid NaN vs. str issues
    hfe_df = hfe_df.loc[hfe_df["MANE_SELECT"] == mane_transcript, :].copy()

    # 3) Normalize Zygosity, keep only homozygotes
    hfe_df["Zygosity"] = hfe_df["Zygosity"].fillna("").str.lower()
    mask_hom = hfe_df["Zygosity"].str.contains("hom", na=False)
    hfe_hom = hfe_df[mask_hom].copy()

    # 4) Decode any percent‐escapes in HGVSp, then look for “C282Y”
    hfe_hom["HGVSp_VEP"] = (
        hfe_hom["HGVSp_VEP"]
        .fillna("")
        .astype(str)
        .apply(unquote)
    )
    filtered_hfe = hfe_hom[
        hfe_hom["HGVSp_VEP"].str.contains("C282Y", case=False, na=False)
    ].copy()

    return filtered_hfe
