# ryr2_filter.py

import pandas as pd

def filter_ryr2(vep_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns all RYR2 P/LP variants whose ClinVar disease name mentions a CPVT term.
    Expects vep_df to have at least: SYMBOL, ClinVar_CLNSIG, ClinVar_CLNDN, Consequence, etc.
    """
    gene = "RYR2"
    # 1) Subset to RYR2
    ryr2_df = vep_df[vep_df["SYMBOL"] == gene].copy()
    
    #print("columns inside vep for ryr2: ")
    #print(ryr2_df.columns.to_list())
    # 2) Normalize ClinSig + clinvar_trait
    ryr2_df["ClinVar_CLNSIG"] = ryr2_df["ClinVar_CLNSIG"].fillna("").str.lower()
    ryr2_df["ClinVar_CLNDN"] = ryr2_df["ClinVar_CLNDN"].fillna("").str.lower()
    
    # 3) Keep only P/LP
    def is_plp(clnstr):
        if isinstance(clnstr, list):
            terms = clnstr
        else:
            terms = re.split(r"[,/]", str(clnstr))

        for t in terms:
            if t.strip().lower() in ("pathogenic", "likely_pathogenic"):
                return True
        return False

    
    mask_clinsig = ryr2_df["ClinVar_CLNSIG"].str.split(",").apply(is_plp)
    ryr2_clin = ryr2_df.loc[mask_clinsig, :].copy()
 
    #print("columns inside vep for ryr2_clin: ")
    #print(ryr2_clin.columns.to_list())

    # 4) CPVT keywords
    cpvt_terms = [
        "catecholaminergic polymorphic ventricular tachycardia",
        "cpvt",
        "polymorphic ventricular tachycardia",
        "ventricular tachycardia"
    ]
    
    # 5) Keep only those whose clinvar_trait mentions a CPVT term
    mask_cpvt = ryr2_clin["ClinVar_CLNDN"].apply(
        lambda s: any(term in s for term in cpvt_terms)
    )
    filtered_ryr2 = ryr2_clin[mask_cpvt].copy()
    
    # 6) If no rows matched, return an empty DataFrame with the same columns as vep_df
    if filtered_ryr2.empty:
        return pd.DataFrame(columns=vep_df.columns)

    # 7) Otherwise, return the filtered rows (all columns)
    return filtered_ryr2
    print("Finished checking ryr2")
