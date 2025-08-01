#!/usr/bin/env python3
"""
pvs1_pipeline.py

1) Load a VEP‐annotated/exonplugin‐augmented file and filter to LoF
2) Merge with ClinGen HI/TS PMIDs
3) Merge with gene‐level LoF mechanism info
4) Compute PVS1 strength & produce final filtered table
"""

import argparse
import pandas as pd
import numpy as np

# ─── helper functions ─────────────────────────────────────────────────────────

lof_terms = {
    'stop_gained','frameshift_variant','splice_acceptor_variant',
    'splice_donor_variant','start_lost','transcript_ablation',
    'protein_altering_variant','stop_lost','transcript_amplification',
    'feature_elongation','feature_truncation','NMD_transcript_variant'
}

def is_lof(consequence_str: str) -> bool:
    """Return True if any term in consequence_str is a canonical LoF."""
    return any(c.strip() in lof_terms for c in consequence_str.split(','))

def collapse_pmids(row: pd.Series, cols: list) -> pd.Series:
    """Collapse numeric PMIDs in the given columns into a semicolon list."""
    pmids = (
        pd.to_numeric(row[cols], errors='coerce')
          .dropna()
          .astype(int)
          .astype(str)
    )
    return ";".join(pmids) if not pmids.empty else np.nan

def nmd_escapes(exon, total_exons, dist_to_end) -> bool:
    """True if variant is in last exon or ≤50 bp of penultimate exon."""
    try:
        exon = int(exon)
        total = int(total_exons)
        d = float(dist_to_end)
    except:
        return False
    if exon == total:
        return True
    if exon == total - 1 and d <= 50:
        return True
    return False

def gene_lof_intolerant(r: pd.Series) -> bool:
    """True if LOEUF<0.35 or pLI>0.9 or HI_score in {3,30}."""
    loeu = r.get('LOEUF')
    pli  = r.get('pLI_gene_value')
    hi   = r.get('HI_score')
    return (
        (pd.notna(loeu) and loeu < 0.35) or
        (pd.notna(pli) and pli > 0.9)        or
        (pd.notna(hi)  and hi in (3,30))
    )

def pct_protein_lost(r: pd.Series) -> float:
    """(% of protein C‐terminus lost)."""
    L = r.get('Protein_Length')
    P = r.get('AA_Position')
    if pd.isna(L) or pd.isna(P) or L == 0:
        return np.nan
    return (L - P) / L * 100.0

def pvs1_label(r: pd.Series) -> str:
    """
    Assign PVS1 strength based on consequence, NMD, protein loss,
    spliceAI, and LOF tolerance.
    """
    # not LoF at all?
    if r.get('LOF_Class','negative') == 'negative':
        return "NotApplicable"

    cons = r.get('Consequence','')
    ds   = r.get('SpliceAI_DS', 0.0)

    # check canonical splice + cryptic
    if cons in {'splice_acceptor_variant','splice_donor_variant'} or \
       (cons == 'splice_region_variant' and ds >= 0.5):
        label = "V2_Strong"
    # start‐lost
    elif cons == 'start_lost':
        label = "Moderate"
    # other truncations
    elif is_lof(cons):
        if r.get('escapes_nmd',False):
            lost_pct = pct_protein_lost(r)
            if pd.isna(lost_pct) or lost_pct >= 10:
                label = "V1_Strong"
            else:
                label = "Moderate"
        else:
            label = "V2_Strong"
    else:
        return "NotApplicable"

    # downgrade if gene is LoF‐tolerant
    if not r.get('gene_lof_intolerant', False):
        downgrade = {"V2_Strong":"V1_Strong","V1_Strong":"Moderate","Moderate":"Supporting"}
        label = downgrade.get(label, label)

    return label

# ─── main pipeline ────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="PVS1 pipeline for LoF variants")
    p.add_argument("--vep-file",    required=True, help="VEP+dbNSFP+exonplugin .txt")
    p.add_argument("--clingen-hi",  required=True, help="ClinGen HI/TS .tsv")
    p.add_argument("--gene-data",   required=True, help="Gene‐level LOF mechanism CSV")
    p.add_argument("--output-file", required=True, help="Final CSV output path")
    args = p.parse_args()

    # ── 1) load VEP‐annotated table & filter to LoF ─────────────────────────
    # skip '##' lines, keep header at first '#'
    vep_hdr = 0
    with open(args.vep_file) as fh:
        for i, line in enumerate(fh):
            if not line.startswith('##'):
                cols = line.lstrip('#').rstrip('\n').split('\t')
                vep_hdr = i
                break
    vep = pd.read_csv(
        args.vep_file, sep='\t', header=vep_hdr, low_memory=False
    )
    vep.columns = [c.lstrip('#') for c in vep.columns]
    # replace '-' → NaN except key cols
    excl = {"Location","#Uploaded_variation","Allele","UPLOADED_ALLELE"}
    torep = [c for c in vep.columns if c not in excl]
    vep[torep] = vep[torep].replace('-', np.nan)
    df_lof = vep[vep['Consequence'].apply(is_lof)].copy()

    # ── 2) merge ClinGen HI/TS PMIDs ───────────────────────────────────────────
    clg = pd.read_csv(args.clingen_hi, sep='\t')
    hap_cols  = [c for c in clg if c.startswith("Haploinsufficiency PMID")]
    trip_cols = [c for c in clg if c.startswith("Triplosensitivity PMID")]
    clg["Haploinsufficiency_PMIDs"] = clg.apply(collapse_pmids, axis=1, cols=hap_cols)
    clg["Triplosensitivity_PMIDs"]  = clg.apply(collapse_pmids, axis=1, cols=trip_cols)
    clg = clg.drop(columns=hap_cols + trip_cols)

    # Aggregate ClinGen table to one row per gene, combining fields
    def combine_series_semis(series):
        items = set()
        for cell in series.dropna():
            items.update(cell.split(';'))
        return ';'.join(sorted(items)) if items else np.nan

    clg = clg.groupby('Gene Symbol', as_index=False).agg({
        'Haploinsufficiency Score':     lambda s: pd.to_numeric(s, errors='coerce').max(),
        'Triplosensitivity Score':      lambda s: pd.to_numeric(s, errors='coerce').max(),
        'Haploinsufficiency Disease ID':lambda s: ';'.join(sorted(set(s.dropna()))),
        'Triplosensitivity Disease ID': lambda s: ';'.join(sorted(set(s.dropna()))),
        'Genomic Location':             lambda s: ';'.join(sorted(set(s.dropna()))),
        'Haploinsufficiency_PMIDs':     lambda s: combine_series_semis(s),
        'Triplosensitivity_PMIDs':      lambda s: combine_series_semis(s),
    })

    mg1 = df_lof.merge(
        clg,
        left_on='SYMBOL', right_on='Gene Symbol',
        how='left', validate='many_to_one'
    ).drop_duplicates()

    # ── 3) merge disease mechanism info ───────────────────────────────────────
    mech = pd.read_csv(args.gene_data)
    df = mg1.merge(
        mech,
        left_on=['SYMBOL','Haploinsufficiency Disease ID'],
        right_on=['Gene','MONDO'],
        how='left', validate='many_to_many'
    )

    # ── 4) parse and compute annotations ───────────────────────────────────────
    # SpliceAI: DS is max of fields 1–4
    ds = df['SpliceAI_pred'].str.split('|', expand=True).iloc[:,1:5]
    ds = ds.apply(pd.to_numeric, errors='coerce')
    df['SpliceAI_DS'] = ds.max(axis=1)

    # Exon/Total_Exons
    ex = df['EXON'].str.split('/', expand=True).astype(float)
    df['Exon']       = ex.iloc[:,0]
    df['Total_Exons']= ex.iloc[:,1]

    # Protein_position → AA_Position & Protein_Length
    pp = df['Protein_position'].fillna('')
    df['AA_Position']     = pd.to_numeric(pp.str.extract(r'^(\d+)')[0], errors='coerce')
    df['Protein_Length']  = pd.to_numeric(pp.str.extract(r'/(\d+)$')[0], errors='coerce')

    # NearestExonJB → Dist_abs, Which_end, Exon_len, then Distance_to_exon_end
    jb = df['NearestExonJB'].str.split('|', expand=True)
    df['Dist_abs']    = pd.to_numeric(jb.iloc[:,1], errors='coerce')
    df['Which_end']   = jb.iloc[:,2]
    df['Exon_len']    = pd.to_numeric(jb.iloc[:,3], errors='coerce')
    df['Distance_to_exon_end'] = np.where(
        df['Which_end'].str.lower()=='start',
        df['Exon_len'] - df['Dist_abs'],
        df['Dist_abs']
    )

    # NMD escape
    df['escapes_nmd'] = df.apply(
        lambda r: nmd_escapes(r['Exon'], r['Total_Exons'], r['Distance_to_exon_end']),
        axis=1
    )

    # LOF_Class default + numeric coercions
    df['LOF_Class'] = df.get('LOF_Class', False).fillna(False)
    for col in ['LOEUF','pLI_gene_value','HI_score']:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # gene LOF intolerance
    df['gene_lof_intolerant'] = df.apply(gene_lof_intolerant, axis=1)

    # assign PVS1 strength
    df['Strength'] = df.apply(pvs1_label, axis=1)

    # MAX_AF filter + strength + intolerance
    df['MAX_AF'] = pd.to_numeric(df.get('MAX_AF', np.nan), errors='coerce')
    lof = df[
        df['Strength'].isin(['V2_Strong','V1_Strong']) &
        df['gene_lof_intolerant'] &
        ((df['MAX_AF'].isna()) | (df['MAX_AF'] < 0.001))
    ].copy()

    # prioritize ClinSig + HI description
    mask = (
        (lof.get('Haploinsufficiency Description','') ==
           'Sufficient evidence for dosage pathogenicity')
        & lof['CLIN_SIG'].str.contains(
            r'\b(?:likely_pathogenic|pathogenic)\b', case=False, na=False
        )
    )
    top = lof[mask]
    rest= lof[~mask]
    lof_sorted = pd.concat([top, rest], ignore_index=True)

    # ── 5) write out ────────────────────────────────────────────────────────────
    lof_sorted.to_csv(args.output_file, index=False)
    print(f"Wrote {len(lof_sorted)} variants to {args.output_file}")

if __name__ == "__main__":
    main()
