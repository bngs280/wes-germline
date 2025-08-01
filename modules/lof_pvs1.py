import pandas as pd
import numpy as np
import re
import argparse

# ─────────── Constants & Configuration ───────────
LOF_TERMS = {
    "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant",
    "start_lost", "stop_lost",
    "transcript_ablation", "transcript_amplification",
    "feature_elongation", "feature_truncation",
    "NMD_transcript_variant",
    "protein_altering_variant",
    "splice_region_variant"
}
SPLICEAI_THRESHOLD = 0.5
NMD_PENULTIMATE_EXON_BP = 50
MAX_AF_THRESHOLD = 0.001

# ─────────── Helper Functions ───────────

def validate_columns(df, required_cols):
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def parse_spliceai(df):
    # Split into at most 5 fields: [ gene, DS_AG, DS_AL, DS_DG, DS_DL | rest... ]
    cols = (
        df['SpliceAI_pred']
          .fillna('0|0|0|0|0')           # if entirely missing, pretend “0” for each DS
          .str.split('|', n=5, expand=True)
    )
    # Grab just the DS_AG–DS_DL slice, replace any “” or “.” with “0”
    ds = cols.iloc[:, 1:5].replace({'': '0', '.': '0'})
    # Convert to float
    ds = ds.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    ds.columns = ['DS_AG','DS_AL','DS_DG','DS_DL']

    # Attach back and compute the max
    df[['DS_AG','DS_AL','DS_DG','DS_DL']] = ds
    df['SpliceAI_DS'] = ds.max(axis=1)

def split_consequences(s):
    return set(c.strip() for c in str(s).split(','))


def escapes_nmd(exon, total_exons, dist_to_end):
    try:
        exon = float(exon)
        total_exons = float(total_exons)
        dist_to_end = float(dist_to_end)
    except Exception:
        return False
    if exon == total_exons:
        return True
    if exon == total_exons - 1 and dist_to_end <= NMD_PENULTIMATE_EXON_BP:
        return True
    return False


def gene_lof_intolerant(df):
    return (
        (df['LOEUF'].astype(float) < 0.35) |
        (df['pLI_gene_value'].astype(float) > 0.9) |
        (df['HI_score'] == 3)
    )


def pct_protein_lost(row):
    plen = row.Protein_Length
    pos = row.AA_Position
    if pd.isna(plen) or pd.isna(pos) or plen == 0:
        return None
    return (plen - pos) / plen * 100.0


def pvs1_label(r):
    # 1. Check if loss-of-function context
    if not (r.Consequence in LOF_TERMS or \
            (r.Consequence == 'splice_region_variant' and r.SpliceAI_DS >= SPLICEAI_THRESHOLD)):
        return 'NotApplicable'

    # 2. Base strength
    if r.Consequence == 'start_lost':
        strength = 'Moderate'
    elif r.Consequence in {'splice_acceptor_variant', 'splice_donor_variant'}:
        strength = 'V2_Strong'
    else:
        if not r.escapes_nmd:
            strength = 'V2_Strong'
        else:
            lost_pct = pct_protein_lost(r)
            strength = 'V1_Strong' if lost_pct is None or lost_pct >= 10 else 'Moderate'

    # 3. Downgrade if gene is LoF-tolerant
    if not r.gene_lof_intolerant:
        downgrade_map = {
            'V2_Strong': 'V1_Strong',
            'V1_Strong': 'Moderate',
            'Moderate': 'Supporting'
        }
        strength = downgrade_map.get(strength, strength)

    return strength

# ─────────── Main Pipeline ───────────

def run_pvs1_pipeline(input_vep_csv, gene_mechanism_csv): # clingen_vcf_csv):
    # 1) Load data
    df_vep = pd.read_csv(input_vep_csv,sep="\t",encoding='latin-1')
    mech = pd.read_csv(gene_mechanism_csv)
    #clingen = pd.read_csv(clingen_vcf_csv)
    df_vep.rename(columns={'Haploinsufficiency Score': 'HI_score'}, inplace=True)
    
    # 2) Validate
    required = [
        'Consequence', 'SpliceAI_pred', 'EXON', 'Protein_position',
        'NearestExonJB', 'LOEUF', 'pLI_gene_value', 'HI_score', 'MAX_AF',
        'SYMBOL', 'Haploinsufficiency Disease ID'
    ]
    validate_columns(df_vep, required)

    # 3) Merge mechanisms
    merged = df_vep.merge(
        mech,
        left_on=['SYMBOL', 'Haploinsufficiency Disease ID'],
        right_on=['Gene', 'MONDO'],
        how='left',
        suffixes=('_vcf', '_clinsum')
    )
    df = merged

    # 4) Parse and annotate
    parse_spliceai(df)

    # Split exon info
    df[['Exon', 'Total_Exons']] = df['EXON'].str.extract(r'^(\d+)[^/]*\/(\d+)').astype(float)

    # Protein positions
    pp = df['Protein_position'].fillna('')
    df['AA_Position'] = pd.to_numeric(pp.str.extract(r'^(\d+)')[0], errors='coerce')
    df['Protein_Length'] = pd.to_numeric(pp.str.extract(r'\/(\d+)$')[0], errors='coerce')

    # Nearest exon junction
    jb = df['NearestExonJB'].fillna('').str.split('|', expand=True)
    jb = jb.reindex(columns=range(4), fill_value=np.nan)
    jb.columns = ['Exon_ID', 'Dist_abs', 'Which_end', 'Exon_len']
    df[['Exon_ID','Dist_abs','Which_end','Exon_len']] = jb
    df[['Dist_abs', 'Exon_len']] = jb[['Dist_abs', 'Exon_len']].apply(pd.to_numeric, errors='coerce')
    df['Distance_to_exon_end'] = np.where(
        df['Which_end'].str.lower() == 'start',
        df['Exon_len'] - df['Dist_abs'],
        df['Dist_abs']
    )

    # 5) Compute NMD and gene intolerance
    df['escapes_nmd'] = escapes_nmd(
        df['Exon'], df['Total_Exons'], df['Distance_to_exon_end']
    )
    df['gene_lof_intolerant'] = gene_lof_intolerant(df)

    # 6) PVS1 labeling
    df['Strength'] = df.apply(pvs1_label, axis=1)

    # 7) Filter high-confidence LoF
    df['MAX_AF'] = pd.to_numeric(df['MAX_AF'], errors='coerce')
    lof = df[
        df['Strength'].isin(['V2_Strong', 'V1_Strong']) &
        df['gene_lof_intolerant'] &
        ((df['MAX_AF'].isna()) | (df['MAX_AF'] < MAX_AF_THRESHOLD))
    ].copy()

    # 8) Sort and return
    return lof.sort_values(['Strength', 'MAX_AF', 'LOEUF'], ascending=[True, True, True])

def main():
    parser = argparse.ArgumentParser(
        description="Run the PVS1 pipeline and output high-confidence LoF calls"
    )
    parser.add_argument(
        "vep_input",
        help="Path to the VEP-annotated TSV (e.g. merged_vcf_vep.txt)"
    )
    parser.add_argument(
        "mech_input",
        help="Path to the gene mechanism CSV (e.g. gene_data_with_lof_class.csv)"
    )
    parser.add_argument(
        "output",
        help="Where to write the high-confidence LoF CSV"
    )
    args = parser.parse_args()

    lof_df = run_pvs1_pipeline(args.vep_input, args.mech_input)
    lof_df.to_csv(args.output, index=False)
    print(f"Wrote {len(lof_df)} high-confidence LoF lines to {args.output}")

if __name__ == "__main__":
    main()

# Example usage:
# high_conf_lof = run_pvs1_pipeline('annotated_vep.csv', 'gene_data_with_lof_class.csv', 'clingen_with_vcf.csv')
# high_conf_lof.to_csv('high_confidence_lof_variants.csv', index=False)
# python pvs1_pipeline.py merged_vcf_vep.txt gene_data_with_lof_class.csv high_conf_lof.csv
# python pvs1_pipeline.py -h
