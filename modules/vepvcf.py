# #!/usr/bin/env python3
# import argparse
# import gzip
# import sys
# import re

# import numpy as np
# import pandas as pd

# def parse_args():
#     p = argparse.ArgumentParser(
#         description="Merge VEP-annotated table, raw VCF, and ClinGen HI/TS data"
#     )
#     p.add_argument("--vep-file",      required=True,
#                    help="Annotated VEP TSV (with ## headers)")
#     p.add_argument("--vcf-file",      required=True,
#                    help="Raw VCF (.vcf or .vcf.gz)")
#     p.add_argument("--clingen-file",  required=True,
#                    help="ClinGen gene curation TSV")
#     p.add_argument("--merged-vcf-vep",required=True,
#                    help="Output path for merged VCF+VEP TSV")
#     p.add_argument("--output-file",   required=True,
#                    help="Final output CSV path")
#     p.add_argument("--sample-sex",    choices=["male","female"],
#                    default="male", help="Sex for hemizygous calls")
#     return p.parse_args()

# def read_vep(vep_file):
#     # capture the VEP headers
#     headers = []
#     with open(vep_file, "r") as f:
#         for line in f:
#             if line.startswith("##"):
#                 headers.append(line)
#             elif line.startswith("#"):
#                 headers.append(line)
#                 break
#     meta = headers[:-1]
#     df = pd.read_csv(
#         vep_file, sep="\t",
#         skiprows=len(meta),
#         low_memory=False,
#         encoding="latin1",  # handle special characters
#     )
#     # strip leading '#' from column names
#     df.columns = [c.lstrip("#") for c in df.columns]
#     # replace '-' with NaN except in key columns
#     exclude = {"Location", "Uploaded_variation", "Allele", "UPLOADED_ALLELE"}
#     cols = [c for c in df.columns if c not in exclude]
#     df.loc[:, cols] = df.loc[:, cols].replace("-", np.nan)
#     # adjust Location to CHROM:POS or CHROM:POS-1
#     def adjust_location(r):
#         try:
#             allele = str(r.get("Allele","")).strip()
#             loc    = str(r["Location"]).strip()
#             uv     = str(r["Uploaded_variation"]).strip()
#             chrom, pos = loc.split(":")
#             left = int(pos.split("-")[0])
#             up   = int(uv.split("_")[1])
#             if len(allele)==1 and allele.upper() in "ACGT":
#                 return loc if "-" not in pos else f"{chrom}:{pos.split('-')[0]}"
#             return f"{chrom}:{left-1}" if left != up-1 else f"{chrom}:{left}"
#         except:
#             return r.get("Location","")
#     df["Adjusted_Location"] = df.apply(adjust_location, axis=1)
#     return df, meta

# def read_clingen(path):
#     df = pd.read_csv(path, sep="\t", low_memory=False, dtype=str)
#     # Identify PMID columns
#     haplo_pm = [c for c in df.columns if c.startswith('Haploinsufficiency PMID')]
#     triplo_pm = [c for c in df.columns if c.startswith('Triplosensitivity PMID')]

#     # Collapse per-row PMIDs into a semicolon-separated string
#     def collapse_pmids(row, cols):
#         pmids = []
#         for col in cols:
#             val = row.get(col)
#             if pd.notna(val):
#                 try:
#                     pmids.append(str(int(float(val))))
#                 except (ValueError, TypeError):
#                     pass
#         return ';'.join(sorted(set(pmids))) if pmids else np.nan

#     df['Haploinsufficiency_PMIDs'] = df.apply(lambda r: collapse_pmids(r, haplo_pm), axis=1)
#     df['Triplosensitivity_PMIDs']  = df.apply(lambda r: collapse_pmids(r, triplo_pm), axis=1)

#     # Helper to combine semicolon-separated lists across genes
#     def combine_series_semis(series):
#         items = set()
#         for cell in series.dropna():
#             items.update(cell.split(';'))
#         return ';'.join(sorted(items)) if items else np.nan

#     # Aggregate ClinGen table to one row per gene, combining fields
#     df_agg = df.groupby('Gene Symbol', as_index=False).agg({
#         'Haploinsufficiency Score':     lambda s: pd.to_numeric(s, errors='coerce').max(),
#         'Triplosensitivity Score':      lambda s: pd.to_numeric(s, errors='coerce').max(),
#         'Haploinsufficiency Disease ID':lambda s: ';'.join(sorted(set(s.dropna()))),
#         'Triplosensitivity Disease ID': lambda s: ';'.join(sorted(set(s.dropna()))),
#         'Genomic Location':             lambda s: ';'.join(sorted(set(s.dropna()))),
#         'Haploinsufficiency_PMIDs':     lambda s: combine_series_semis(s),
#         'Triplosensitivity_PMIDs':      lambda s: combine_series_semis(s),
#     })
#     return df_agg

# def read_vcf(path):
#     import gzip
#     vcf_headers = []
#     vcf_lines   = []
#     # Use gzip.open for .gz files, open otherwise
#     opener = gzip.open if path.endswith(".gz") else open
#     with opener(path, "rt", encoding="utf-8", errors="replace") as f:
#         for line in f:
#             if line.startswith("##"):
#                 vcf_headers.append(line)
#             elif line.startswith("#CHROM"):
#                 vcf_headers.append(line)
#                 # grab the real column names (including your sample name)
#                 vcf_columns = line.rstrip().split("\t")
#             else:
#                 vcf_lines.append(line.rstrip("\n").split("\t"))

#     # build the DataFrame using the real VCF columns
#     #df = pd.DataFrame(rows, columns=header_cols)
#     fixed_columns      = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
#     df = pd.DataFrame(vcf_lines, columns=fixed_columns + ["SAMPLE_VALUES"])
#     print(df)

#     #df = pd.DataFrame(rows, columns=fixed_columns+[sample_column_name])
#     required_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']
#     if not all(col in df.columns for col in required_cols):
#         print("Error: VCF missing required columns.", file=sys.stderr)
#         sys.exit(1)

#     df["CHROM_POS"] = df["#CHROM"] + ":" + df["POS"].astype(str)
#     format_index = df.columns.get_loc("FORMAT")
#     sample_column_name = df.columns[format_index + 1]  # Dynamically detect sample column
#     print(sample_column_name)
#     #format_fields = df['FORMAT'].iloc[0].split(':')  # Extract field names from FORMAT
#     #df["CHROM_POS"] = df["#CHROM"] + ":" + df["POS"]
#     # split FORMAT and sample columns
#     fmt_splits  = df["FORMAT"].str.split(":", expand=False)
#     geno_splits = df[sample_column_name].str.split(":", expand=False)
#     print("format_fields")
#     print(fmt_splits)
#     lengths = fmt_splits.str.len().to_numpy()
#     rows_i  = np.repeat(df.index.values, lengths)
#     keys    = [k for sub in fmt_splits  for k  in sub]
#     vals    = [v for sub in geno_splits for v in sub]
#     long    = pd.DataFrame({"row": rows_i, "key": keys, "val": vals})
#     sample_df = (long.pivot(index="row", columns="key", values="val")
#                    .reindex(index=df.index)
#                    .fillna("."))
#     format_fields = sample_df.columns.tolist()
#     print("format_fields")
#     print(format_fields)

#     sample_df.columns = format_fields 
#     vcf = pd.concat([df.drop(columns=['FORMAT', sample_column_name]), sample_df], axis=1)
#     #vcf = pd.concat([df.drop(columns=["FORMAT","SAMPLE_VALUES"]), sample_df], axis=1)
#     # depth
#     if   'DP' in vcf.columns:
#         pass
#     elif 'DPI' in vcf.columns:
#         vcf['DP'] = vcf['DPI']
#     elif 'AD' in vcf.columns:
#         ad = vcf['AD'].str.split(",", expand=True).iloc[:,0:2].astype(float)
#         vcf['DP'] = ad.sum(axis=1)
#     if 'DPI' in vcf.columns:
#         vcf.drop(columns=['DPI'], inplace=True)
#     # GQ
#     if 'GQ' not in vcf.columns:
#         def calc_gq(r):
#             pl = None
#             if 'PL' in r and pd.notna(r['PL']):
#                 raw = list(map(int, str(r['PL']).split(",")))
#                 base= min(raw)
#                 pl  = [p-base for p in raw]
#             elif 'GL' in r and pd.notna(r['GL']):
#                 raw = [round(-10*float(x)) for x in str(r['GL']).split(",")]
#                 base= min(raw)
#                 pl  = [p-base for p in raw]
#             if pl is None: return np.nan
#             gt = str(r['GT']).replace("|","/")
#             idx= {'0/0':0,'0/1':1,'1/1':2}.get(gt)
#             if idx is None: return np.nan
#             return min(p for i,p in enumerate(pl) if i!=idx)
#         vcf['GQ'] = vcf.apply(calc_gq, axis=1)
#     # PL from GL
#     if 'PL' not in vcf.columns and 'GL' in vcf.columns:
#         def gl_to_pl(gl):
#             raw = [round(-10*float(x)) for x in str(gl).split(",")]
#             base= min(raw)
#             return ",".join(str(p-base) for p in raw)
#         vcf['PL'] = vcf['GL'].apply(gl_to_pl)
#     # VAF
#     if 'VAF' not in vcf.columns and 'AD' in vcf.columns and 'DP' in vcf.columns:
#         def calc_vaf(ad, dp):
#             try:
#                 a0,a1 = map(int, ad.split(","))
#                 return a1/(a0+a1) if (a0+a1)>0 else np.nan
#             except:
#                 return np.nan
#         vcf['VAF'] = vcf.apply(lambda r: calc_vaf(r['AD'], r['DP']), axis=1)
#     # coerce numeric
#     for c in ('DP','GQ','VAF'):
#         if c in vcf.columns:
#             vcf[c] = pd.to_numeric(vcf[c], errors="coerce")
#     return vcf

# def classify_zygosity(gt, vaf, chrom, sex):
#     alleles = re.split(r"[\/|]", str(gt)) if pd.notna(gt) else []
#     ref = alleles.count("0")
#     alt = len([a for a in alleles if a not in {"0",".",""}])
#     ploidy = len(alleles)
#     hap = (sex=="male" and chrom.replace("chr","") in {"X","Y"})
#     if hap:
#         if alt>=1: return "Hemizygous"
#         if ref>=1: return "Hemizygous"
#     else:
#         if alt==0 and ref>=1: return "Homozygous"
#         if alt==ploidy and ploidy>0: return "Homozygous"
#         if alt>=1 and ref>=1: return "Heterozygous"
#     if pd.notna(vaf):
#         if vaf>=0.85: return "Homozygous"
#         if vaf<0.15:  return "Homozygous"
#         if 0.35<=vaf<=0.65: return "Heterozygous"
#     return "Unclassified"

# def main():
#     args = parse_args()

#     # 1) Read VEP
#     annotated_df, vep_meta = read_vep(args.vep_file)

#     # 2) Read ClinGen
#     read_clin = read_clingen(args.clingen_file)

#     # 3) Merge VEP + ClinGen
#     clingen_with_vcf = annotated_df.merge(
#         read_clin,
#         left_on="SYMBOL",
#         right_on="Gene Symbol",
#         how="left",
#         suffixes=('_vcf','_clingen'),
#         validate='many_to_one'
#     ).drop(columns=['Gene ID','cytoBand','Genomic Location','Date Last Evaluated'], errors='ignore') \
#      .drop_duplicates()

#     # 4) Read & process VCF
#     vcf_processed = read_vcf(args.vcf_file)
#     vcf_processed['Zygosity'] = [
#     classify_zygosity(gt, vaf, chrom, args.sample_sex)
#     for gt, vaf, chrom in zip(
#         vcf_processed.get('GT', []),
#         vcf_processed.get('VAF', []),
#         vcf_processed.get('#CHROM', [])
#     )]

#     # 5) Merge everything
#     merged_df = pd.merge(
#         clingen_with_vcf,
#         vcf_processed,
#         left_on='Adjusted_Location',
#         right_on='CHROM_POS',
#         how='left'
#     )
#     merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].replace(r'^\s*$', np.nan, regex=True)
#     merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].str.replace(
#         r'\b(MONDO):\1:', r'\1:', regex=True)
    
#     def clean_cell(s):
#         if pd.isna(s):
#             return s
#         parts = re.split(r'[,\|]+', s)
#         seen = set()
#         out = []
#         for p in parts:
#             p = p.strip()
#             if not p or p in seen:
#                 continue
#             seen.add(p)
#             out.append(p)
#         return "|".join(out)
#     merged_df['ClinVar_CLNDISDB_clean'] = merged_df['ClinVar_CLNDISDB'].apply(clean_cell)

#     pos = merged_df.columns.get_loc('ClinVar_CLNVI')
#     patterns = {
#         'MONDO': r'MONDO:([^|]+)',
#         'MedGen': r'MedGen:([^|]+)',
#         'OMIM':  r'OMIM:(\d+)',
#         'Orpha': r'Orphanet:(\d+)',
#         'HPO':   r'Human_Phenotype_Ontology:(HP:\d+)',
#     }

#     for col, pat in patterns.items():
#         merged_df[col] = (
#             merged_df['ClinVar_CLNDISDB_clean']
#               .str.findall(pat)
#               .apply(lambda x: ",".join(x) if isinstance(x, list) and x else np.nan)
#         )
    
#     for col in reversed(list(patterns.keys())):
#         merged_df.insert(pos, col, merged_df.pop(col))
#     merged_df.drop(columns=['ClinVar_CLNDISDB','ClinVar_CLNDISDB_clean'], inplace=True, errors=True)

#     # 6) Write merged VCF+VEP TSV
#     with open(args.merged_vcf_vep, "w", newline="") as f:
#         merged_df.to_csv(f, sep="\t", index=False, lineterminator="\n")

#     # 7) Subset & rename for final CSV
#     desired_columns = [
#         "Adjusted_Location", "SYMBOL", "Gene", "REF", "ALT", "Zygosity", "Codons",
#         "Existing_variation", "HGVSc", "HGVSp", "LOVD", "Amino_acids",
#         "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "IMPACT", "DP", "GQ",
#         "PL", "QUAL", "FILTER", "PHENOTYPES", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT",
#         "ClinVar_CLNDN", "ClinVar_CLNVI", "MONDO","OMIM","Orpha","MedGen","HPO", "PUBMED",
#         "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score",
#         "Consequence", "MutationTaster_pred",
#         "MutationTaster_score", "SIFT4G_score",
#         "SIFT4G_pred", "REVEL_score", "GERP++_NR", "GERP++_RS",
#         "GERP++_RS_rankscore", "gnomAD4.1_joint_SAS_AF", "RegeneronME_ALL_AF",
#         "RegeneronME_E_ASIA_AF", "RegeneronME_SAS_AF"
#     ]
#     #"LRT_pred", "LRT_score","ada_score", "rf_score", 
#     merged_df= merged_df[desired_columns]

#     # initial rename
#     merged_df = merged_df.rename(columns={
#         'SYMBOL': 'SYMBOL (Gene Name)',
#         'Gene':    'GENE (GENE ID)',
#         'Existing_variation': 'SNPs/Rsid'
#     })

#     # Zygosity label
#     zmap = {'Heterozygous':'HET','Homozygous':'HOM'}
#     merged_df.insert(
#         merged_df.columns.get_loc("Zygosity")+1,
#         "Zygosity_label",
#         merged_df["Zygosity"].map(zmap)
#     )

#     # CLINSIG label
#     clinvar_map = {
#         "likely_benign": "LB",
#         "benign": "Ben",
#         "likely_pathogenic": "LP",
#         "pathogenic": "PAT",
#         "uncertain_significance": "VUS",
#         "conflicting_interpretations_of_pathogenicity": "Conflicting"
#     }
#     def map_clinsig(c):
#         if pd.isna(c): return "NA"
#         pat = r'([^/,]+)'
#         def rep(m):
#             term = m.group(0).strip().lower()
#             return clinvar_map.get(term, term)
#         return re.sub(pat, rep, c)
#     merged_df.insert(
#         merged_df.columns.get_loc('ClinVar_CLNSIG')+1,
#         "CLINSIG_label",
#         merged_df['ClinVar_CLNSIG'].apply(map_clinsig)
#     )

    
#     # ACMG DISEASE ID
#     def merge_prioritized(r):
#         if pd.notna(r['OMIM']):    return f"OMIM:{r['OMIM']}"
#         elif pd.notna(r['Orpha']): return f"Orpha:{r['Orpha']}"
#         elif pd.notna(r['MONDO']): return f"Orpha:{r['MONDO']}"
#         elif pd.notna(r['MedGen']):    return f"MedGen:{r['MedGen']}"
#         return ""
#     merged_df['ACMG DISEASE ID (OMIM/ORPHA ID)'] = merged_df.apply(merge_prioritized, axis=1)

#     # final rename
#     merged_df = merged_df.rename(columns={
#         'Gene': 'GENE (GENE ID)',
#         'SYMBOL': 'SYMBOL (Gene Name)',
#         'Existing_variation': 'SNPs/Rsid',
#         'clinvar_clnsig': 'CLINVAR CLNSIG',
#         'Adjusted_Location': 'LOCATION',
#         'Consequence': 'EFFECT'
#     })

#     # 8) Write final CSV
#     with open(args.output_file, "w", newline="") as f:
#         merged_df.to_csv(f, sep=",", index=False, lineterminator="\n")

#     print(f"Done. Merged TSV: {args.merged_vcf_vep}", file=sys.stderr)
#     print(f"Final CSV:      {args.output_file}",      file=sys.stderr)

# if __name__ == "__main__":
#     main()

######### with ACMG intigration #######
#!/usr/bin/env python3
import argparse
import gzip
import sys
import re

import numpy as np
import pandas as pd
from tqdm import tqdm  # For progress tracking
import time
import requests
def parse_args():
    p = argparse.ArgumentParser(
        description="Merge VEP-annotated table, raw VCF, and ClinGen HI/TS data"
    )
    p.add_argument("--vep-file",      required=True,
                   help="Annotated VEP TSV (with ## headers)")
    p.add_argument("--vcf-file",      required=True,
                   help="Raw VCF (.vcf or .vcf.gz)")
    p.add_argument("--acmg-file",      required=True,
                   help="TSV file (.tsv)")
    p.add_argument("--clingen-file",  required=True,
                   help="ClinGen gene curation TSV")
    p.add_argument("--merged-vcf-vep",required=True,
                   help="Output path for merged VCF+VEP TSV")
    p.add_argument("--output-file",   required=True,
                   help="Final output CSV path")
    p.add_argument("--sample-sex",    choices=["male","female"],
                   default="male", help="Sex for hemizygous calls")
    return p.parse_args()

def read_vep(vep_file):
    # capture the VEP headers
    headers = []
    with open(vep_file, "r") as f:
        for line in f:
            if line.startswith("##"):
                headers.append(line)
            elif line.startswith("#"):
                headers.append(line)
                break
    meta = headers[:-1]
    df = pd.read_csv(
        vep_file, sep="\t",
        skiprows=len(meta),
        low_memory=False,
        encoding="latin1",  # handle special characters
    )
    # strip leading '#' from column names
    df.columns = [c.lstrip("#") for c in df.columns]
    # replace '-' with NaN except in key columns
    exclude = {"Location", "Uploaded_variation", "Allele", "UPLOADED_ALLELE"}
    cols = [c for c in df.columns if c not in exclude]
    df.loc[:, cols] = df.loc[:, cols].replace("-", np.nan)
    # 2. Enhanced PubMed ID fetcher
    def get_pmids_from_litvar(rsid):
        """Fetches PMIDs for a given RSID with robust error handling"""
        if not isinstance(rsid, str) or not rsid.startswith('rs'):
            return "NA"
        
        url = f"https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/litvar%40{rsid}%23%23/publications"
        try:
            time.sleep(0.3)  # Rate limiting protection
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                pmids = [str(p) for p in data.get("pmids", []) if str(p).isdigit()]
                return "|".join(pmids) if pmids else "NA"
            return "NA"
        except Exception as e:
            print(f"Error fetching {rsid}: {str(e)}")
            return "NA"

    # 3. Process Existing_variation column (column index 12)
    rsid_column = df.columns[12]  # Using index 12 as specified
    print(f"Processing {len(df)} variants from column '{rsid_column}'")
    
    # Initialize progress bar
    tqdm.pandas(desc="Fetching PubMed IDs")
    
    # Extract first RSID from comma-separated values and fetch PMIDs
    df['pubmed_id'] = df[rsid_column].progress_apply(
        lambda x: get_pmids_from_litvar(str(x).split(',')[0].strip()) if pd.notna(x) else "NA"
    )

    # 4. Show results summary
    success_count = (df['pubmed_id'] != "NA").sum()
    print(f"\nFound PubMed IDs for {success_count}/{len(df)} variants ({success_count/len(df):.1%})")
    
    #return df

    # adjust Location to CHROM:POS or CHROM:POS-1
    def adjust_location(r):
        try:
            allele = str(r.get("Allele","")).strip()
            loc    = str(r["Location"]).strip()
            uv     = str(r["Uploaded_variation"]).strip()
            chrom, pos = loc.split(":")
            left = int(pos.split("-")[0])
            up   = int(uv.split("_")[1])
            if len(allele)==1 and allele.upper() in "ACGT":
                return loc if "-" not in pos else f"{chrom}:{pos.split('-')[0]}"
            return f"{chrom}:{left-1}" if left != up-1 else f"{chrom}:{left}"
        except:
            return r.get("Location","")
    df["Adjusted_Location"] = df.apply(adjust_location, axis=1)
    return df, meta

def read_acmg(acmg_file):
    df = pd.read_csv(acmg_file, sep="\t", low_memory=False, dtype=str)
    df = df[['chromosome', 'position', 'refAllele', 'altAllele', 'acmgClassification', 'rationale']]
    df = df.rename(columns={'chromosome': '#CHROM', 'position': 'POS', 'refAllele': 'REF', 'altAllele': 'ALT', 'acmgClassification': 'ACMGC_Classification', 'rationale': 'ACMG_Criteria'})
    #print(df)
    return df

def read_clingen(path):
    df = pd.read_csv(path, sep="\t", low_memory=False, dtype=str)
    # Identify PMID columns
    haplo_pm = [c for c in df.columns if c.startswith('Haploinsufficiency PMID')]
    triplo_pm = [c for c in df.columns if c.startswith('Triplosensitivity PMID')]

    # Collapse per-row PMIDs into a semicolon-separated string
    def collapse_pmids(row, cols):
        pmids = []
        for col in cols:
            val = row.get(col)
            if pd.notna(val):
                try:
                    pmids.append(str(int(float(val))))
                except (ValueError, TypeError):
                    pass
        return ';'.join(sorted(set(pmids))) if pmids else np.nan

    df['Haploinsufficiency_PMIDs'] = df.apply(lambda r: collapse_pmids(r, haplo_pm), axis=1)
    df['Triplosensitivity_PMIDs']  = df.apply(lambda r: collapse_pmids(r, triplo_pm), axis=1)

    # Helper to combine semicolon-separated lists across genes
    def combine_series_semis(series):
        items = set()
        for cell in series.dropna():
            items.update(cell.split(';'))
        return ';'.join(sorted(items)) if items else np.nan

    # Aggregate ClinGen table to one row per gene, combining fields
    df_agg = df.groupby('Gene Symbol', as_index=False).agg({
        'Haploinsufficiency Score':     lambda s: pd.to_numeric(s, errors='coerce').max(),
        'Triplosensitivity Score':      lambda s: pd.to_numeric(s, errors='coerce').max(),
        'Haploinsufficiency Disease ID':lambda s: ';'.join(sorted(set(s.dropna()))),
        'Triplosensitivity Disease ID': lambda s: ';'.join(sorted(set(s.dropna()))),
        'Genomic Location':             lambda s: ';'.join(sorted(set(s.dropna()))),
        'Haploinsufficiency_PMIDs':     lambda s: combine_series_semis(s),
        'Triplosensitivity_PMIDs':      lambda s: combine_series_semis(s),
    })
    return df_agg

def read_vcf(path):
    import gzip
    vcf_headers = []
    vcf_lines   = []
    # Use gzip.open for .gz files, open otherwise
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##"):
                vcf_headers.append(line)
            elif line.startswith("#CHROM"):
                vcf_headers.append(line)
                # grab the real column names (including your sample name)
                vcf_columns = line.rstrip().split("\t")
            else:
                vcf_lines.append(line.rstrip("\n").split("\t"))

    # build the DataFrame using the real VCF columns
    #df = pd.DataFrame(rows, columns=header_cols)
    fixed_columns      = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    df = pd.DataFrame(vcf_lines, columns=fixed_columns + ["SAMPLE_VALUES"])
    print(df)

    #df = pd.DataFrame(rows, columns=fixed_columns+[sample_column_name])
    required_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']
    if not all(col in df.columns for col in required_cols):
        print("Error: VCF missing required columns.", file=sys.stderr)
        sys.exit(1)

    df["CHROM_POS"] = df["#CHROM"] + ":" + df["POS"].astype(str)
    format_index = df.columns.get_loc("FORMAT")
    sample_column_name = df.columns[format_index + 1]  # Dynamically detect sample column
    print(sample_column_name)
    #format_fields = df['FORMAT'].iloc[0].split(':')  # Extract field names from FORMAT
    #df["CHROM_POS"] = df["#CHROM"] + ":" + df["POS"]
    # split FORMAT and sample columns
    fmt_splits  = df["FORMAT"].str.split(":", expand=False)
    geno_splits = df[sample_column_name].str.split(":", expand=False)
    print("format_fields")
    print(fmt_splits)
    lengths = fmt_splits.str.len().to_numpy()
    rows_i  = np.repeat(df.index.values, lengths)
    keys    = [k for sub in fmt_splits  for k  in sub]
    vals    = [v for sub in geno_splits for v in sub]
    long    = pd.DataFrame({"row": rows_i, "key": keys, "val": vals})
    sample_df = (long.pivot(index="row", columns="key", values="val")
                   .reindex(index=df.index)
                   .fillna("."))
    format_fields = sample_df.columns.tolist()
    print("format_fields")
    print(format_fields)

    sample_df.columns = format_fields 
    vcf = pd.concat([df.drop(columns=['FORMAT', sample_column_name]), sample_df], axis=1)
    #vcf = pd.concat([df.drop(columns=["FORMAT","SAMPLE_VALUES"]), sample_df], axis=1)
    # depth
    if   'DP' in vcf.columns:
        pass
    elif 'DPI' in vcf.columns:
        vcf['DP'] = vcf['DPI']
    elif 'AD' in vcf.columns:
        ad = vcf['AD'].str.split(",", expand=True).iloc[:,0:2].astype(float)
        vcf['DP'] = ad.sum(axis=1)
    if 'DPI' in vcf.columns:
        vcf.drop(columns=['DPI'], inplace=True)
    # GQ
    if 'GQ' not in vcf.columns:
        def calc_gq(r):
            pl = None
            if 'PL' in r and pd.notna(r['PL']):
                raw = list(map(int, str(r['PL']).split(",")))
                base= min(raw)
                pl  = [p-base for p in raw]
            elif 'GL' in r and pd.notna(r['GL']):
                raw = [round(-10*float(x)) for x in str(r['GL']).split(",")]
                base= min(raw)
                pl  = [p-base for p in raw]
            if pl is None: return np.nan
            gt = str(r['GT']).replace("|","/")
            idx= {'0/0':0,'0/1':1,'1/1':2}.get(gt)
            if idx is None: return np.nan
            return min(p for i,p in enumerate(pl) if i!=idx)
        vcf['GQ'] = vcf.apply(calc_gq, axis=1)
    # PL from GL
    if 'PL' not in vcf.columns and 'GL' in vcf.columns:
        def gl_to_pl(gl):
            raw = [round(-10*float(x)) for x in str(gl).split(",")]
            base= min(raw)
            return ",".join(str(p-base) for p in raw)
        vcf['PL'] = vcf['GL'].apply(gl_to_pl)
    # VAF
    if 'VAF' not in vcf.columns and 'AD' in vcf.columns and 'DP' in vcf.columns:
        def calc_vaf(ad, dp):
            try:
                a0,a1 = map(int, ad.split(","))
                return a1/(a0+a1) if (a0+a1)>0 else np.nan
            except:
                return np.nan
        vcf['VAF'] = vcf.apply(lambda r: calc_vaf(r['AD'], r['DP']), axis=1)
    # coerce numeric
    for c in ('DP','GQ','VAF'):
        if c in vcf.columns:
            vcf[c] = pd.to_numeric(vcf[c], errors="coerce")
    return vcf

def classify_zygosity(gt, vaf, chrom, sex):
    alleles = re.split(r"[\/|]", str(gt)) if pd.notna(gt) else []
    ref = alleles.count("0")
    alt = len([a for a in alleles if a not in {"0",".",""}])
    ploidy = len(alleles)
    hap = (sex=="male" and chrom.replace("chr","") in {"X","Y"})
    if hap:
        if alt>=1: return "Hemizygous"
        if ref>=1: return "Hemizygous"
    else:
        if alt==0 and ref>=1: return "Homozygous"
        if alt==ploidy and ploidy>0: return "Homozygous"
        if alt>=1 and ref>=1: return "Heterozygous"
    if pd.notna(vaf):
        if vaf>=0.85: return "Homozygous"
        if vaf<0.15:  return "Homozygous"
        if 0.35<=vaf<=0.65: return "Heterozygous"
    return "Unclassified"

def main():
    args = parse_args()

    # 1) Read VEP
    annotated_df, vep_meta = read_vep(args.vep_file)

    # ) Read ACMG tsv file
    acmg_f = read_acmg(args.acmg_file)

    # 2) Read ClinGen
    read_clin = read_clingen(args.clingen_file)

    # 3) Merge VEP + ClinGen
    clingen_with_vcf = annotated_df.merge(
        read_clin,
        left_on="SYMBOL",
        right_on="Gene Symbol",
        how="left",
        suffixes=('_vcf','_clingen'),
        validate='many_to_one'
    ).drop(columns=['Gene ID','cytoBand','Genomic Location','Date Last Evaluated'], errors='ignore') \
     .drop_duplicates()

    # 4) Read & process VCF
    vcf_processed1 = read_vcf(args.vcf_file)
    vcf_processed1['Zygosity'] = [
    classify_zygosity(gt, vaf, chrom, args.sample_sex)
    for gt, vaf, chrom in zip(
        vcf_processed1.get('GT', []),
        vcf_processed1.get('VAF', []),
        vcf_processed1.get('#CHROM', [])
    )]
    #print(vcf_processed1)

    # Merge DataFrames on the key columns
    vcf_processed = pd.merge(
        vcf_processed1,
        acmg_f[['#CHROM', 'POS', 'REF', 'ALT', 'ACMGC_Classification', 'ACMG_Criteria']],
        on=['#CHROM', 'POS', 'REF', 'ALT'],
        how='left'  # Use 'inner' if you only want rows with matches
    )
    #print(vcf_processed)
    # If you want to keep only rows where a match was found (drop NaN in ACMG columns)
    #vcf_processed = vcf_processed.dropna(subset=['ACMGC_Classification', 'ACMG_Criteria'])

    # 5) Merge everything
    merged_df = pd.merge(
        clingen_with_vcf,
        vcf_processed,
        left_on='Adjusted_Location',
        right_on='CHROM_POS',
        how='left'
    )
    merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].replace(r'^\s*$', np.nan, regex=True)
    merged_df['ClinVar_CLNDISDB'] = merged_df['ClinVar_CLNDISDB'].str.replace(
        r'\b(MONDO):\1:', r'\1:', regex=True)
    
    def clean_cell(s):
        if pd.isna(s):
            return s
        parts = re.split(r'[,\|]+', s)
        seen = set()
        out = []
        for p in parts:
            p = p.strip()
            if not p or p in seen:
                continue
            seen.add(p)
            out.append(p)
        return "|".join(out)
    merged_df['ClinVar_CLNDISDB_clean'] = merged_df['ClinVar_CLNDISDB'].apply(clean_cell)

    pos = merged_df.columns.get_loc('ClinVar_CLNVI')
    patterns = {
        'MONDO': r'MONDO:([^|]+)',
        'MedGen': r'MedGen:([^|]+)',
        'OMIM':  r'OMIM:(\d+)',
        'Orpha': r'Orphanet:(\d+)',
        'HPO':   r'Human_Phenotype_Ontology:(HP:\d+)',
    }

    for col, pat in patterns.items():
        merged_df[col] = (
            merged_df['ClinVar_CLNDISDB_clean']
              .str.findall(pat)
              .apply(lambda x: ",".join(x) if isinstance(x, list) and x else np.nan)
        )
    
    for col in reversed(list(patterns.keys())):
        merged_df.insert(pos, col, merged_df.pop(col))
    merged_df.drop(columns=['ClinVar_CLNDISDB','ClinVar_CLNDISDB_clean'], inplace=True, errors=True)

    # 6) Write merged VCF+VEP TSV
    with open(args.merged_vcf_vep, "w", newline="") as f:
        merged_df.to_csv(f, sep="\t", index=False, lineterminator="\n")

    # 7) Subset & rename for final CSV
    desired_columns = [
        "Adjusted_Location", "SYMBOL", "Gene", "REF", "ALT", "Zygosity", "Codons",
        "Existing_variation", "HGVSc", "HGVSp", "LOVD", "Amino_acids",
        "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "IMPACT", "DP", "GQ",
        "PL", "QUAL", "FILTER", "PHENOTYPES", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT",
        "ClinVar_CLNDN", "ClinVar_CLNVI", "MONDO","OMIM","Orpha","MedGen","HPO", "PUBMED",
        "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score",
        "Consequence", "MutationTaster_pred",
        "MutationTaster_score", "SIFT4G_score",
        "SIFT4G_pred", "REVEL_score", "GERP++_NR", "GERP++_RS",
        "GERP++_RS_rankscore", "gnomAD4.1_joint_SAS_AF", "RegeneronME_ALL_AF",
        "RegeneronME_E_ASIA_AF", "RegeneronME_SAS_AF", "ACMGC_Classification", "ACMG_Criteria", "pubmed_id"
    ]
    #"LRT_pred", "LRT_score","ada_score", "rf_score", 
    merged_df= merged_df[desired_columns]

    # initial rename
    merged_df = merged_df.rename(columns={
        'SYMBOL': 'SYMBOL (Gene Name)',
        'Gene':    'GENE (GENE ID)',
        'Existing_variation': 'SNPs/Rsid'
    })

    # Zygosity label
    zmap = {'Heterozygous':'HET','Homozygous':'HOM'}
    merged_df.insert(
        merged_df.columns.get_loc("Zygosity")+1,
        "Zygosity_label",
        merged_df["Zygosity"].map(zmap)
    )

    # CLINSIG label
    clinvar_map = {
        "likely_benign": "LB",
        "benign": "Ben",
        "likely_pathogenic": "LP",
        "pathogenic": "PAT",
        "uncertain_significance": "VUS",
        "conflicting_interpretations_of_pathogenicity": "Conflicting"
    }
    def map_clinsig(c):
        if pd.isna(c): return "NA"
        pat = r'([^/,]+)'
        def rep(m):
            term = m.group(0).strip().lower()
            return clinvar_map.get(term, term)
        return re.sub(pat, rep, c)
    merged_df.insert(
        merged_df.columns.get_loc('ClinVar_CLNSIG')+1,
        "CLINSIG_label",
        merged_df['ClinVar_CLNSIG'].apply(map_clinsig)
    )

    
    # ACMG DISEASE ID
    def merge_prioritized(r):
        if pd.notna(r['OMIM']):    return f"OMIM:{r['OMIM']}"
        elif pd.notna(r['Orpha']): return f"Orpha:{r['Orpha']}"
        elif pd.notna(r['MONDO']): return f"Orpha:{r['MONDO']}"
        elif pd.notna(r['MedGen']):    return f"MedGen:{r['MedGen']}"
        return ""
    merged_df['ACMG DISEASE ID (OMIM/ORPHA ID)'] = merged_df.apply(merge_prioritized, axis=1)

    # final rename
    merged_df = merged_df.rename(columns={
        'Gene': 'GENE (GENE ID)',
        'SYMBOL': 'SYMBOL (Gene Name)',
        'Existing_variation': 'SNPs/Rsid',
        'clinvar_clnsig': 'CLINVAR CLNSIG',
        'Adjusted_Location': 'LOCATION',
        'Consequence': 'EFFECT'
    })

    # 8) Write final CSV
    with open(args.output_file, "w", newline="") as f:
        merged_df.to_csv(f, sep=",", index=False, lineterminator="\n")

    print(f"Done. Merged TSV: {args.merged_vcf_vep}", file=sys.stderr)
    print(f"Final CSV:      {args.output_file}",      file=sys.stderr)

if __name__ == "__main__":
    main()
