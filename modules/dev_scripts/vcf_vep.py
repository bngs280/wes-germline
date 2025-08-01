##procces_vcf_vep and merge
import pandas as pdimport sys
import gzip
import numpy as np
import re

def main():
    # Read input file paths from arguments
    vep_file = sys.argv[1]      # Annotated VEP output file
    vcf_file = sys.argv[2]      # Raw VCF file
    output_file = sys.argv[3]   # Final output file path
    merged_vcf_vep = sys.argv[4] # save vcf + vep file

    # ---------------------- Read VEP Output ----------------------
    vep_headers = []
    with open(vep_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                vep_headers.append(line)
            elif line.startswith('#'):
                vep_headers.append(line)
                break

    vep_metadata_headers = vep_headers[:-1]
    # Read annotated VEP file, skipping headers
    annotated_df = pd.read_csv(vep_file, sep='\t', skiprows=len(vep_metadata_headers), low_memory=False)
    annotated_df.columns = [col.lstrip('#') for col in annotated_df.columns]
    
    # List the columns that you want to exclude from replacement
    cols_to_exclude = ["Location", "#Uploaded_variation", "Allele", "UPLOADED_ALLELE"]

    # Create a list of columns where you want to replace '-' with np.nan
    cols_to_replace = [col for col in annotated_df.columns if col not in cols_to_exclude]

    # Replace '-' with np.nan only in those columns.
    annotated_df[cols_to_replace] = annotated_df[cols_to_replace].replace('-', np.nan)

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

    # ---------------------- Read Raw VCF File ----------------------
    vcf_headers = []
    vcf_lines = []

    # Read raw VCF file, preserving headers
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                vcf_headers.append(line)  # Preserve metadata headers
            elif line.startswith('#CHROM'):
                vcf_headers.append(line)  # Preserve column headers
                vcf_columns = line.strip().split('\t')  # Extract column names
            else:
                vcf_lines.append(line.strip().split('\t'))

    fixed_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    # Convert VCF data to DataFrame
    vcf_df = pd.DataFrame(vcf_lines, columns=fixed_columns + ["SAMPLE_VALUES"])

    # ---------------------- Process VCF Data ----------------------
    # Ensure necessary columns exist
    required_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'FORMAT']
    if not all(col in vcf_df.columns for col in required_cols):
        print("Error: VCF missing required columns.", file=sys.stderr)
        sys.exit(1)

    vcf_df["CHROM_POS"] = vcf_df["#CHROM"] + ":" + vcf_df["POS"].astype(str)

    # Detect the sample column dynamically (it is always after FORMAT)
    format_index = vcf_df.columns.get_loc("FORMAT")
    sample_column_name = vcf_df.columns[format_index + 1]  # Dynamically detect sample column
    print(sample_column_name)

    # Vectorized split of FORMAT and sample columns into lists
    fmt_splits  = vcf_df['FORMAT'].str.split(':')
    geno_splits = vcf_df[sample_column_name].str.split(':')

    # Compute how many fields in each row
    lengths = fmt_splits.str.len().to_numpy()

    # Build a flat index of row-numbers repeated for each field
    rows = np.repeat(vcf_df.index.to_numpy(), lengths)

    # Flatten all keys and values
    keys   = [k   for sub in fmt_splits  for k   in sub]
    values = [v   for sub in geno_splits for v   in sub]

    # Make a long DataFrame and pivot it into wide form
    df_long     = pd.DataFrame({'row': rows, 'key': keys, 'val': values})
    sample_data = (df_long
                   .pivot(index='row', columns='key', values='val')
                   .reindex(index=vcf_df.index)      # ensure same order
                   .fillna('.')                       # or leave as NaN
                  )

    # Get your final list of FORMAT fields
    format_fields = sample_data.columns.tolist()

    sample_data.columns = format_fields  # Assign correct column names
    print("sample_data.columns")
    print(sample_data.columns)
    
    # Merge processed sample data back with the fixed columns
    vcf_processed = pd.concat([vcf_df.drop(columns=['FORMAT', sample_column_name]), sample_data], axis=1)

    # ---------------------  Calculate DP if missing -----------------------
    if   'DP'  in vcf_processed.columns:
        depth_col = 'DP'
    elif 'DPI' in vcf_processed.columns:
        depth_col = 'DPI'
    else:
        depth_col = None

    if depth_col == 'DPI':                       # copy Isaac/Strelka depth
        vcf_processed['DP'] = vcf_processed['DPI']
    elif depth_col is None and 'AD' in vcf_processed.columns:    # derive
        vcf_processed['DP'] = (
            vcf_processed['AD']
            .str.split(',', expand=True)[[0, 1]]
            .astype(float).sum(axis=1)
        )

    # Remove DPI column if it exists
    vcf_processed.drop(columns=[c for c in ('DPI',) if c in vcf_processed.columns],
                       inplace=True)

    # ---------------------- Calculate PL if missing ----------------------
    def gl_to_pl(gl):
        if pd.isna(gl):
            return np.nan
        raw = [round(-10 * float(g)) for g in str(gl).split(',')]
        return ','.join(str(p - min(raw)) for p in raw)

    if 'PL' not in vcf_processed.columns and 'GL' in vcf_processed.columns:
        print("Adding PL value")
        vcf_processed['PL'] = vcf_processed['GL'].apply(gl_to_pl)

    # ---------------------- Calculate GQ if Missing -----------------------
    if 'GQ' not in vcf_processed.columns:
        def calc_gq(row):
            """
            Reconstruct GQ from PL (preferred) or GL.
            GQ = smallest PL among all genotypes other than the chosen GT.
            """
            pl = None

            # 1) Try PL first – already on Phred scale
            if 'PL' in format_fields and pd.notna(row.get('PL')):
                pl_raw = list(map(int, str(row['PL']).split(',')))
                min_pl = min(pl_raw)
                pl     = [p - min_pl for p in pl_raw]

            # 2) Otherwise derive PL from GL
            elif 'GL' in format_fields and pd.notna(row.get('GL')):
                gl      = list(map(float, str(row['GL']).split(',')))
                pl_raw  = [round(-10 * g) for g in gl]       # Phred-scale
                min_pl  = min(pl_raw)
                pl      = [p - min_pl for p in pl_raw]       # re-zero so best genotype = 0

            if pl is None:
                return np.nan   # can't compute GQ without likelihoods

            gt = str(row['GT']).replace('|', '/')
            # biallelic ordering per VCF spec: 0/0, 0/1, 1/1
            gt_index = {'0/0': 0, '0/1': 1, '1/1': 2}.get(gt)
            if gt_index is None or gt_index >= len(pl):
                return np.nan

            # GQ = min PL of all non-best genotypes
            return min(p for i, p in enumerate(pl) if i != gt_index)

        vcf_processed['GQ'] = vcf_processed.apply(calc_gq, axis=1)

    # ---------------------- Calculate VAF if Missing ----------------------
    if 'VAF' not in vcf_processed.columns and 'AD' in vcf_processed.columns and 'DP' in vcf_processed.columns:
        def calculate_vaf(ad, dp):
            try:
                ad_values = list(map(int, ad.split(',')))
                if len(ad_values) > 1:
                    return float(ad_values[1]) / float(sum(ad_values))
            except:
                return np.nan
            return np.nan

        vcf_processed['VAF'] = vcf_processed.apply(lambda row: calculate_vaf(row['AD'], row['DP']), axis=1)
    
    # Ensure VAF is numeric
    if 'VAF' in vcf_processed.columns:
        vcf_processed['VAF'] = pd.to_numeric(vcf_processed['VAF'], errors='coerce')

    # ---------------------- Classify Zygosity ----------------------
    def classify_zygosity(gt, vaf):
        if pd.isna(vaf):  # Handle missing VAF values
            return 'Unclassified'
        try:
            vaf = float(vaf)  # Ensure vaf is float
            if gt == '1/1' or (pd.isna(gt) and vaf >= 0.85):
                return 'Homozygous'
            elif gt == '0/1' or (pd.isna(gt) and 0.35 <= vaf <= 0.65):
                return 'Heterozygous'
            elif gt == '0/0' or (pd.isna(gt) and vaf < 0.01):
                return 'Reference'
            elif vaf < 0.35:
                return 'Low VAF'
            else:
                return 'Unclassified'
        except (ValueError, TypeError):
            return 'Unclassified'

    vcf_processed['Zygosity'] = vcf_processed.apply(lambda row: classify_zygosity(row['GT'], row['VAF']), axis=1)
    
    # Convert numeric fields
    numeric_fields = ['GQ', 'DP', 'VAF']
    for field in numeric_fields:
        if field in vcf_processed.columns:
            vcf_processed[field] = pd.to_numeric(vcf_processed[field], errors='coerce')

    # ---------------------- Merge VEP and VCF Data ----------------------
    vcf_processed["#CHROM_POS"] = vcf_processed["#CHROM"] + ":" + vcf_processed["POS"].astype(str)
    print("Columns in annotated_df:", annotated_df.columns.tolist())
    print("Columns in vcf_processed:", vcf_processed.columns.tolist())

    # Merge VEP annotations with processed VCF data
    merged_df = pd.merge(annotated_df, vcf_processed, left_on=['Adjusted_Location'], right_on=['#CHROM_POS'], how='left')

    # Save merged VCF+VEP data
    with open(merged_vcf_vep, 'w') as f:
        merged_df.to_csv(f, sep='\t', index=False)

    # Select and rename desired columns
    desired_columns = [
        "Adjusted_Location", "SYMBOL", "Gene", "REF", "ALT", "Zygosity", "Codons", "Existing_variation", 'HGVSc_VEP', 'HGVSp_VEP', "LOVD", "Amino_acids", "MANE_SELECT", "MANE_PLUS_CLINICAL", "SpliceRegion", "IMPACT", "DP", "GQ", "PL","QUAL", "FILTER",
        "PHENOTYPES", "CLIN_SIG", "clinvar_review", "clinvar_Orphanet_id", 'clinvar_MedGen_id', "clinvar_OMIM_id", "PUBMED", "CADD_PHRED", "CADD_RAW", "AlphaMissense_rankscore", "AlphaMissense_score", "Consequence", "LRT_pred", "LRT_score", "MutationTaster_pred", "MutationTaster_score", "ada_score", "rf_score", "SIFT4G_score", "SIFT4G_pred", "REVEL_score", "GERP++_NR" ,"GERP++_RS",
        "GERP++_RS_rankscore", "gnomAD_exomes_AF", "gnomAD_exomes_AFR_AF", "gnomAD_exomes_AMR_AF", "gnomAD_exomes_EAS_AF", "gnomAD_exomes_FIN_AF", "gnomAD_exomes_NFE_AF", "gnomAD_exomes_SAS_AF", "1000Gp3_AF", "1000Gp3_AFR_AF", "1000Gp3_AMR_AF", "1000Gp3_EAS_AF", "1000Gp3_EUR_AF", "1000Gp3_SAS_AF", "ESP6500_AA_AF", "ESP6500_EA_AF", "ExAC_SAS_AF", "SpliceAI_pred"]

    merged_df = merged_df[desired_columns]

    merged_df = merged_df.rename(columns={'HGNC_ID': 'SYMBOL (Gene Name)', 'Gene': 'GENE (GENE ID)', 'Existing_variation': 'SNPS/RSID'})
    
    # Add Zygosity label
    zygosity_mapping = {
        'Heterozygous': 'HET',
        'Homozygous': 'HOM'
    }
    merged_df.insert(merged_df.columns.get_loc("Zygosity") + 1,
            "Zygosity_label",
            merged_df["Zygosity"].map(zygosity_mapping))

    # Add CLINSIG label
    clinvar_mapping = {
        "likely_benign": "LB",
        "benign": "Ben",
        "likely_pathogenic": "LP",
        "pathogenic": "PAT",
        "uncertain_significance": "VUS",
        "conflicting_interpretations_of_pathogenicity": "Conflicting"
    }

    def map_clinsig(clinsig):
        if pd.isna(clinsig):  # Handle missing values
            return "NA"
        pattern = r'([^/,]+)'   # Matches any sequence of characters except / and ,
        def replace_match(match):
            term = match.group(0).strip().lower()  #Extract term and normalize
            return clinvar_mapping.get(term, term)
        mapped_clinsig = re.sub(pattern, replace_match, clinsig)
        return mapped_clinsig

    merged_df.insert(merged_df.columns.get_loc('CLIN_SIG') + 1, "CLINSIG_label", merged_df['CLIN_SIG'].apply(map_clinsig))

    # Add ACMG DISEASE ID
    def merge_prioritized(row):
        # Check for OMIM (highest priority)
        if pd.notnull(row['clinvar_OMIM_id']):
            return f"OMIM:{row['clinvar_OMIM_id']}"
        elif pd.notnull(row['clinvar_Orphanet_id']):
            return f"Orpha:{row['clinvar_Orphanet_id']}"
        elif pd.notnull(row['clinvar_MedGen_id']):
            return f"MedGen:{row['clinvar_MedGen_id']}"
        return ""

    merged_df['ACMG DISEASE ID (OMIM/ORPHA ID)'] = merged_df.apply(merge_prioritized, axis=1)
    
    # Rename additional columns
    merged_df = merged_df.rename(columns={
                'Gene': 'GENE (GENE ID)',
                'SYMBOL': 'SYMBOL (Gene Name)',
                'Existing_variation': 'SNPs/Rsid',
                'CLIN_SIG': 'CLINVAR CLNSIG',
                'Adjusted_Location': 'LOCATION',
                'Consequence': 'EFFECT'
            })

    # ---------------------- Write Final Output ----------------------
    # Get current column names as a list
    cols = merged_df.columns.tolist()
    # Prefix the first column with '#' if it's not already prefixed
    if not cols[0].startswith('#'):
        cols[0] = '#' + cols[0]
    merged_df.columns = cols

    with open(output_file, 'w') as f:
        merged_df.to_csv(f, sep=',', index=False)

    print(f"Final annotated VCF saved at: {output_file}")
    print(f"Merged VCF+VEP data saved at: {merged_vcf_vep}")

if __name__ == "__main__":
    main()
