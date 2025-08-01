import argparse
import os
import pandas as pd

def load_vcf(vcf_file):
    # Read the VCF file into a pandas DataFrame, skipping lines starting with '##'
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    
    # Set column names based on standard VCF columns and two sample columns (tumor and normal)
    vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL']
    
    return vcf_df

def load_tsv(tsv_file):
    # Read the TSV file into a pandas DataFrame
    tsv_df = pd.read_csv(tsv_file, sep=',')
    
    return tsv_df

def extract_tlod(info_field):
    # Check if the info_field is a string (i.e., not NaN)
    if pd.isna(info_field):
        return None
    # Split the INFO field on ';' and search for the TLOD value
    info_parts = info_field.split(';')
    for part in info_parts:
        if part.startswith('TLOD='):
            return part.split('=')[1]
    return None

def split_sample_field(sample_field):
    # Split the SAMPLE column into subfields and handle missing fields
    sample_values = sample_field.split(':') if pd.notna(sample_field) else []
    
    # Define the expected number of columns (excluding 'PGT' and 'PID')
    expected_columns = 8  # This includes GT, AD, AF, DP, F1R2, F2R1, FAD, SB
    
    # Pad with NaN if fewer fields are present
    sample_values.extend([None] * (expected_columns - len(sample_values)))
    
    return sample_values[:expected_columns]

def merge_vcf_tsv(vcf_file, tsv_file, output_file):
    # Extract the initial Sample ID from filenames
    vcf_sample_id = os.path.basename(vcf_file).split("_")[0]
    tsv_sample_id = os.path.basename(tsv_file).split("_")[0]

    # Check if the Sample IDs match
    if vcf_sample_id != tsv_sample_id:
        print(f"Sample IDs do not match: {vcf_sample_id} vs {tsv_sample_id}")
        return

    # Load VCF and TSV files into DataFrames
    vcf_df = load_vcf(vcf_file)
    tsv_df = load_tsv(tsv_file)
    
    # Merge the VCF and TSV files on matching columns
    merged_df = pd.merge(
        tsv_df,
        vcf_df,
        left_on=['Chr', 'Position', 'Ref', 'Alt'],
        right_on=['#CHROM', 'POS', 'REF', 'ALT'],
        how='left'
    )
    
    # Split the TUMOR and NORMAL columns into subfields and handle missing values
    tumor_fields = merged_df['TUMOR'].apply(split_sample_field)
    normal_fields = merged_df['NORMAL'].apply(split_sample_field)
    
    # Create new columns from the sample fields
    sample_columns = ['GT', 'AD', 'AF', 'DP', 'F1R2', 'F2R1', 'FAD', 'SB']
    tumor_df = pd.DataFrame(tumor_fields.tolist(), columns=[f'TUMOR_{col}' for col in sample_columns])
    normal_df = pd.DataFrame(normal_fields.tolist(), columns=[f'NORMAL_{col}' for col in sample_columns])
    
    # Concatenate the new columns with the merged DataFrame
    merged_df = pd.concat([merged_df, tumor_df, normal_df], axis=1)
    
    # Extract the TLOD value from the INFO column
    merged_df['TLOD'] = merged_df['INFO'].apply(extract_tlod)
    
    # Save the merged DataFrame to a TSV file
    merged_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # Argument parser for command line input
    parser = argparse.ArgumentParser(description='Merge VCF and TSV files based on matching columns and Sample ID in filenames')
    
    parser.add_argument('--vcf', required=True, help='Input VCF file with tumor-normal samples')
    parser.add_argument('--tsv', required=True, help='Input TSV file')
    parser.add_argument('--output', required=True, help='Output merged TSV file')
    
    args = parser.parse_args()
    
    # Call the merge function
    merge_vcf_tsv(args.vcf, args.tsv, args.output)
