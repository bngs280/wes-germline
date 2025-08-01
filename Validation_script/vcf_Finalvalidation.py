import os
import argparse
import subprocess
import pysam
import glob
import pandas as pd
import gzip

def validate_vcf(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    if not file_path.endswith('_filtered_PASS.vcf'):
        raise ValueError(f"The file {file_path} is not a valid VCF file.")
    
    if os.path.getsize(file_path) == 0:
        raise ValueError(f"The file {file_path} is empty.")

    with open(file_path, 'rt') as vcf_file:
        lines = vcf_file.readlines()
        
        # Check header lines
        header_lines = [line for line in lines if line.startswith('##')]
        if not header_lines:
            raise ValueError(f"The VCF file {file_path} does not contain header lines starting with '##'.")
        
        # Check column headers
        column_line = [line for line in lines if line.startswith('#CHROM')]
        if not column_line:
            raise ValueError(f"The VCF file {file_path} does not contain the required column headers.")
        
        required_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        column_headers = column_line[0].strip().split('\t')
        
        if not all(col in column_headers for col in required_columns):
            raise ValueError(f"The VCF file {file_path} does not contain all the required columns: {', '.join(required_columns)}")
        
        # Check for variant data
        variant_lines = [line for line in lines if not line.startswith('#')]
        if not variant_lines:
            raise ValueError(f"The VCF file {file_path} does not contain variant data.")
    
    print(f"Output VCF file {file_path} passed validation.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate a VCF file.")
    parser.add_argument("--vcf", type=str, required=True, help="Path to the VCF file to validate.")

    args = parser.parse_args()

    try:
        validate_vcf(args.vcf)
    except (FileNotFoundError, ValueError) as e:
        print(f"Validation failed: {str(e)}")
