import os
import argparse
import subprocess
import pysam
import glob
import pandas as pd
import gzip
# Define chromosome lengths for hg19 and hg38

HG19_LENGTHS = {
    'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276,
    'chr5': 180915260, 'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022,
    'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516, 'chr12': 133851895,
    'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 'chr16': 90354753,
    'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
    'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566, 'chrM': 16571
}

HG38_LENGTHS = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415, 'chrM': 16569
}

def validate_vcf(file_path, user_genome_version):
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
        
        # Check genome version (hg19 vs hg38)
        contig_lines = [line for line in lines if line.startswith('##contig')]
        vcf_genome_version = None
        genome_lengths = HG38_LENGTHS if user_genome_version == 'hg38' else HG19_LENGTHS
        for contig_line in contig_lines:
            # Extract the chromosome name and length from contig lines
            if "length=" in contig_line:
                contig_info = contig_line.split('<')[1].strip('>')
                contig_info_dict = {k: v for k, v in (item.split('=') for item in contig_info.split(','))}
                contig_name = contig_info_dict.get('ID')
                # Clean length to remove any '>' or newlines before converting
                length_str = contig_info_dict.get('length').strip('>\n')
                length = int(length_str)  # Ensure this is an integer
                # Check if it matches the user-specified genome version
                if contig_name in genome_lengths and genome_lengths[contig_name] == length:
                    vcf_genome_version = user_genome_version
                    break
        if vcf_genome_version == user_genome_version:
            print(f"VCF file {file_path} matches the user-specified genome version: {user_genome_version}.")
        else:
            raise ValueError(f"VCF file {file_path} does not match the user-specified genome version: {user_genome_version}.")
    print(f"Output VCF file {file_path} passed validation.")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate a VCF file.")
    parser.add_argument("--vcf", type=str, required=True, help="Path to the VCF file to validate.")
    parser.add_argument("--genome_version", type=str, required=True, choices=["hg19", "hg38"], help="Specify the genome version (hg19 or hg38).")

    args = parser.parse_args()

    try:
        validate_vcf(args.vcf, args.genome_version)
    except (FileNotFoundError, ValueError) as e:
        print(f"Validation failed: {str(e)}")
