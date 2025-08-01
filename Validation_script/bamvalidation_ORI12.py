import os
import argparse
import pysam

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

def validate_bam(file_path, user_genome_version):
    """
    Validates the given BAM file based on several criteria:
    - The file must exist.
    - The file name must end with '_sorted.bam'.
    - The file must not be empty.
    - The BAM header must contain certain tags: @HD, @SQ, @RG, and @PG.
    - The BAM file must contain aligned reads.
    - The BAM file must match the user-provided genome assembly (hg19 or hg38).
    """
    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")

    # Check if the file is a valid sorted BAM file
    if not file_path.endswith('_sorted.bam'):
        raise ValueError(f"The file {file_path} is not a valid BAM file. It must end with '_sorted.bam'.")

    # Check if the file is empty
    if os.path.getsize(file_path) == 0:
        raise ValueError(f"The file {file_path} is empty.")

    # Open the BAM file using pysam
    bamfile = pysam.AlignmentFile(file_path, "rb")

    # Validate BAM header for required tags
    header = bamfile.header
    required_tags = ['HD', 'SQ', 'RG', 'PG']
    for tag in required_tags:
        if tag not in header:
            raise ValueError(f"The BAM file {file_path} does not contain the required @{tag} tag in the header.")

    # Check if the BAM file contains aligned reads
    if bamfile.count(until_eof=True) == 0:
        raise ValueError(f"The BAM file {file_path} does not contain aligned reads.")

    # Infer genome version by checking @SQ header (sequence dictionary)
    genome_lengths = HG38_LENGTHS if user_genome_version == 'hg38' else HG19_LENGTHS
    bam_genome_version = None

    for sq in header['SQ']:
        chrom = sq['SN']  # Chromosome name
        length = sq['LN']  # Chromosome length
        if chrom in genome_lengths and genome_lengths[chrom] == length:
            bam_genome_version = user_genome_version
            break

    if bam_genome_version == user_genome_version:
        print(f"BAM file {file_path} matches the user-specified genome version: {user_genome_version}.")
    else:
        raise ValueError(f"BAM file {file_path} does not match the user-specified genome version: {user_genome_version}.")

    # Close the BAM file
    bamfile.close()

    print(f"BAM file {file_path} passed validation.")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Validate a BAM file.")
    parser.add_argument('--bam', type=str, required=True, help="Path to the BAM file to validate.")
    parser.add_argument('--assembly_version', type=str, required=True, choices=["hg19", "hg38"], help="Specify the genome version (hg19 or hg38).")

    # Parse the arguments
    args = parser.parse_args()

    try:
        validate_bam(args.bam, args.assembly_version)
    except (FileNotFoundError, ValueError) as e:
        print(f"Validation failed: {str(e)}")
