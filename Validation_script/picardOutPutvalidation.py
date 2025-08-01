import os
import argparse
import pysam

def validate_bam(file_path):
    """
    Validates the given BAM file based on several criteria:
    - The file must exist.
    - The file name must end with '_sorted.bam'.
    - The file must not be empty.
    - The BAM header must contain certain tags: @HD, @SQ, @RG, and @PG.
    - The BAM file must contain aligned reads.
    """
    # Check if the file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    # Check if the file is a valid sorted BAM file
    if not file_path.endswith('_sorted_md.bam'):
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
    
    # Close the BAM file
    bamfile.close()
    
    print(f"BAM file {file_path} passed validation.")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Validate a BAM file.")
    parser.add_argument('--bam', type=str, required=True, help="Path to the BAM file to validate.")
    
    # Parse the arguments
    args = parser.parse_args()

    try:
        validate_bam(args.bam)
    except (FileNotFoundError, ValueError) as e:
        print(f"Validation failed: {str(e)}")
