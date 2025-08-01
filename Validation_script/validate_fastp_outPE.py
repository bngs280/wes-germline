import os
import gzip
import argparse
import sys  # Add sys module to handle exits


def validate_fastq(file_path):
    """Validate a FASTQ file."""
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        sys.exit(1)
    
    if not file_path.endswith('.fastq.gz'):
        print(f"Error: The file {file_path} is not a valid fastq.gz file.")
        sys.exit(1)
    
    if os.path.getsize(file_path) == 0:
        print(f"Error: The file {file_path} is empty.")
        sys.exit(1)
    
    with gzip.open(file_path, 'rt') as read_obj:
        first_char = read_obj.read(1)
    
    if not first_char:
        print(f"Error: The file {file_path} is empty.")
        sys.exit(1)
    
    with gzip.open(file_path, 'rt') as f:
        line_count = 0
        for line in f:
            line_count += 1
            if line_count % 4 == 1 and not line.startswith('@'):
                print(
                    f"Error: File {file_path} is not in proper FASTQ format at read {line_count // 4 + 1}. "
                    f"Line does not start with '@': {line.strip()}"
                )
                sys.exit(1)
            if line_count % 4 == 2 and line.strip() == '':
                print(
                    f"Error: File {file_path} is not in proper FASTQ format at read {line_count // 4 + 1}. "
                    "Sequence line is empty."
                )
                sys.exit(1)
            if line_count % 4 == 3 and not line.startswith('+'):
                print(
                    f"Error: File {file_path} is not in proper FASTQ format at read {line_count // 4 + 1}. "
                    f"Line does not start with '+': {line.strip()}"
                )
                sys.exit(1)
            if line_count % 4 == 0 and line.strip() == '':
                print(
                    f"Error: File {file_path} is not in proper FASTQ format at read {line_count // 4}. "
                    "Quality line is empty."
                )
                sys.exit(1)
        
        if line_count % 4 != 0:
            print(f"Error: File {file_path} is not in proper FASTQ format. It has {line_count % 4} extra lines.")
            sys.exit(1)
    
    return line_count // 4


def check_paired_files(file_r1, file_r2):
    """Validate and compare paired FASTQ files."""
    reads_r1 = validate_fastq(file_r1)
    reads_r2 = validate_fastq(file_r2)
    
    if reads_r1 != reads_r2:
        print(f"Error: The number of reads in {file_r1} and {file_r2} are not equal.")
        sys.exit(1)
    
    print(f"Validation successful for input files {file_r1} and {file_r2}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate and compare paired FASTQ files.")
    parser.add_argument("--input1", required=True, help="Path to the R1 FASTQ file.")
    parser.add_argument("--input2", required=True, help="Path to the R2 FASTQ file.")
    
    args = parser.parse_args()
    
    try:
        check_paired_files(args.input1, args.input2)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

