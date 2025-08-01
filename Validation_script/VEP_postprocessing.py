import argparse
import pandas as pd

# Function to read and process the VCF file
def process_vep_file(input_file, output_file):
    # Initialize an empty list to store the lines
    vcf_data = []

    # Open the VCF file and read lines
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                if line.startswith('#Uploaded_variation'):
                    # Extract column names from the header line
                    column_names = line.strip().split('\t')
                continue  # Skip all header lines including the column header
            vcf_data.append(line.strip().split('\t'))

    # Create a DataFrame from the collected data with extracted column names
    vep = pd.DataFrame(vcf_data, columns=column_names)

    # Save the DataFrame to the specified output file
    vep.to_csv(output_file, sep='\t', index=False)
    print(f"File has been processed and saved to: {output_file}")

# Main function to handle command-line arguments
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process VEP annotated VCF files.")
    parser.add_argument('--input', required=True, help="Input VEP annotated VCF file")
    parser.add_argument('--output', required=True, help="Output TSV file")
    
    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Process the input VEP file and generate the output
    process_vep_file(args.input, args.output)

# Execute the main function if this script is run
if __name__ == "__main__":
    main()
