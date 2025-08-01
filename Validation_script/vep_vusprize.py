import pandas as pd
import argparse

def main(input_file, output_file):
    # Read the TSV file into a DataFrame
    # Skip lines until the line starting with '#Uploaded_variation' is found
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('#Uploaded_variation'):
                header_line = line.strip().split('\t')
                break

    # Read the data into a DataFrame, skipping lines before the header
    df = pd.read_csv(input_file, sep='\t', comment='#', header=None, names=header_line)

    # Rename the column from 'Eigen-phred_coding' to 'Eigen-pred_coding'
    df = df.rename(columns={'Eigen-phred_coding': 'Eigen-pred_coding'})
    #df = df.drop()
    df = df.drop_duplicates(subset=['#Uploaded_variation'])

    # Save the rearranged DataFrame to the output file
    df.to_csv(output_file, sep='\t', index=False)

    print(f"File saved to {output_file}")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process a TSV file and extract/reorder columns.")
    parser.add_argument('--input_file', type=str, help="Path to the input TSV file.")
    parser.add_argument('--output_file', type=str, help="Path to save the output TSV file.")

    # Parse arguments
    args = parser.parse_args()

    # Call the main function
    main(args.input_file, args.output_file)
