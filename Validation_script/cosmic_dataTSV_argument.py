import pandas as pd
import argparse

def merge_cosmic_data(cosmic_file, input_file, output_file):
    # Load the CosmicCoding TSV file and drop duplicates on the matching columns
    cosmic_df = pd.read_csv(cosmic_file, sep='\t', usecols=['CHROM', 'POS', 'REF', 'ALT', 'COSMIC_ID', 'COSMIC_STATUS']).drop_duplicates(subset=['CHROM', 'POS', 'REF', 'ALT'])
    
    # Load your input TSV file with the variants
    tsv_df = pd.read_csv(input_file, sep='\t')
    
    # Merge the two DataFrames based on the matching columns, only bringing COSMIC_ID and COSMIC_STATUS
    merged_df = pd.merge(
        tsv_df,
        cosmic_df[['CHROM', 'POS', 'REF', 'ALT', 'COSMIC_ID', 'COSMIC_STATUS']], 
        how='left', 
        left_on=['Chr', 'Position', 'Ref', 'Alt'], 
        right_on=['CHROM', 'POS', 'REF', 'ALT']
    )
    
    # Drop the matching columns from CosmicCoding after merge to keep only COSMIC_ID and COSMIC_STATUS
    merged_df = merged_df.drop(columns=['CHROM', 'POS', 'REF', 'ALT'])
    
    # Save the merged data back to a TSV file
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"Output saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Merge CosmicCoding data with variant data.')
    
    # Add arguments for cosmic file, input file, and output file
    parser.add_argument('--cosmic_file', type=str, help='Path to the CosmicCoding TSV file', required=True)
    parser.add_argument('--input_file', type=str, help='Path to the input TSV file containing variants', required=True)
    parser.add_argument('--output_file', type=str, help='Path to save the merged output TSV file', required=True)
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the function with the arguments
    merge_cosmic_data(args.cosmic_file, args.input_file, args.output_file)

if __name__ == '__main__':
    main()
