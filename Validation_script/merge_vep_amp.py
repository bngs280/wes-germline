import pandas as pd
import argparse
import os

# Function to retain only the first value from comma-separated entries
def retain_first_value(row, column_name):
    return row[column_name].split(",")[0] if pd.notnull(row[column_name]) else row[column_name]

# Main function to merge VEP and CancerVar files based on Sample ID in filenames
def merge_annotations(vep_file, amp_file, output_file):
    # Extract the initial Sample ID from filenames
    vep_sample_id = os.path.basename(vep_file).split("_")[0]
    amp_sample_id = os.path.basename(amp_file).split(".")[0]

    # Check if the Sample IDs match
    if vep_sample_id != amp_sample_id:
        print(f"Sample IDs do not match: {vep_sample_id} vs {amp_sample_id}")
        return

    # Load VEP and AMP (CancerVar) files
    vep = pd.read_csv(vep_file, sep="\t")
    cancevar = pd.read_csv(amp_file, sep="\t")

    # Split '#Uploaded_variation' into components
    vep[['Chr', 'Position', 'Ref_Alt']] = vep['#Uploaded_variation'].str.split('_', expand=True)

    # Further split the 'Ref_Alt' column based on '/'
    vep[['Ref', 'Alt']] = vep['Ref_Alt'].str.split('/', expand=True, n=1)

    # Drop the original 'Ref_Alt' column
    vep.drop('Ref_Alt', axis=1, inplace=True)

    # Convert Position column to integer type
    vep['Position'] = vep['Position'].astype(int)

    # Rename columns in the AMP (CancerVar) file for consistency
    cancevar.rename(columns={'#Chr': 'Chr', 'Start': 'Position', ' CancerVar: CancerVar and Evidence ': 'Classification'}, inplace=True)

    # Add "chr" prefix to the Chr column in AMP (CancerVar)
    cancevar['Chr'] = 'chr' + cancevar['Chr'].astype(str)

    # Merge VEP with AMP (CancerVar) on Chr, Position, Ref, and Alt columns
    merged_df = vep.merge(cancevar[['Chr', 'Position', 'Ref', 'Alt', 'Classification']],
                          on=['Chr', 'Position', 'Ref', 'Alt'],
                          how='left')

    # Add a new column 'Classification' in the VEP file where CancerVar matches
    merged_df['Classification'] = merged_df['Classification']

    # Extract the classification part between '#' and 'EVS' using regular expression
    merged_df['AMP_ACMG_Classification'] = merged_df['Classification'].str.extract(r'#(.*?)\s*EVS')

    # Drop duplicate rows based on '#Uploaded_variation'
    merged_df = merged_df.drop_duplicates(subset=['#Uploaded_variation'])

    # Apply the function to relevant columns to retain the first value from comma-separated entries
    merged_df["HGVSc_ANNOVAR"] = merged_df.apply(retain_first_value, column_name="HGVSc_ANNOVAR", axis=1)
    merged_df["HGVSp_ANNOVAR"] = merged_df.apply(retain_first_value, column_name="HGVSp_ANNOVAR", axis=1)

    # Save the final merged result to the output file
    merged_df.to_csv(output_file, index=False)
    print(f"File saved to {output_file}")

if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description='Merge VEP and AMP (CancerVar) annotations based on Sample ID.')
    parser.add_argument('--vep', required=True, help='Path to the VEP annotation file')
    parser.add_argument('--amp', required=True, help='Path to the AMP (CancerVar) file')
    parser.add_argument('--output', required=True, help='Path to the output merged file')

    # Parse the arguments
    args = parser.parse_args()

    # Call the merge function
    merge_annotations(args.vep, args.amp, args.output)
