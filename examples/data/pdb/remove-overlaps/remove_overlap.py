import os
import pandas as pd
import numpy as np

class DataTransformer:
    def __init__(self, input_files, output_file):
        self.input_files = input_files
        self.output_file = output_file

    def couple_triples(self):
        """Read, transform, and write data from input files to output file."""
        transformed_data = []

        for input_file in self.input_files:
            # Read the CSV file into a DataFrame, ignoring the last column
            df = pd.read_csv(input_file, usecols=[0, 1, 2, 3, 4, 5])
            print(f"Data read from {input_file}")

            # Transform the data
            for _, row in df.iterrows():
                protein1, residue1, protein2, residue2, protein3, residue3 = row
                transformed_data.append([protein1, residue1, protein2, residue2])
                transformed_data.append([protein2, residue2, protein3, residue3])
                transformed_data.append([protein1, residue1, protein3, residue3])
            print(f"Data transformation complete for {input_file}")

        # Convert the transformed data back to a DataFrame
        transformed_data = np.array(transformed_data)
        transformed_df = pd.DataFrame(transformed_data, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])

        # Check for duplicate rows
        duplicates = transformed_df[transformed_df.duplicated(keep=False)]
        if not duplicates.empty:
            print("Duplicate rows found:")
            print(duplicates)
        else:
            print("No duplicate rows found.")

        # Remove duplicate rows and keep only unique rows
        unique_df = transformed_df.drop_duplicates()
        num_unique_rows = unique_df.shape[0]
        num_removed_rows = transformed_df.shape[0] - num_unique_rows

        print(f"Number of unique rows: {num_unique_rows}")
        print(f"Number of rows removed: {num_removed_rows}")

        # Ensure the output directory exists
        os.makedirs('full-proteasome-data', exist_ok=True)

        # Write the unique transformed data to a new CSV file
        output_path = os.path.join('full-proteasome-data', self.output_file)
        unique_df.to_csv(output_path, index=False)
        print(f"Transformed data has been written to {output_path}")

    def merge_and_clean_doubles(self, input_file1, input_file2, output_file):
        """Read, merge, clean, and write data from two input files to output file."""
        # Read the CSV files into DataFrames, ignoring the last column
        df1 = pd.read_csv(input_file1, usecols=[0, 1, 2, 3])
        df2 = pd.read_csv(input_file2, usecols=[0, 1, 2, 3])
        print(f"Data read from {input_file1} and {input_file2}")

        # Merge the two DataFrames
        combined_df = pd.concat([df1, df2])
        print("Data merged")

        # Check for duplicate rows
        duplicates = combined_df[combined_df.duplicated(keep=False)]
        if not duplicates.empty:
            print("Duplicate rows found:")
            print(duplicates)
        else:
            print("No duplicate rows found.")

        # Remove duplicate rows and keep only unique rows
        unique_df = combined_df.drop_duplicates()
        num_unique_rows = unique_df.shape[0]
        num_removed_rows = combined_df.shape[0] - num_unique_rows

        print(f"Number of unique rows: {num_unique_rows}")
        print(f"Number of rows removed: {num_removed_rows}")

        # Ensure the output directory exists
        os.makedirs('full-proteasome-data', exist_ok=True)

        # Write the unique transformed data to a new CSV file
        output_path = os.path.join('full-proteasome-data', output_file)
        unique_df.to_csv(output_path, index=False)
        print(f"Cleaned data has been written to {output_path}")

    def find_and_remove_common_rows(self, file1, file2):
        """Find and remove common rows between two CSV files from the first file."""
        # Read the CSV files into DataFrames
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        print(f"Data read from {file1} and {file2}")

        # Ensure columns have the same data type
        df1 = df1.astype(str)
        df2 = df2.astype(str)
        print("df1 : ", df1)
        print("df2 : ", df2)
        
        # Identify common rows
        common_rows = pd.merge(df1, df2, how='inner')
        if not common_rows.empty:
            print("Common rows found:")
            print(common_rows)
        else:
            print("No common rows found.")
        print("Common rows: ", common_rows)

# Example usage
if __name__ == "__main__":
    input_files = ['non_ambiguous_triple_xls.csv', 'ambiguous_triple_xls.csv']
    output_file = 'combined_coupled_triples.csv'

    transformer = DataTransformer(input_files, output_file)
    transformer.couple_triples()

    # New function usage
    input_file1 = 'ambiguous_double_xls.csv'
    input_file2 = 'non_ambiguous_double_xls.csv'
    output_file_doubles = 'cleaned_doubles.csv'

    transformer.merge_and_clean_doubles(input_file1, input_file2, output_file_doubles)

    # Find and remove common rows between the two final dataframes
    final_file1 = os.path.join('full-proteasome-data', 'combined_coupled_triples.csv')
    final_file2 = os.path.join('full-proteasome-data', 'cleaned_doubles.csv')
    transformer.find_and_remove_common_rows(final_file1, final_file2)