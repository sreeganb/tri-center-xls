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
            #print(f"Data read from {input_file}")

            # Transform the data
            for _, row in df.iterrows():
                protein1, residue1, protein2, residue2, protein3, residue3 = row
                transformed_data.append([protein1, residue1, protein2, residue2])
                transformed_data.append([protein2, residue2, protein3, residue3])
                transformed_data.append([protein1, residue1, protein3, residue3])
            #print(f"Data transformation complete for {input_file}")

        # Convert the transformed data back to a DataFrame
        transformed_data = np.array(transformed_data)
        transformed_df = pd.DataFrame(transformed_data, columns=['Protein1', 'Residue1', 'Protein2', 'Residue2'])

        # Check for duplicate rows
        duplicates = transformed_df[transformed_df.duplicated(keep=False)]
        if not duplicates.empty:
            print("Duplicate rows found:")
            #print(duplicates)
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
        #print(f"Transformed data has been written to {output_path}")

    def merge_and_clean_doubles(self, input_file1, input_file2, output_file):
        """Read, merge, clean, and write data from two input files to output file."""
        # Read the CSV files into DataFrames, ignoring the last column
        df1 = pd.read_csv(input_file1, usecols=[0, 1, 2, 3])
        df2 = pd.read_csv(input_file2, usecols=[0, 1, 2, 3])
        #print(f"Data read from {input_file1} and {input_file2}")

        # Merge the two DataFrames
        combined_df = pd.concat([df1, df2])
        #print("Data merged")

        # Check for duplicate rows
        duplicates = combined_df[combined_df.duplicated(keep=False)]
        if not duplicates.empty:
            print("Duplicate rows found:")
            #print(duplicates)
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
        #print(f"Cleaned data has been written to {output_path}")

    def find_and_remove_common_rows(self, file1, file2, output_file1, output_file2):
        """Find and remove common rows between two CSV files from both files."""
        # Read the CSV files into DataFrames
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        #print(f"Data read from {file1} and {file2}")

        # Ensure columns have the same data type
        df1 = df1.astype(str)
        df2 = df2.astype(str)
        #print("df1 : ", df1)
        #print("df2 : ", df2)
        
        # Identify common rows
        common_rows = pd.merge(df1, df2, how='inner')
        if not common_rows.empty:
            print("Common rows found:")
            #print(common_rows)
        else:
            print("No common rows found.")
        print("Common rows: ", common_rows)

        # Remove common rows from both df1 and df2
        df1_unique = df1[~df1.apply(tuple, 1).isin(common_rows.apply(tuple, 1))]
        df2_unique = df2[~df2.apply(tuple, 1).isin(common_rows.apply(tuple, 1))]

        # Ensure the output directory exists
        os.makedirs('full-proteasome-data', exist_ok=True)

        # Write the unique data to new CSV files
        output_path1 = os.path.join('full-proteasome-data', output_file1)
        output_path2 = os.path.join('full-proteasome-data', output_file2)
        df1_unique.to_csv(output_path1, index=False)
        df2_unique.to_csv(output_path2, index=False)
        #print(f"Unique data from {file1} has been written to {output_path1}")
        #print(f"Unique data from {file2} has been written to {output_path2}")

        # Call filter_proteins with the unique DataFrames
        self.filter_proteins(df1_unique, df2_unique)

    def filter_proteins(self, df_triple_links, df_double_links):
        """
        Filters the DataFrames to include only rows where Protein1, Protein2, and Protein3 (for triple links)
        or Protein1 and Protein2 (for double links) are exclusively from the allowed set.

        Args:
            df_triple_links (DataFrame): The DataFrame containing triple links.
            df_double_links (DataFrame): The DataFrame containing double links.

        Returns:
            tuple: A tuple containing the filtered triple links DataFrame and double links DataFrame.
        """
        allowed_proteins = {'Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn2'}

        def filter_dataframe(df, columns):
            if all(col in df.columns for col in columns):
                return df[df[columns].apply(lambda x: set(x).issubset(allowed_proteins), axis=1)]
            else:
                print(f"Warning: DataFrame does not contain the required columns: {columns}")
                return pd.DataFrame(columns=df.columns)

        # Print the number of rows in the original DataFrames
        print(f"Original number of rows in triple links: {df_triple_links.shape[0]}")
        print(f"Original number of rows in double links: {df_double_links.shape[0]}")

        # Filter the DataFrames
        filtered_triple_links = filter_dataframe(df_triple_links, ['Protein1', 'Protein2', 'Protein3'])
        filtered_double_links = filter_dataframe(df_double_links, ['Protein1', 'Protein2'])

        # Print the number of rows in the filtered DataFrames
        print(f"Filtered number of rows in triple links: {filtered_triple_links.shape[0]}")
        print(f"Filtered number of rows in double links: {filtered_double_links.shape[0]}")

        # Calculate and print the number of rows filtered out
        triple_links_filtered_out = df_triple_links.shape[0] - filtered_triple_links.shape[0]
        double_links_filtered_out = df_double_links.shape[0] - filtered_double_links.shape[0]
        print(f"Number of rows filtered out from triple links: {triple_links_filtered_out}")
        print(f"Number of rows filtered out from double links: {double_links_filtered_out}")

        # Ensure the directory exists
        output_dir = 'base_subcomplex_data'
        os.makedirs(output_dir, exist_ok=True)

        # Save filtered DataFrames to CSV files
        filtered_triple_links.to_csv(os.path.join(output_dir, 'filtered-triple-links.csv'), index=False)
        filtered_double_links.to_csv(os.path.join(output_dir, 'filtered-double-links.csv'), index=False)

        return filtered_triple_links, filtered_double_links
    #--------------------------------------------------------------------------


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
    transformer.find_and_remove_common_rows(final_file1, final_file2, 'final_combined_coupled_triples.csv', 'final_cleaned_doubles.csv')
    
    # Filter proteins
    #df_triple_links = pd.read_csv('base_subcomplex_data/filtered-triple-links.csv')
    #df_double_links = pd.read_csv('base_subcomplex_data/filtered-double-links.csv')
    df_triple_links = pd.read_csv('full-proteasome-data/combined_coupled_triples.csv')
    df_double_links = pd.read_csv('full-proteasome-data/final_cleaned_doubles.csv')
    transformer.filter_proteins(df_triple_links, df_double_links)