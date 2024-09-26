import pandas as pd

def find_common_rows(file1, file2):
    """Find and print common rows between two datasets based on Residue1 and Residue2 columns."""
    # Read the CSV files into DataFrames
    df1 = pd.read_csv(file1, delimiter=',')
    df2 = pd.read_csv(file2, delimiter=',')

    # Print column names for debugging
    print("Columns in double XLs:", df1.columns)
    print("Columns in triple XLs:", df2.columns)

    # Extract Residue1 and Residue2 columns and convert each row to a set of values
    df1_sets = df1[['Residue1', 'Residue2']].apply(lambda row: frozenset(row), axis=1)
    df2_sets = df2[['Residue1', 'Residue2']].apply(lambda row: frozenset(row), axis=1)

    # Identify common rows based on set equality
    common_rows_indices = df1_sets.isin(df2_sets)
    common_rows = df1[common_rows_indices]

    # Count common rows in df2 that are also in df1
    common_count = df2_sets.isin(df1_sets).sum()

    # Calculate total number of rows in df2
    total_rows_df2 = len(df2)

    # Calculate percentage of overlap
    percentage_overlap = (common_count / total_rows_df2) * 100

    if not common_rows.empty:
        print("Common rows found:")
        print(common_rows.to_markdown(index=False))
    else:
        print("No common rows found.")
    
    print(f"Number of rows in triple XLs that are also present in double XLs: {common_count}")
    print(f"Total number of rows in triple XLs: {total_rows_df2}")
    print(f"Percentage of overlap with double XLs: {percentage_overlap:.2f}%")

    # Remove the rows in df1 that are also in df2
    df1_filtered = df1[~common_rows_indices]

    # Calculate the number of rows in the resulting df1 and the number of rows that were removed
    rows_removed = len(df1) - len(df1_filtered)
    rows_remaining = len(df1_filtered)

    # Print the resulting df1
#    print("Resulting double XLs after removing common rows:")
#    print(df1_filtered.to_markdown(index=False))

    # Print the counts
    print(f"Number of rows removed from double XLs: {rows_removed}")
    print(f"Number of rows remaining in double XLs: {rows_remaining}")

    # Write the resulting df1 to a CSV file
    output_file = 'base_subcomplex_data/filtered_double_xls.csv'
    df1_filtered.to_csv(output_file, index=False)
#    print(f"Filtered double XLs saved to {output_file}")

    # Print the rows in triple XLs that are also present in double XLs
    common_rows_in_df2 = df2[df2_sets.isin(df1_sets)]
#    print("Rows in triple XLs that are also present in double XLs:")
#    print(common_rows_in_df2.to_markdown(index=False))

    # Print the rows in triple XLs that are not present in double XLs
    non_common_rows_in_df2 = df2[~df2_sets.isin(df1_sets)]
#    print("Rows in triple XLs that are not present in double XLs:")
#    print(non_common_rows_in_df2.to_markdown(index=False))

# Example usage
if __name__ == "__main__":
    file1 = 'base_subcomplex_data/double.csv'
    file2 = 'base_subcomplex_data/triple.csv'

    find_common_rows(file1, file2)
