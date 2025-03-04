import os
import pandas as pd

class CrosslinkProcessor:
    def __init__(self, input_file):
        self.input_file = input_file
        self.data = None
        self.processed_data = None
        self.protein_c = None
        self.no_protein_c = None
        self.intra_crosslinks = None
        self.inter_crosslinks = None
        self.partial_crosslinks = None
        self.output_dir = 'output_data'
        os.makedirs(self.output_dir, exist_ok=True)

    def read_data(self):
        # Read the CSV file and select specified columns
        self.data = pd.read_csv(self.input_file)
        self.processed_data = self.data[['XL A', 'Protein A', 'XL B', 'Protein B', 'XL C', 'Protein C']].fillna('')

    def process_protein_c(self):
        # Split dataframe based on 'Protein C' being empty or not
        self.no_protein_c = self.processed_data[self.processed_data['Protein C'] == ''].reset_index(drop=True)
        self.protein_c = self.processed_data[self.processed_data['Protein C'] != ''].reset_index(drop=True)
        # Drop duplicates
        self.protein_c.drop_duplicates(inplace=True)

    def segregate_crosslinks(self):
        # Intra protein crosslinks: all proteins are the same
        intra_mask = (
            (self.protein_c['Protein A'] == self.protein_c['Protein B']) &
            (self.protein_c['Protein B'] == self.protein_c['Protein C'])
        )
        self.intra_crosslinks = self.protein_c[intra_mask].reset_index(drop=True)

        # Inter protein crosslinks: all proteins are different
        inter_mask = (
            (self.protein_c['Protein A'] != self.protein_c['Protein B']) &
            (self.protein_c['Protein A'] != self.protein_c['Protein C']) &
            (self.protein_c['Protein B'] != self.protein_c['Protein C'])
        )
        self.inter_crosslinks = self.protein_c[inter_mask].reset_index(drop=True)

        # Partial crosslinks: two proteins are the same, one is different
        partial_mask = ~intra_mask & ~inter_mask
        self.partial_crosslinks = self.protein_c[partial_mask].reset_index(drop=True)

    def write_output(self):
        # Desired column order
        desired_order = ['Protein A', 'XL A', 'Protein B', 'XL B', 'Protein C', 'XL C']

        # Reorder columns in protein_c dataframe
        self.protein_c = self.protein_c[desired_order]

        # Reorder columns in no_protein_c dataframe
        self.no_protein_c = self.no_protein_c[['Protein A', 'XL A', 'Protein B', 'XL B']]

        # Reorder columns in intra, inter, and partial crosslinks dataframes
        self.intra_crosslinks = self.intra_crosslinks[desired_order]
        self.inter_crosslinks = self.inter_crosslinks[desired_order]
        self.partial_crosslinks = self.partial_crosslinks[desired_order]

        # Save dataframes to CSV files
        self.protein_c.to_csv(os.path.join(self.output_dir, 'protein_c.csv'), index=False)
        self.no_protein_c.to_csv(os.path.join(self.output_dir, 'no_protein_c.csv'), index=False)
        self.intra_crosslinks.to_csv(os.path.join(self.output_dir, 'intra_protein_crosslinks.csv'), index=False)
        self.inter_crosslinks.to_csv(os.path.join(self.output_dir, 'inter_protein_crosslinks.csv'), index=False)
        self.partial_crosslinks.to_csv(os.path.join(self.output_dir, 'partial_crosslinks.csv'), index=False)
    
        # Combine 'Protein A', 'Protein B', and 'Protein C' into a single series
        all_proteins = pd.concat([
            self.protein_c['Protein A'],
            self.protein_c['Protein B'],
            self.protein_c['Protein C']
        ], ignore_index=True)

        # Remove empty strings if any
        all_proteins = all_proteins[all_proteins != '']

        # Compute frequency counts
        protein_counts = all_proteins.value_counts()

        # Save frequency counts to CSV
        protein_counts.to_csv(os.path.join(self.output_dir, 'protein_counts.csv'), header=['Count'])
    
    def process_all(self):
        self.read_data()
        self.process_protein_c()
        self.segregate_crosslinks()
        self.write_output()

if __name__ == '__main__':
    processor = CrosslinkProcessor('input_data/tsto-supplementary.csv')
    processor.process_all()