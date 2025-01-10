from Bio.PDB import MMCIFParser, PDBIO, Select

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

def extract_chain_with_renaming(cif_file, chain_id, output_pdb, new_chain_id):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)

    # Rename the chain ID if it matches
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                chain.id = new_chain_id

    # Save the structure with the renamed chain
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=ChainSelect(new_chain_id))

# Input CIF file, desired chain ID, and output PDB file
#cif_file = "8k2c.cif"  # Replace with your CIF file
cif_file = "4v6x.cif"  # Replace with your CIF file
#chain_id = "CB"             # Replace with your desired chain ID
chain_id = "Ah"             # Replace with your desired chain ID
#new_chain_id = "S"              # New single-character chain ID
new_chain_id = "A"              # New single-character chain ID
#output_pdb = "chain_CB.pdb"
output_pdb = "chain_B_4v6x.pdb"


# Run the extraction with renaming
extract_chain_with_renaming(cif_file, chain_id, output_pdb, new_chain_id)

print(f"Chain {chain_id} renamed to {new_chain_id} and saved to {output_pdb}")
