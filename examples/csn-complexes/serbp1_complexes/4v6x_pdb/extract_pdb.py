from Bio.PDB import MMCIFParser, PDBIO, Select

class MultiChainSelect(Select):
    def __init__(self, allowed_chains):
        """
        allowed_chains: A set (or list) of chain IDs that should be included in the output.
        """
        self.allowed_chains = set(allowed_chains)

    def accept_chain(self, chain):
        return chain.id in self.allowed_chains

def extract_multiple_chains_with_renaming(cif_file, chain_map, output_pdb):
    """
    Extracts multiple chains from a CIF file, renames them according to chain_map,
    and writes them all into a single PDB file.

    Args:
        cif_file (str): Path to the input .cif file.
        chain_map (dict): A mapping of {old_chain_id: new_chain_id} for renaming.
        output_pdb (str): Path to the output .pdb file.
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)

    # Rename each chain if it's in chain_map, keep the rest as is (or remove them if desired)
    for model in structure:
        for chain in model:
            if chain.id in chain_map:
                chain.id = chain_map[chain.id]

    # We'll select only the new chain IDs from chain_map
    new_ids = set(chain_map.values())

    # Save the structure with only the renamed chains
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=MultiChainSelect(new_ids))

    print(f"Renamed chains {list(chain_map.keys())} -> {list(chain_map.values())} and saved to {output_pdb}")

# Example usage
if __name__ == "__main__":
    cif_file = "4v6x.cif"
    
    # Suppose you have multiple chains you want extracted, mapping old IDs -> new IDs
    # e.g. chain_map = {"Ah": "A", "Bh": "B", "Ch": "C"}
    chain_map = {
        "Ah": "S",
        "AC": "A",
        "AD": "B",
        "Ae": "C",
        "Af": "D",
        "AK": "E",
        "AL": "F",
        "AM": "G",
        "CN": "H",
        "Cg": "I",
        "CG": "J",
        "CZ": "K",
        "AP": "L",
        "Ac": "M",
        "AB": "N",
        "B2": "U",
        "AJ": "O",
        "Aa": "P",
        "A5": "Q",
        "Az": "R"
    }
    
    output_pdb = "combined_chains_4v6x.pdb"
    
    extract_multiple_chains_with_renaming(cif_file, chain_map, output_pdb)
