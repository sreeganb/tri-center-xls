from Bio.PDB import MMCIFParser, PDBIO, Select

class MultiChainSelect(Select):
    def __init__(self, allowed_chains):
        self.allowed_chains = set(allowed_chains)

    def accept_chain(self, chain):
        return chain.id in self.allowed_chains

def process_chunk(cif_file, chain_group, output_name):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_file)
    
    # Rename chains
    for model in structure:
        for chain in model:
            if chain.id in chain_group["map"]:
                chain.id = chain_group["map"][chain.id]
    
    # Reset atom serial numbers
    current_serial = 1
    for model in structure:
        for chain in model:
            if chain.id in chain_group["new_ids"]:
                for residue in chain:
                    for atom in residue:
                        atom.serial_number = current_serial
                        current_serial += 1

    # Save chunk
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_name, select=MultiChainSelect(chain_group["new_ids"]))
    print(f"Saved {output_name} with {current_serial-1} atoms")

# Split your chain_map into groups
chain_groups = [
    {
        "map": {
            "Ah": "S", "AC": "A", "AD": "B", "Ae": "C", "Af": "D",
            "AK": "E", "AL": "F", "AM": "G", "CN": "H", "Cg": "I",
            "CG": "J", "CZ": "K", "AP": "L", "Ac": "M", "AB": "N",
            "B2": "U", "AJ": "O", "Aa": "P", "Az": "R"

        },
        "new_ids": {"S", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "U"},
        "output": "group1_4v6x.pdb"
    },
    {
        "map": {
            "A5": "Q"
        },
        "new_ids": { "Q"},
        "output": "group2_4v6x.pdb"
    }
]

for group in chain_groups:
    process_chunk("4v6x.cif", group, group["output"])