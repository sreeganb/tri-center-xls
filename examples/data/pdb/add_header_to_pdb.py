#!/usr/bin/env python3.11

def add_header_to_pdb(pdb_file, proteins, chain_ids, output_file):
  """
  Adds a header to a PDB file with protein names for each chain ID and creates a new PDB file.

  Args:
    pdb_file: Path to the input PDB file.
    proteins: A list of protein names.
    chain_ids: A list of chain IDs.
    output_file: Path to the output PDB file.
  """

  if len(proteins) != len(chain_ids):
    raise ValueError("Number of proteins and chain IDs must be equal.")

  chain_info = dict(zip(chain_ids, proteins))

  with open(pdb_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for chain_id, protein_name in chain_info.items():
      f_out.write(f"HEADER    {protein_name} CHAIN {chain_id}\n")
    for line in f_in:
      f_out.write(line)

# Example usage:
proteins = ['Rpt1', 'Rpt2', 'Rpt3', 'Rpt4', 'Rpt5', 'Rpt6', 'Rpn2']
chain_ids = ['v', 'w', 'y', 'z', '0', 'x', '1']
input_pdb = 'base_proteasome.pdb'
output_pdb = 'modified_base_proteasome.pdb'
add_header_to_pdb(input_pdb, proteins, chain_ids, output_pdb)
