import IMP
import IMP.atom
import IMP.rmf
import RMF
import os
import numpy as np

# Map (molecule_name, copy_index) -> desired chain ID
MOLECULE_COPY_TO_CHAIN = {
    ('DDI1', 0): 'A',
    ('DDI1', 1): 'B',
    ('DDI2', 0): 'C',
    ('DDI2', 1): 'D',
}

# FASTA files for each molecule (residue index 1-based)
FASTA_FILES = {
    'DDI1': './input_data/ddi1.fasta',
    'DDI2': './input_data/ddi2.fasta',
}

# Standard 1-letter to 3-letter amino acid code
AA_1TO3 = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
    'X': 'UNK', '-': 'UNK',
}

def _read_fasta(fasta_path):
    """Return dict {residue_index (1-based) -> resname3}."""
    seq_dict = {}
    if not os.path.exists(fasta_path):
        return seq_dict
    with open(fasta_path) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                continue
            seq += line.strip()
    for i, aa in enumerate(seq, start=1):
        seq_dict[i] = AA_1TO3.get(aa.upper(), 'UNK')
    return seq_dict

def _get_chain_parent(h):
    """Climb parents to find a Chain node; return Chain or None."""
    p = h.get_parent()
    while p and p != IMP.atom.Hierarchy():
        try:
            return IMP.atom.Chain(p)
        except Exception:
            pass
        p = p.get_parent()
    return None

def _try_set_chain_id(hier):
    """Assign chain IDs to Chain nodes using MOLECULE_COPY_TO_CHAIN where missing."""
    # Iterate molecule nodes
    mol_nodes = IMP.atom.get_by_type(hier, IMP.atom.MOLECULE_TYPE)
    for mol_h in mol_nodes:
        mol = IMP.atom.Molecule(mol_h)
        mol_name = mol.get_name()
        # Copy index is stored on hierarchy as a "copy index" attribute in PMI
        copy_idx = None
        try:
            copy_idx = IMP.atom.Copy(mol_h).get_copy_index()
        except Exception:
            # Fallback: many PMI builds store as "copy_index" attribute
            if mol_h.get_has_attribute("copy_index"):
                copy_idx = mol_h.get_attribute("copy_index")
        # Default to 0 if not present
        if copy_idx is None:
            copy_idx = 0
        chain_id = MOLECULE_COPY_TO_CHAIN.get((mol_name, int(copy_idx)))
        if not chain_id:
            continue
        # Assign to existing Chain nodes under this molecule
        chain_nodes = [c for c in IMP.atom.get_by_type(mol_h, IMP.atom.CHAIN_TYPE)]
        if chain_nodes:
            for ch in chain_nodes:
                try:
                    IMP.atom.Chain(ch).set_id(chain_id)
                except Exception:
                    pass
        else:
            # Create a Chain node if none exists, and reparent residues under it
            ch = IMP.atom.Hierarchy.setup_particle(IMP.Particle(mol_h.get_model()), name="chain_"+chain_id)
            ch.set_parent(mol_h)
            try:
                IMP.atom.Chain.setup_particle(ch.get_particle()).set_id(chain_id)
            except Exception:
                pass
            # Move residues under this chain
            for res in list(IMP.atom.get_by_type(mol_h, IMP.atom.RESIDUE_TYPE)):
                # Only move direct children without a chain parent
                if _get_chain_parent(res) is None:
                    IMP.atom.Hierarchy(res).set_parent(ch)

def _set_residue_types_from_fasta(hier, fasta_map):
    """Fill residue types using FASTA sequences by molecule name and residue index."""
    # fasta_map: {molecule_name: {resid (1-based): resname3}}
    mol_nodes = IMP.atom.get_by_type(hier, IMP.atom.MOLECULE_TYPE)
    for mol_h in mol_nodes:
        mol = IMP.atom.Molecule(mol_h)
        mol_name = mol.get_name()
        if mol_name not in fasta_map:
            continue
        seq_dict = fasta_map[mol_name]
        res_nodes = IMP.atom.get_by_type(mol_h, IMP.atom.RESIDUE_TYPE)
        for res_h in res_nodes:
            res = IMP.atom.Residue(res_h)
            resid = res.get_index()
            # Skip if type already set to a known standard residue
            try:
                rts = res.get_residue_type().get_string()
                if rts != "UNK":
                    continue
            except Exception:
                pass
            resname = seq_dict.get(resid)
            if not resname:
                continue
            # Set residue type
            try:
                rt = IMP.atom.get_residue_type(resname)
                res.set_residue_type(rt)
            except Exception:
                pass

def rmf_to_pdbs(rmf_in, output_dir=None, frames=None, sel_state=0):
    """
    Convert RMF frames to PDB files with fixed residue names and chain IDs.
    """
    # Setup output directory
    if output_dir is None:
        rmf_name = os.path.splitext(os.path.basename(rmf_in))[0]
        output_dir = os.path.join("pdbs", rmf_name)
    os.makedirs(output_dir, exist_ok=True)

    # Load RMF file
    m = IMP.Model()
    f = RMF.open_rmf_file_read_only(rmf_in)
    n_frames = f.get_number_of_frames()
    print(f"Total frames in RMF: {n_frames}")

    # Get hierarchy
    root = IMP.rmf.create_hierarchies(f, m)[0]

    # Select specific state if needed
    states = IMP.atom.get_by_type(root, IMP.atom.STATE_TYPE)
    hier = states[sel_state] if states else root

    # Load FASTA sequences
    fasta_map = {}
    for mol_name, fasta_path in FASTA_FILES.items():
        fasta_map[mol_name] = _read_fasta(fasta_path)
        print(f"Loaded FASTA for {mol_name}: {len(fasta_map[mol_name])} residues")

    # Assign chains and residue names once
    _try_set_chain_id(hier)
    _set_residue_types_from_fasta(hier, fasta_map)

    # Determine frames to export
    if frames is None:
        frames = range(n_frames)
    print(f"Exporting {len(frames)} frames to {output_dir}")

    # Export each frame
    for i in frames:
        if i % 100 == 0:
            print(f"Processing frame {i}/{n_frames}")
        IMP.rmf.load_frame(f, RMF.FrameID(i))
        pdb_path = os.path.join(output_dir, f"frame_{i}.pdb")
        IMP.atom.write_pdb(hier, pdb_path)

    print(f"Done! Exported {len(frames)} PDB files to {output_dir}")
    del f, m

# Example usage
if __name__ == "__main__":
    # Export every 100th frame with labeling
    rmf_to_pdbs("./input_data/A_models_cluster2_0_aligned.rmf3", frames=range(0, 10000, 100))