#!/usr/bin/env python

import os, re, argparse
import pandas as pd

# Map (Protein, Copy) -> Chain
protein_copy_to_chain = {
    ('DDI1', 0): 'A', ('DDI1', 1): 'B',
    ('DDI2', 0): 'C', ('DDI2', 1): 'D',
}
# Reverse map: Chain -> Protein
chain_to_protein = {v: k[0] for k, v in protein_copy_to_chain.items()}

ncopies = {'DDI1': 2, 'DDI2': 0, 'DDI3': 0}
# the number of copies for the DDI proteins
valid_proteins = {p for p, _ in protein_copy_to_chain.keys()}

def ensure_dir(d): 
    if d and not os.path.exists(d): os.makedirs(d)

def to_int_res(x): 
    s = re.sub(r'\D+', '', str(x)); return int(s) if s else None

def load_valid_residues(pdb_file):
    """
    Parses a PDB file to find all valid (Protein, ResidueNum) pairs
    based on the ATOM records and chain_to_protein mapping.
    """
    if not pdb_file: return None
    valid_residues = set()
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    resnum_str = line[22:26].strip()
                    protein = chain_to_protein.get(chain_id)
                    
                    if protein and resnum_str:
                        try:
                            resnum = int(resnum_str)
                            valid_residues.add((protein, resnum))
                        except ValueError:
                            continue # Malformed resnum
        print(f"Loaded {len(valid_residues)} unique (Protein, Residue) pairs from {pdb_file}")
        if not valid_residues:
             print(f"Warning: No valid residues found. Check PDB chain IDs (e.g., A, B) match mappings.")
        return valid_residues
    except FileNotFoundError:
        print(f"Warning: PDB file not found at {pdb_file}. No residue validation will be performed.")
        return None

def copy_pairs_bifunc(p1, r1, p2, r2):
    # Rules: different proteins -> all copy pairs; same protein+same residue -> different copies only; same protein+diff residue -> intra and inter
    if p1 != p2:
        return [(i, j) for i in range(ncopies.get(p1,0)) for j in range(ncopies.get(p2,0))]
    if r1 == r2:
        return [(0,1), (1,0)]
    return [(0,0), (1,1), (0,1), (1,0)]

def copy_triples_trifunc(proteins, residues):
    # Pairwise options
    def opts(pA, rA, pB, rB):
        if pA != pB: return [(0,0),(0,1),(1,0),(1,1)]
        if rA == rB: return [(0,1),(1,0)]
        return [(0,0),(1,1),(0,1),(1,0)]
    c12, c23, c13 = opts(proteins[0],residues[0],proteins[1],residues[1]), \
                    opts(proteins[1],residues[1],proteins[2],residues[2]), \
                    opts(proteins[0],residues[0],proteins[2],residues[2])
    valid = set()
    for (a,b) in c12:
        for (b2,c) in c23:
            if b != b2: continue
            for (a2,c2) in c13:
                if a==a2 and c==c2: valid.add((a,b,c))
    return sorted(valid)

def chains_for_pair(p1,c1,p2,c2):
    ch1 = protein_copy_to_chain.get((p1,c1))
    ch2 = protein_copy_to_chain.get((p2,c2))
    return (ch1, ch2) if ch1 and ch2 else (None, None)

def add_link(rows, ch1, r1, ch2, r2, atom='CA'):
    rows.append([ch1, r1, atom, ch2, r2, atom])

def dedupe_links(df):
    key = df.apply(lambda r: tuple(sorted([(r.Chain1, r.Residue1, r.Atom1),
                                           (r.Chain2, r.Residue2, r.Atom2)])), axis=1)
    df2 = df.copy()
    df2['k'] = key
    df2 = df2.drop_duplicates('k').drop(columns='k')
    return df2

def process_bifunc_csv(input_csv, output_csv, remove_duplicates=False, valid_residues=None):
    if not input_csv or not os.path.exists(input_csv): return None
    
    # Load CSV and check if it's "naive" (missing copy columns)
    df = pd.read_csv(input_csv)
    is_naive = 'Copy1' not in df.columns or 'Copy2' not in df.columns
    
    if is_naive:
        print(f"Processing bifunctional (naive): {input_csv}. Will generate copy combinations.")
        # Re-read using only first 4 columns and force names
        df = pd.read_csv(input_csv, usecols=[0,1,2,3])
        df.columns = ['Protein1','Residue1','Protein2','Residue2']
    else:
        print(f"Processing bifunctional (copy-specific): {input_csv}. Will use provided Copy1/Copy2.")
        # Ensure correct types for copy columns if they exist
        df['Copy1'] = df['Copy1'].astype(int)
        df['Copy2'] = df['Copy2'].astype(int)
    
    # --- Standard Filtering ---
    df = df[(df.Protein1.isin(valid_proteins)) & (df.Protein2.isin(valid_proteins))].copy()
    df['Residue1'] = df['Residue1'].map(to_int_res)
    df['Residue2'] = df['Residue2'].map(to_int_res)
    df = df.dropna(subset=['Residue1','Residue2'])
    
    # --- New PDB Validation Step ---
    if valid_residues:
        initial_count = len(df)
        df = df[df.apply(lambda r: (r.Protein1, int(r.Residue1)) in valid_residues, axis=1)]
        df = df[df.apply(lambda r: (r.Protein2, int(r.Residue2)) in valid_residues, axis=1)]
        print(f"  PDB Validation: Retained {len(df)} of {initial_count} links after checking residues.")

    rows = []
    for _, row in df.iterrows():
        p1,p2 = row.Protein1, row.Protein2
        r1,r2 = int(row.Residue1), int(row.Residue2)
        
        if is_naive:
            # "Blow up" logic from original script
            for c1,c2 in copy_pairs_bifunc(p1,r1,p2,r2):
                ch1,ch2 = chains_for_pair(p1,c1,p2,c2)
                if ch1 and ch2: add_link(rows, ch1, r1, ch2, r2)
        else:
            # "Copy-specific" logic using columns from file
            c1, c2 = int(row.Copy1), int(row.Copy2)
            ch1, ch2 = chains_for_pair(p1, c1, p2, c2)
            if ch1 and ch2: add_link(rows, ch1, r1, ch2, r2)

    if not rows:
        print(f"Wrote bifunctional: {output_csv} (0 links)")
        return None

    out = pd.DataFrame(rows, columns=['Chain1','Residue1','Atom1','Chain2','Residue2','Atom2'])
    if remove_duplicates: out = dedupe_links(out)
    ensure_dir(os.path.dirname(output_csv))
    out.to_csv(output_csv, index=False)
    print(f"Wrote bifunctional: {output_csv} ({len(out)} links)")
    return out

def process_trifunc_csv(input_csv, output_csv, remove_duplicates=False, valid_residues=None):
    if not input_csv or not os.path.exists(input_csv): return None

    # Load CSV and check if it's "naive"
    df = pd.read_csv(input_csv)
    is_naive = 'Copy1' not in df.columns or 'Copy2' not in df.columns or 'Copy3' not in df.columns
    
    if is_naive:
        print(f"Processing trifunctional (naive): {input_csv}. Will generate copy combinations.")
        # Re-read using only first 6 columns and force names
        df = pd.read_csv(input_csv, usecols=range(6))
        df.columns = ['Protein1','Residue1','Protein2','Residue2','Protein3','Residue3']
    else:
        print(f"Processing trifunctional (copy-specific): {input_csv}. Will use provided Copy1/2/3.")
        df['Copy1'] = df['Copy1'].astype(int)
        df['Copy2'] = df['Copy2'].astype(int)
        df['Copy3'] = df['Copy3'].astype(int)
        
    # --- Standard Filtering ---
    df = df[(df.Protein1.isin(valid_proteins)) & (df.Protein2.isin(valid_proteins)) & (df.Protein3.isin(valid_proteins))].copy()
    for col in ['Residue1','Residue2','Residue3']: df[col] = df[col].map(to_int_res)
    df = df.dropna(subset=['Residue1','Residue2','Residue3'])

    # --- New PDB Validation Step ---
    if valid_residues:
        initial_count = len(df)
        df = df[df.apply(lambda r: (r.Protein1, int(r.Residue1)) in valid_residues, axis=1)]
        df = df[df.apply(lambda r: (r.Protein2, int(r.Residue2)) in valid_residues, axis=1)]
        df = df[df.apply(lambda r: (r.Protein3, int(r.Residue3)) in valid_residues, axis=1)]
        print(f"  PDB Validation: Retained {len(df)} of {initial_count} links after checking residues.")
        
    rows = []
    for _, row in df.iterrows():
        p = [row.Protein1, row.Protein2, row.Protein3]
        r = [int(row.Residue1), int(row.Residue2), int(row.Residue3)]
        
        if is_naive:
            # "Blow up" logic from original script
            for c1,c2,c3 in copy_triples_trifunc(p, r):
                ch1 = protein_copy_to_chain.get((p[0], c1))
                ch2 = protein_copy_to_chain.get((p[1], c2))
                ch3 = protein_copy_to_chain.get((p[2], c3))
                if not (ch1 and ch2 and ch3): continue
                add_link(rows, ch1, r[0], ch2, r[1])  # 1-2
                add_link(rows, ch2, r[1], ch3, r[2])  # 2-3
                add_link(rows, ch3, r[2], ch1, r[0])  # 3-1
        else:
            # "Copy-specific" logic using columns from file
            c = [int(row.Copy1), int(row.Copy2), int(row.Copy3)]
            ch1 = protein_copy_to_chain.get((p[0], c[0]))
            ch2 = protein_copy_to_chain.get((p[1], c[1]))
            ch3 = protein_copy_to_chain.get((p[2], c[2]))
            if not (ch1 and ch2 and ch3): continue
            add_link(rows, ch1, r[0], ch2, r[1])  # 1-2
            add_link(rows, ch2, r[1], ch3, r[2])  # 2-3
            add_link(rows, ch3, r[2], ch1, r[0])  # 3-1

    if not rows:
        print(f"Wrote trifunctional: {output_csv} (0 links)")
        return None

    out = pd.DataFrame(rows, columns=['Chain1','Residue1','Atom1','Chain2','Residue2','Atom2'])
    if remove_duplicates: out = dedupe_links(out)
    ensure_dir(os.path.dirname(output_csv))
    out.to_csv(output_csv, index=False)
    print(f"Wrote trifunctional: {output_csv} ({len(out)} links)")
    return out

def main():
    ap = argparse.ArgumentParser(description="Generate pairwise links from bi/tri-functional XLs (ChimeraX/PyMOL format).")
    ap.add_argument('--bifunc', default='input_data/reduced_ddi_bifunctional.csv', help='Bifunctional CSV (usecols: Protein1,Residue1,Protein2,Residue2)')
    ap.add_argument('--trifunc', default='input_data/reduced_ddi_trifunctional.csv', help='Trifunctional CSV (Protein1,Residue1,Protein2,Residue2,Protein3,Residue3)')
    
    # --- NEW ARGUMENT ---
    ap.add_argument('--pdb', default=None, help='(Optional) PDB file to validate residues against. If provided, links with missing residues will be skipped.')
    
    ap.add_argument('--outdir', default='pymol_data', help='Output directory')
    ap.add_argument('--dedupe', action='store_true', help='Remove duplicate undirected links')
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # --- Load valid residues from PDB (if provided) ---
    valid_residues = load_valid_residues(args.pdb)
    
    bi_out  = os.path.join(args.outdir, 'bifunc_pairs.csv')
    tri_out = os.path.join(args.outdir, 'trifunc_pairs.csv')

    # --- Pass the valid_residues set to the processing functions ---
    process_bifunc_csv(args.bifunc, bi_out, remove_duplicates=args.dedupe, valid_residues=valid_residues)
    process_trifunc_csv(args.trifunc, tri_out, remove_duplicates=args.dedupe, valid_residues=valid_residues)

    print("Done.")

if __name__ == '__main__':
    main()