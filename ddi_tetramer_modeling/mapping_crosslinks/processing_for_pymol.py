import os, re, argparse
import pandas as pd

# Map (Protein, Copy) -> Chain
protein_copy_to_chain = {
    ('DDI1', 0): 'A', ('DDI1', 1): 'B',
    ('DDI2', 0): 'C', ('DDI2', 1): 'D',
}
ncopies = {'DDI1': 2, 'DDI2': 2}
valid_proteins = {p for p, _ in protein_copy_to_chain.keys()}

def ensure_dir(d): 
    if d and not os.path.exists(d): os.makedirs(d)

def to_int_res(x): 
    s = re.sub(r'\D+', '', str(x)); return int(s) if s else None

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

def process_bifunc_csv(input_csv, output_csv, remove_duplicates=False):
    if not input_csv or not os.path.exists(input_csv): return None
    df = pd.read_csv(input_csv, usecols=[0,1,2,3])  # only first 4 columns
    df.columns = ['Protein1','Residue1','Protein2','Residue2']
    df = df[(df.Protein1.isin(valid_proteins)) & (df.Protein2.isin(valid_proteins))].copy()
    df['Residue1'] = df['Residue1'].map(to_int_res)
    df['Residue2'] = df['Residue2'].map(to_int_res)
    df = df.dropna(subset=['Residue1','Residue2'])
    rows = []
    for _, row in df.iterrows():
        p1,p2 = row.Protein1, row.Protein2
        r1,r2 = int(row.Residue1), int(row.Residue2)
        for c1,c2 in copy_pairs_bifunc(p1,r1,p2,r2):
            ch1,ch2 = chains_for_pair(p1,c1,p2,c2)
            if ch1 and ch2: add_link(rows, ch1, r1, ch2, r2)
    out = pd.DataFrame(rows, columns=['Chain1','Residue1','Atom1','Chain2','Residue2','Atom2'])
    if remove_duplicates: out = dedupe_links(out)
    ensure_dir(os.path.dirname(output_csv))
    out.to_csv(output_csv, index=False)
    print(f"Wrote bifunctional: {output_csv} ({len(out)} links)")
    return out

def process_trifunc_csv(input_csv, output_csv, remove_duplicates=False):
    if not input_csv or not os.path.exists(input_csv): return None
    df = pd.read_csv(input_csv, header=0, names=['Protein1','Residue1','Protein2','Residue2','Protein3','Residue3'])
    df = df[(df.Protein1.isin(valid_proteins)) & (df.Protein2.isin(valid_proteins)) & (df.Protein3.isin(valid_proteins))].copy()
    for col in ['Residue1','Residue2','Residue3']: df[col] = df[col].map(to_int_res)
    df = df.dropna(subset=['Residue1','Residue2','Residue3'])
    rows = []
    for _, row in df.iterrows():
        p = [row.Protein1, row.Protein2, row.Protein3]
        r = [int(row.Residue1), int(row.Residue2), int(row.Residue3)]
        for c1,c2,c3 in copy_triples_trifunc(p, r):
            ch1 = protein_copy_to_chain.get((p[0], c1))
            ch2 = protein_copy_to_chain.get((p[1], c2))
            ch3 = protein_copy_to_chain.get((p[2], c3))
            if not (ch1 and ch2 and ch3): continue
            add_link(rows, ch1, r[0], ch2, r[1])  # 1-2
            add_link(rows, ch2, r[1], ch3, r[2])  # 2-3
            add_link(rows, ch3, r[2], ch1, r[0])  # 3-1
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
    ap.add_argument('--outdir', default='pymol_data', help='Output directory')
    ap.add_argument('--dedupe', action='store_true', help='Remove duplicate undirected links')
    args = ap.parse_args()

    ensure_dir(args.outdir)
    bi_out  = os.path.join(args.outdir, 'bifunc_pairs.csv')
    tri_out = os.path.join(args.outdir, 'trifunc_pairs.csv')

    process_bifunc_csv(args.bifunc, bi_out, remove_duplicates=args.dedupe)
    process_trifunc_csv(args.trifunc, tri_out, remove_duplicates=args.dedupe)

    print("Done.")

if __name__ == '__main__':
    main()