#!/usr/bin/env python3
"""
Compute per-frame crosslink satisfaction from an RMF3 using IMP/PMI.
Outputs:
  - distances_by_frame.csv: frame-indexed table of per-XL min distances
  - satisfied_counts.csv: per-frame count of XLs with distance <= cutoff
  - time series plot of satisfied/violated XLs per frame (PDF)

Crosslink CSV formats supported:
  A) PMI-style:
     id, Prot A, Residue A, Copy A, Prot B, Residue B, Copy B
     (Copy columns optional)
  B) Chain-style:
     Chain1, Residue1, Chain2, Residue2
     Requires CHAIN_MAP below to map chain -> (molecule_name, copy_index)
"""

import os
import sys
import argparse
from collections import defaultdict
import itertools

import RMF
import IMP
import IMP.core
import IMP.atom
import IMP.rmf

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Customize as needed for Chain-style CSVs (B)
CHAIN_MAP = {
    # chain_id: (molecule_name_in_PMI, copy_index_int)
    'A': ('DDI1', 0),
    'B': ('DDI1', 1),
    'C': ('DDI2', 0),
    'D': ('DDI2', 1),
}

def _as_int_or_none(x):
    if x is None or x == "" or (isinstance(x, float) and np.isnan(x)):
        return None
    return int(x)

def _select_particles(hier, mol, resi, copy_idx=None):
    copy_idx = _as_int_or_none(copy_idx)
    sel = IMP.atom.Selection(
        hier,
        molecule=mol,
        residue_index=int(resi),
        copy_index=copy_idx,
    ).get_selected_particles()
    # Keep only particles with coordinates
    return [p for p in sel if IMP.core.XYZ.get_is_setup(p)]

def _read_crosslinks(csv_path):
    """
    Returns:
      xls: list of dict entries with keys:
           id, ProtA, ResA, CopyA, ProtB, ResB, CopyB
    """
    df = pd.read_csv(csv_path)
    cols = {c.lower(): c for c in df.columns}

    entries = []
    if {'prot a', 'residue a', 'prot b', 'residue b'}.issubset(cols.keys()):
        # PMI-style
        for _, row in df.iterrows():
            entries.append(dict(
                id=str(row.get(cols.get('id'), f'xl_{_}')),
                ProtA=row[cols['prot a']],
                ResA=int(row[cols['residue a']]),
                CopyA=row.get(cols.get('copy a'), None),
                ProtB=row[cols['prot b']],
                ResB=int(row[cols['residue b']]),
                CopyB=row.get(cols.get('copy b'), None),
            ))
    elif {'chain1', 'residue1', 'chain2', 'residue2'}.issubset(cols.keys()):
        # Chain-style, map to PMI molecules/copies
        for i, row in df.iterrows():
            c1 = str(row[cols['chain1']]).strip()
            c2 = str(row[cols['chain2']]).strip()
            if c1 not in CHAIN_MAP or c2 not in CHAIN_MAP:
                # skip if we cannot map; user should adjust CHAIN_MAP
                continue
            m1, k1 = CHAIN_MAP[c1]
            m2, k2 = CHAIN_MAP[c2]
            entries.append(dict(
                id=str(row.get(cols.get('id'), f'xl_{i}')),
                ProtA=m1, ResA=int(row[cols['residue1']]), CopyA=k1,
                ProtB=m2, ResB=int(row[cols['residue2']]), CopyB=k2,
            ))
    else:
        raise ValueError("CSV must be PMI-style (Prot/Residue/Copy) or Chain-style (Chain/Residue).")

    if len(entries) == 0:
        raise ValueError("No crosslinks parsed. Check CSV columns or CHAIN_MAP.")
    return entries

def _build_pairs_for_hierarchy(hier, xls_entries):
    """
    For each XL entry, build the list of particle pairs present in this hierarchy.
    Returns:
      dict: id -> list of (particle, particle)
    """
    id_to_pairs = {}
    for e in xls_entries:
        s1 = _select_particles(hier, e['ProtA'], e['ResA'], e.get('CopyA', None))
        s2 = _select_particles(hier, e['ProtB'], e['ResB'], e.get('CopyB', None))
        pairs = [(p0, p1) for p0, p1 in itertools.product(s1, s2) if p0 != p1]
        if pairs:
            id_to_pairs.setdefault(e['id'], []).extend(pairs)
    return id_to_pairs

def _min_distance_for_pairs(pairs):
    dmin = np.inf
    for p0, p1 in pairs:
        d = IMP.core.get_distance(IMP.core.XYZ(p0), IMP.core.XYZ(p1))
        if d < dmin:
            dmin = d
    return dmin if np.isfinite(dmin) else np.nan

def compute_distances_per_frame(rmf_path, xls_entries, step=1):
    """
    Returns:
      distances_df: DataFrame shape (n_frames, n_xls) with min distance per XL per frame
      frame_indices: list of frame indices used
    """
    mdl = IMP.Model()
    rh = RMF.open_rmf_file_read_only(rmf_path)
    hier = IMP.rmf.create_hierarchies(rh, mdl)[0]
    n_frames = rh.get_number_of_frames()
    if n_frames == 0:
        raise RuntimeError("RMF has zero frames.")
    ids = [e['id'] for e in xls_entries]

    # Pre-select particle pairs per XL for this hierarchy
    id_to_pairs = _build_pairs_for_hierarchy(hier, xls_entries)
    # Keep only XLs that exist in the representation
    valid_ids = [xlid for xlid in ids if xlid in id_to_pairs]
    if not valid_ids:
        raise RuntimeError("No crosslinks match particles in this RMF hierarchy.")

    frame_indices = list(range(0, n_frames, step))
    data = {xlid: [] for xlid in valid_ids}

    for fr in frame_indices:
        IMP.rmf.load_frame(rh, RMF.FrameID(fr))
        for xlid in valid_ids:
            dmin = _min_distance_for_pairs(id_to_pairs[xlid])
            data[xlid].append(dmin)

    distances_df = pd.DataFrame(data, index=frame_indices)
    distances_df.index.name = 'frame'
    return distances_df, frame_indices

def plot_time_series(distances_df, cutoff, outdir, prefix="xl"):
    """
    Plot time series of satisfied vs violated crosslinks per frame.
    """
    os.makedirs(outdir, exist_ok=True)
    
    # Count satisfied and violated per frame
    satisfied = (distances_df <= float(cutoff)).sum(axis=1)
    violated = (distances_df > float(cutoff)).sum(axis=1)
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(satisfied.index, satisfied.values, 
             label=f'Satisfied (≤ {cutoff} Å)', 
             color='green', linewidth=2, marker='o', markersize=3)
    plt.plot(violated.index, violated.values, 
             label=f'Violated (> {cutoff} Å)', 
             color='red', linewidth=2, marker='s', markersize=3)
    
    plt.xlabel('Frame Number', fontsize=12)
    plt.ylabel('Number of Crosslinks', fontsize=12)
    plt.title('Crosslink Satisfaction Over Trajectory', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    pdf_path = os.path.join(outdir, f"{prefix}_satisfaction_timeseries.pdf")
    plt.savefig(pdf_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    return pdf_path

def summarize_and_save(distances_df, cutoff, outdir, prefix="xl"):
    os.makedirs(outdir, exist_ok=True)
    
    # Save distances CSV
    dist_csv = os.path.join(outdir, f"{prefix}_distances_by_frame.csv")
    distances_df.to_csv(dist_csv)

    # Per-frame satisfied counts
    satisfied = (distances_df <= float(cutoff)).sum(axis=1)
    counts_csv = os.path.join(outdir, f"{prefix}_satisfied_counts.csv")
    satisfied.to_csv(counts_csv, header=['num_satisfied'])

    # Generate time series plot
    pdf_path = plot_time_series(distances_df, cutoff, outdir, prefix)

    return dist_csv, counts_csv, pdf_path

def main():
    ap = argparse.ArgumentParser(description="Per-frame crosslink satisfaction from RMF3.")
    ap.add_argument('--rmf', required=True, help='Path to RMF3 file')
    ap.add_argument('--xls', required=True, help='Crosslink CSV (PMI or Chain style)')
    ap.add_argument('--cutoff', type=float, default=30.0, help='Satisfaction cutoff in Å (default: 30)')
    ap.add_argument('--step', type=int, default=1, help='Frame step (default: 1)')
    ap.add_argument('--outdir', default='analys', help='Output directory (default: analys)')
    ap.add_argument('--prefix', default='xl', help='Output filename prefix (default: xl)')
    args = ap.parse_args()

    xls_entries = _read_crosslinks(args.xls)
    distances_df, frame_indices = compute_distances_per_frame(args.rmf, xls_entries, step=args.step)
    dist_csv, counts_csv, pdf_path = summarize_and_save(distances_df, args.cutoff, args.outdir, prefix=args.prefix)
    
    print(f"Wrote:\n  {dist_csv}\n  {counts_csv}\n  {pdf_path}")
    print(f"Frames analyzed: {len(frame_indices)} (first={frame_indices[0]}, last={frame_indices[-1]})")
    print(f"Satisfied counts head:\n{(distances_df <= args.cutoff).sum(axis=1).head()}")

if __name__ == '__main__':
    main()