# Load necessary ChimeraX modules
from chimerax.core.commands import run
from chimerax.core.errors import UserError
from chimerax.atomic import Residue, Atoms, AtomicStructure
from chimerax.geometry import distance
import statistics

# Define the file paths (hardcoded)
PDB_FILE = "./pdbs/A_models_cluster2_0_aligned/frame_0.pdb"
#CROSSLINK_FILE = "pymol_data/bifunc_pairs.csv"
CROSSLINK_FILE = "pymol_data/trifunc_pairs.csv"
OUTPUT_STATS_FILE = "analys/crosslink_statistics.txt"

# Function to load PDB structure
def load_structure(session, pdb_file):
    run(session, f'open {pdb_file}')

# Function to read crosslink data
def read_crosslink_data(file_path):
    crosslinks = []
    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            if line.strip():
                chain1, res1, atom1, chain2, res2, atom2 = line.strip().split(',')
                crosslinks.append((chain1.strip(), int(res1.strip()), atom1.strip(), 
                                   chain2.strip(), int(res2.strip()), atom2.strip()))
    return crosslinks

# Function to get an atom by chain, residue number, and atom name
def find_atom_by_chain_residue(model, chain, residue_num, atom_name):
    for residue in model.residues:
        if residue.chain_id == chain and residue.number == residue_num:
            for atom in residue.atoms:
                if atom.name == atom_name:
                    return atom
    return None

# Function to visualize crosslinks on the structure
def plot_crosslinks(session, pdb_model, crosslinks, cutoff=30.0):
    violation_count = 0
    satisfied_count = 0
    processed_crosslinks = set()
    distances = []
    missing_count = 0
    
    for crosslink in crosslinks:
        chain1, res1, atom1, chain2, res2, atom2 = crosslink

        # Skip duplicate crosslinks
        if (chain1, res1, atom1, chain2, res2, atom2) in processed_crosslinks or \
           (chain2, res2, atom2, chain1, res1, atom1) in processed_crosslinks:
            continue
        
        processed_crosslinks.add((chain1, res1, atom1, chain2, res2, atom2))

        # Find atoms in the structure by chain, residue number, and atom name
        atom1_obj = find_atom_by_chain_residue(pdb_model, chain1, res1, atom1)
        atom2_obj = find_atom_by_chain_residue(pdb_model, chain2, res2, atom2)
        
        # If both atoms are found, create or modify a link (distance line)
        if atom1_obj and atom2_obj:
            dist = distance(atom1_obj.coord, atom2_obj.coord)
            distances.append(dist)
            
            if dist < cutoff:
                color = 'dark green'
                satisfied_count += 1
            else:
                color = 'dark red'
                violation_count += 1
            
            try:
                run(session, f'distance {atom1_obj.atomspec} {atom2_obj.atomspec} color {color} radius 0.6')
            except UserError:
                run(session, f'distance style {atom1_obj.atomspec} {atom2_obj.atomspec} color {color} radius 0.6')
        else:
            missing_count += 1
    
    # Calculate statistics
    stats = {
        'total_read': len(crosslinks),
        'total_unique': len(processed_crosslinks),
        'satisfied': satisfied_count,
        'violated': violation_count,
        'missing': missing_count,
        'cutoff': cutoff,
    }
    
    if distances:
        stats['min_distance'] = min(distances)
        stats['max_distance'] = max(distances)
        stats['mean_distance'] = statistics.mean(distances)
        stats['median_distance'] = statistics.median(distances)
        stats['stdev_distance'] = statistics.stdev(distances) if len(distances) > 1 else 0.0
    else:
        stats['min_distance'] = None
        stats['max_distance'] = None
        stats['mean_distance'] = None
        stats['median_distance'] = None
        stats['stdev_distance'] = None
    
    return stats, distances

def print_statistics(stats):
    """Print statistics to console."""
    print("\n" + "=" * 60)
    print("CROSSLINK STATISTICS")
    print("=" * 60)
    print(f"\nDistance Cutoff: {stats['cutoff']:.1f} Å\n")
    print(f"Total crosslinks read:        {stats['total_read']}")
    print(f"Unique crosslinks processed:  {stats['total_unique']}")
    print(f"Satisfied (< {stats['cutoff']:.1f} Å):       {stats['satisfied']}")
    print(f"Violated (>= {stats['cutoff']:.1f} Å):      {stats['violated']}")
    print(f"Missing atoms:                {stats['missing']}")
    
    if stats['mean_distance'] is not None:
        print(f"\nDistance Range: {stats['min_distance']:.2f} - {stats['max_distance']:.2f} Å")
        print(f"Mean Distance:  {stats['mean_distance']:.2f} Å")
        print(f"Median Distance: {stats['median_distance']:.2f} Å")
        print(f"Std Dev:        {stats['stdev_distance']:.2f} Å")
        print(f"\nSatisfaction Rate: {100.0 * stats['satisfied'] / (stats['satisfied'] + stats['violated']):.1f}%")
        print(f"Violation Rate:    {100.0 * stats['violated'] / (stats['satisfied'] + stats['violated']):.1f}%")
    
    print("=" * 60 + "\n")

# Main function to run in ChimeraX
def visualize_crosslinks(session, cutoff=30.0):
    # Load the PDB structure
    load_structure(session, PDB_FILE)
    
    # Get the currently opened atomic structure model
    pdb_model = None
    for model in session.models:
        if isinstance(model, AtomicStructure):
            pdb_model = model
            break
    
    if pdb_model is None:
        raise ValueError("No atomic structure model found.")
    
    # Apply cartoon representation and color to the protein
    run(session, 'cartoon')
    run(session, 'color khaki cartoons')
    
    # Read the crosslink data
    crosslinks = read_crosslink_data(CROSSLINK_FILE)
    
    # Plot crosslinks on the structure and collect statistics
    stats, distances = plot_crosslinks(session, pdb_model, crosslinks, cutoff=cutoff)
    
    # Print and write statistics
    print_statistics(stats)
    #write_statistics(stats, distances, OUTPUT_STATS_FILE)
    #print(f"\nStatistics written to: {OUTPUT_STATS_FILE}")

# Example of usage:
# In ChimeraX: visualize_crosslinks(session)
# In ChimeraX with custom cutoff: visualize_crosslinks(session, cutoff=35.0)
