import IMP
import IMP.atom
import IMP.rmf
import RMF
import os
import numpy as np


def rmf_to_pdbs(rmf_in, output_dir=None, frames=None, sel_state=0):
    """
    Convert RMF frames to PDB files.
    
    Parameters:
    - rmf_in: Path to input RMF file
    - output_dir: Output directory (default: pdbs/{rmf_name}/)
    - frames: List/array of frame indices to export (default: all frames)
    - sel_state: State to export (default: 0)
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
    hier = IMP.rmf.create_hierarchies(f, m)[0]
    
    # Select specific state if needed
    states = IMP.atom.get_by_type(hier, IMP.atom.STATE_TYPE)
    if states:
        print(f"Found {len(states)} states, using state {sel_state}")
        hier_to_write = states[sel_state]
    else:
        hier_to_write = hier
    
    # Determine frames to export
    if frames is None:
        frames = range(n_frames)
    
    print(f"Exporting {len(frames)} frames to {output_dir}")
    
    # Export each frame
    for i in frames:
        if i % 100 == 0:
            print(f"Processing frame {i}/{n_frames}")
        
        # Load frame
        IMP.rmf.load_frame(f, RMF.FrameID(i))
        
        # Write PDB
        pdb_path = os.path.join(output_dir, f"frame_{i}.pdb")
        IMP.atom.write_pdb(hier_to_write, pdb_path)
    
    print(f"Done! Exported {len(frames)} PDB files to {output_dir}")
    
    del f, m


# Example usage
if __name__ == "__main__":
    # Export all frames
    #rmf_to_pdbs("./input_data/A_models_cluster2_0_aligned.rmf3")
    
    # Export every 100th frame
    rmf_to_pdbs("./input_data/A_models_cluster2_0_aligned.rmf3", frames=range(0, 10000, 100))

    # Custom output directory
    # rmf_to_pdbs("./input_data/A_models_cluster2_0_aligned.rmf3", output_dir="my_pdbs/run1")