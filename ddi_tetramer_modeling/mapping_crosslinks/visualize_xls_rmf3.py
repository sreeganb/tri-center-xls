import IMP
import IMP.pmi
import IMP.pmi.tools
import IMP.rmf
import IMP.core
import RMF
import pandas as pd
import numpy as np


class CrosslinkVisualizer:
    """Add crosslink restraints to RMF for ChimeraX visualization"""
    
    def __init__(self, rmf_in, crosslink_csv, rmf_out=None, 
                 threshold=30.0, atom_type='CA'):
        self.rmf_in = rmf_in
        self.crosslink_csv = crosslink_csv
        self.rmf_out = rmf_out or rmf_in.replace('.rmf3', '_with_xls.rmf3')
        self.threshold = threshold
        self.atom_type = atom_type
        
        self.chain_to_molecule = {
            'A': ('DDI1', 0), 'B': ('DDI1', 1),
            'C': ('DDI2', 0), 'D': ('DDI2', 1),
        }
        
    def read_crosslinks(self):
        """Read crosslinks from CSV"""
        df = pd.read_csv(self.crosslink_csv)
        crosslinks = []
        
        for _, row in df.iterrows():
            chain1, res1 = row['Chain1'], int(row['Residue1'])
            chain2, res2 = row['Chain2'], int(row['Residue2'])
            
            if chain1 in self.chain_to_molecule and chain2 in self.chain_to_molecule:
                mol1, copy1 = self.chain_to_molecule[chain1]
                mol2, copy2 = self.chain_to_molecule[chain2]
                crosslinks.append((mol1, res1, copy1, mol2, res2, copy2))
        
        print(f"Loaded {len(crosslinks)} crosslinks from {self.crosslink_csv}")
        return crosslinks
    
    def find_particle(self, hier, molecule, residue, copy_index):
        """Find particle for given molecule, residue, copy"""
        # Try with atom type first
        sel = IMP.atom.Selection(
            hier,
            molecule=molecule,
            residue_index=residue,
            atom_type=IMP.atom.AtomType(self.atom_type),
            copy_index=copy_index
        ).get_selected_particles()
        
        if len(sel) == 1:
            return sel[0]
        
        # Try without atom type (coarse-grained)
        sel = IMP.atom.Selection(
            hier,
            molecule=molecule,
            residue_index=residue,
            copy_index=copy_index,
            resolution=1
        ).get_selected_particles()
        
        return sel[0] if len(sel) >= 1 else None
    
    def get_rmf_node(self, rmf_file, particle):
        """Get RMF node for a particle"""
        # Find the corresponding RMF node for this particle
        return rmf_file.get_node_from_association(particle)
    
    def add_crosslink_restraints(self):
        """Add crosslink restraints to RMF using RMF bonds"""
        
        crosslinks = self.read_crosslinks()
        
        # Load input RMF
        m = IMP.Model()
        rmf_in = RMF.open_rmf_file_read_only(self.rmf_in)
        hier = IMP.rmf.create_hierarchies(rmf_in, m)[0]
        
        # Create output RMF
        rmf_out = RMF.create_rmf_file(self.rmf_out)
        IMP.rmf.add_hierarchy(rmf_out, hier)
        
        # Get RMF associations for particles
        IMP.rmf.load_frame(rmf_in, RMF.FrameID(0))
        
        # Create bond categories for satisfied and violated crosslinks
        bond_factory = RMF.BondFactory(rmf_out)
        
        # Prepare crosslink particle pairs
        crosslink_pairs = []
        added_count = 0
        
        for xl in crosslinks:
            mol1, res1, copy1, mol2, res2, copy2 = xl
            
            # Find particles
            p1 = self.find_particle(hier, mol1, res1, copy1)
            p2 = self.find_particle(hier, mol2, res2, copy2)
            
            if p1 and p2 and p1 != p2:
                crosslink_pairs.append((p1, p2, xl))
                added_count += 1
        
        print(f"Found {added_count} valid crosslink pairs")
        
        # Process each frame and add bonds
        n_frames = rmf_in.get_number_of_frames()
        print(f"Processing {n_frames} frames...")
        
        distances_all_frames = []
        
        for frame_idx in range(0, n_frames, 100):  # Every 100th frame
            if frame_idx % 1000 == 0:
                print(f"  Frame {frame_idx}/{n_frames}")
            
            # Load frame from input
            IMP.rmf.load_frame(rmf_in, RMF.FrameID(frame_idx))
            
            # Calculate distances and create bonds
            frame_distances = []
            
            for p1, p2, xl in crosslink_pairs:
                # Calculate distance
                xyz1 = IMP.core.XYZ(p1)
                xyz2 = IMP.core.XYZ(p2)
                dist = IMP.core.get_distance(xyz1, xyz2)
                frame_distances.append(dist)
                
                # Get RMF nodes for these particles
                try:
                    node1 = rmf_out.get_node_from_association(p1)
                    node2 = rmf_out.get_node_from_association(p2)
                    
                    # Create bond
                    bond = bond_factory.get(rmf_out.add_bond(node1, node2))
                    bond.set_bonded_0(0)  # Bond type
                    bond.set_bonded_1(1)
                    
                except:
                    pass  # Node not found in RMF, skip
            
            distances_all_frames.append(frame_distances)
            
            # Save frame to output
            IMP.rmf.save_frame(rmf_out, str(frame_idx))
        
        # Print distance statistics
        distances_all_frames = np.array(distances_all_frames)
        print(f"\nDistance statistics:")
        print(f"  Mean: {np.mean(distances_all_frames):.2f} Å")
        print(f"  Min: {np.min(distances_all_frames):.2f} Å")
        print(f"  Max: {np.max(distances_all_frames):.2f} Å")
        print(f"  Satisfied (<{self.threshold}Å): {(distances_all_frames < self.threshold).sum()}/{distances_all_frames.size}")
        
        print(f"\nSaved output to: {self.rmf_out}")
        print(f"Total crosslinks added: {added_count}")
        print(f"\nTo visualize in ChimeraX:")
        print(f"1. chimerax {self.rmf_out}")
        print(f"2. Bonds should appear automatically")
        print(f"3. Color by distance: color byattribute bfactor")


def main():
    visualizer = CrosslinkVisualizer(
        rmf_in='input_data/A_models_cluster2_0_aligned.rmf3',
        crosslink_csv='pymol_data/bifunc_pairs.csv',
        rmf_out='trajectory_with_crosslinks.rmf3',
        threshold=30.0
    )
    
    visualizer.add_crosslink_restraints()


if __name__ == '__main__':
    main()