#!/usr/bin/env python3
"""
Extract protein and RNA structures from mmCIF files and generate topology files.

This script processes mmCIF and FASTA files to extract protein and RNA chains,
save them as separate PDB files, and create a topology file compatible with
integrative modeling platforms.
"""

import os
import re
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Set
from collections import defaultdict

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select, Structure

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

# Constants
RNA_RESIDUE_NAMES = {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "ADE", "CYT", "GUA", "THY", "URI"}
RIBOSOMAL_PATTERN = re.compile(r"ribosomal protein ([SL])(\d+)([A-Za-z]*)", re.IGNORECASE)

class PDBChainRemapper:
    """Handles remapping of multi-character chain IDs to PDB-compatible single-character IDs."""
    
    def __init__(self):
        self.chain_map = {}
        self.next_id = 0
        self.id_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    
    def get_pdb_chain_id(self, original_id):
        """Map an original chain ID to a PDB-compatible single character ID."""
        # Handle tuple input (some Bio.PDB versions return chain.id as tuple)
        if isinstance(original_id, tuple):
            original_id = original_id[0]
            
        if original_id not in self.chain_map:
            # If the ID is already PDB compatible and not used, keep it
            if (len(str(original_id)) == 1 and 
                str(original_id) in self.id_chars and 
                str(original_id) not in self.chain_map.values()):
                self.chain_map[original_id] = str(original_id)
            else:
                # Assign new ID
                while self.next_id < len(self.id_chars):
                    new_id = self.id_chars[self.next_id]
                    self.next_id += 1
                    if new_id not in self.chain_map.values():
                        self.chain_map[original_id] = new_id
                        break
                else:
                    raise ValueError(f"Ran out of available chain IDs. Current map: {self.chain_map}")
        
        return self.chain_map[original_id]
    
class ChainRenamer(Select):
    """Selector that also renames chains to be PDB compatible."""
    
    def __init__(self, chain_ids, chain_remapper):
        self.chain_ids = {str(cid) if isinstance(cid, tuple) else cid for cid in chain_ids}
        self.remapper = chain_remapper
        self.chain_map = {}  # Maps original ID to new ID
    
    def accept_chain(self, chain):
        chain_id = str(chain.id) if isinstance(chain.id, tuple) else chain.id
        return chain_id in self.chain_ids
    
    def get_atom(self, atom):
        chain = atom.get_parent().get_parent()
        chain_id = str(chain.id) if isinstance(chain.id, tuple) else chain.id
        
        if chain_id not in self.chain_map:
            try:
                self.chain_map[chain_id] = self.remapper.get_pdb_chain_id(chain_id)
                chain.id = self.chain_map[chain_id]
            except Exception as e:
                logger.error(f"Failed to remap chain ID {chain_id}: {e}")
                raise
                
        return atom

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Extract protein/RNA structures from mmCIF files.")
    
    parser.add_argument("--mmcif", required=True, 
                        help="Path to the mmCIF file")
    parser.add_argument("--fasta", required=True, 
                        help="Path to the FASTA file")
    parser.add_argument("--output-dir", default="output_files",
                        help="Directory for output files (default: output_files)")
    parser.add_argument("--input-dir", default="input_files",
                        help="Directory containing input files (default: input_files)")
    parser.add_argument("--rna-fasta", default="RNA.fasta",
                        help="Filename for RNA FASTA reference (default: RNA.fasta)")
    parser.add_argument("--protein-fasta", default="ribosome.fasta",
                        help="Filename for protein FASTA reference (default: ribosome.fasta)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable verbose output")
    
    return parser.parse_args()


class ChainSelector(Select):
    """Select specific chains for output."""
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids


def standardize_molecule_name(name):
    """Convert molecule names to standard format."""
    # Process ribosomal protein names (e.g., "ribosomal protein S3a" to "RPS3A")
    match = RIBOSOMAL_PATTERN.search(name)
    if match:
        subunit = match.group(1)  # S or L
        number = match.group(2)   # Number
        suffix = match.group(3)   # Optional letter suffix
        return f"RP{subunit}{number}{suffix}"
    else:
        # Return original name for non-ribosomal proteins
        return name

def determine_molecule_type(chain):
    """Determine if a chain is protein or RNA based on residue composition."""
    rna_count = 0
    other_count = 0
    
    for residue in chain:
        if residue.id[0] == " ":  # Standard residue (not water, etc.)
            if residue.resname.strip() in RNA_RESIDUE_NAMES:
                rna_count += 1
            else:
                other_count += 1
    
    # If more than 80% of residues are RNA, classify as RNA
    if rna_count > 0 and rna_count >= 0.8 * (rna_count + other_count):
        return "rna"
    else:
        return "protein"


def parse_fasta_headers(fasta_file):
    """Extract chain IDs and molecule names from FASTA headers."""
    chain_info = {}
    
    try:
        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    try:
                        parts = line.strip().split("|")
                        
                        # Extract chain ID - handle different formats
                        chain_part = parts[1].strip()
                        auth_match = re.search(r'\[auth\s+([^\]]+)\]', chain_part)
                        if auth_match:
                            chain_id = auth_match.group(1)
                        else:
                            chain_id = chain_part.split()[1]
                        
                        # Extract molecule name
                        molecule_name = parts[2].strip()
                        
                        # Standardize name for ribosomal proteins
                        std_name = standardize_molecule_name(molecule_name)
                        
                        chain_info[chain_id] = {
                            'molecule_name': std_name,
                            'original_name': molecule_name,
                            'species': parts[3] if len(parts) > 3 else ""
                        }
                        
                    except (IndexError, ValueError) as e:
                        logger.warning(f"Failed to parse FASTA header: {line.strip()} - {e}")
    except Exception as e:
        logger.error(f"Error reading FASTA file: {e}")
        raise
        
    return chain_info


def extract_residue_ranges(chain):
    """Extract continuous residue ranges from a chain."""
    residue_list = [res.id[1] for res in chain if res.id[0] == " "]
    if not residue_list:
        return []
    
    ranges = []
    start = residue_list[0]
    prev = start
    
    for res_id in residue_list[1:]:
        if res_id > prev + 1:  # Gap detected
            ranges.append((start, prev))
            start = res_id
        prev = res_id
        
    ranges.append((start, prev))  # Add the last range
    return ranges


def generate_topology_entry(molecule_name, molecule_type, chain_id, res_range, pdb_file, 
                           fasta_file, fasta_id):
    """Generate a properly formatted topology entry line."""
    start, end = res_range
    
    # Set molecule-specific parameters
    if molecule_type == "rna":
        color = "dim gray"
        bead_size = "10"
        # For RNA, use 'U' as chain ID in topology file
        chain_id_display = "U"
    else:
        color = "medium blue"
        bead_size = "1"
        chain_id_display = chain_id
    
    # Format the topology line with proper spacing
    entry = (
        f"|{molecule_name:<7}|{color:<14}|{fasta_file:<20}|{fasta_id:<9}|"
        f"{pdb_file:<20}|{chain_id_display}|{start},{end:<6}|   |"
        f"{bead_size:<3}|  |1 | | ||"
    )
    
    return entry

# Add this new function
def rename_chains_in_structure(structure, chain_ids_to_keep, remapper):
    """
    Directly rename chains in a structure to be PDB-compatible.
    This must be done before saving the structure to avoid PDB format errors.
    
    Args:
        structure: Bio.PDB Structure object
        chain_ids_to_keep: Set of chain IDs to keep
        remapper: PDBChainRemapper instance
    
    Returns:
        dict: Mapping of old chain IDs to new chain IDs
    """
    chain_map = {}
    
    for model in structure:
        chains_to_process = []
        # First collect all chains to avoid modification during iteration
        for chain in model:
            chain_id = str(chain.id) if isinstance(chain.id, tuple) else chain.id
            if chain_id in chain_ids_to_keep:
                chains_to_process.append(chain)
        
        # Now process the chains
        for chain in chains_to_process:
            chain_id = str(chain.id) if isinstance(chain.id, tuple) else chain.id
            new_id = remapper.get_pdb_chain_id(chain_id)
            chain_map[chain_id] = new_id
            chain.id = new_id
    
    return chain_map

# Add these functions right after the rename_chains_in_structure function
def count_atoms(structure, chain_ids):
    """Count atoms in the specified chains."""
    atom_count = 0
    for model in structure:
        for chain in model:
            chain_id = str(chain.id) if isinstance(chain.id, tuple) else chain.id
            if chain_id in chain_ids:
                atom_count += sum(1 for _ in chain.get_atoms())
    return atom_count

def distribute_chains_to_files(structure, chains, max_atoms_per_file=65000, max_chains_per_file=40):
    """Distribute chains to multiple files based on atom count and chain count limits."""
    chain_groups = []
    current_group = []
    current_atom_count = 0
    current_chain_count = 0
    
    chain_atom_counts = [(chain_id, count_atoms(structure, [chain_id])) for chain_id in chains]
    sorted_chains = sorted(chain_atom_counts, key=lambda x: x[1], reverse=True)
    
    logger.info(f"Total atoms in all chains: {sum(atom_count for _, atom_count in chain_atom_counts)}")
    
    standalone_chains = []
    remaining_chains = []
    
    for chain_id, atom_count in sorted_chains:
        if atom_count > max_atoms_per_file / 2:
            standalone_chains.append((chain_id, atom_count))
            logger.info(f"Chain {chain_id} has {atom_count} atoms - will be placed in its own file")
        else:
            remaining_chains.append((chain_id, atom_count))
    
    for chain_id, _ in standalone_chains:
        chain_groups.append([chain_id])
    
    for chain_id, atom_count in remaining_chains:
        if (current_atom_count + atom_count > max_atoms_per_file or 
            current_chain_count >= max_chains_per_file):
            if current_group:
                chain_groups.append(current_group)
            current_group = []
            current_atom_count = 0
            current_chain_count = 0
        
        current_group.append(chain_id)
        current_atom_count += atom_count
        current_chain_count += 1
    
    if current_group:
        chain_groups.append(current_group)
    
    for i, group in enumerate(chain_groups):
        group_atom_count = sum(count_atoms(structure, [chain_id]) for chain_id in group)
        logger.info(f"Group {i+1}: {len(group)} chains, {group_atom_count} atoms")
    
    return chain_groups

class PDBAtomRenumberer(Select):
    """Renumber atoms sequentially when writing a PDB file."""
    
    def __init__(self, chain_ids=None, max_atoms=99999):
        self.chain_ids = set(chain_ids) if chain_ids is not None else None
        self.atom_count = 0
        self.max_atoms = max_atoms
        self.skipped_atoms = 0
    
    def accept_chain(self, chain):
        chain_id = chain.id if not isinstance(chain.id, tuple) else chain.id[0]
        return self.chain_ids is None or chain_id in self.chain_ids
    
    def accept_atom(self, atom):
        chain = atom.get_parent().get_parent()
        if self.accept_chain(chain):
            if self.atom_count >= self.max_atoms:
                self.skipped_atoms += 1
                return False
            
            self.atom_count += 1
            atom.serial_number = self.atom_count
            return True
        return False
    
class CombinedPDBWriter(Select):
    """Both select chains and renumber atoms within PDB format limits."""
    
    def __init__(self, chain_ids, max_atoms=65000):
        self.chain_ids = set(chain_ids)
        self.atom_count = 0
        self.max_atoms = max_atoms
        self.skipped_atoms = 0
        self.processed_chains = set()
        
    def accept_chain(self, chain):
        return chain.id in self.chain_ids
        
    def accept_residue(self, residue):
        chain = residue.get_parent()
        return self.accept_chain(chain)
        
    def accept_atom(self, atom):
        chain = atom.get_parent().get_parent()
        if self.accept_chain(chain):
            if self.atom_count + 1 > self.max_atoms:
                self.skipped_atoms += 1
                return False
            
            chain_id = chain.id
            if chain_id not in self.processed_chains and self.atom_count > 0.85 * self.max_atoms: #reduced to 85%
                self.skipped_atoms += 1
                return False
            
            self.processed_chains.add(chain_id)
            self.atom_count += 1
            atom.serial_number = self.atom_count
            return True
        return False
    
def extract_protein_structures(mmcif_file, fasta_file, output_dir, rna_fasta="RNA.fasta", 
                              protein_fasta="ribosome.fasta"):
    """
    Extract protein and RNA structures from mmCIF file and create output files.
    Split proteins into multiple files if needed to handle PDB chain ID limits.
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse mmCIF file
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", mmcif_file)
        model = structure[0]  # Use first model
        logger.info(f"Successfully parsed mmCIF file: {mmcif_file}")
    except Exception as e:
        logger.error(f"Failed to parse mmCIF file: {e}")
        raise
    
    # Parse FASTA headers
    try:
        chain_info = parse_fasta_headers(fasta_file)
        logger.info(f"Found {len(chain_info)} chains in FASTA file")
    except Exception as e:
        logger.error(f"Failed to parse FASTA file: {e}")
        raise
    
    # Process chains
    protein_chains = []  # Use list to maintain original order
    rna_chains = []
    chain_mapping = {}
    
    # First pass: determine molecule types
    for old_chain_id, info in chain_info.items():
        try:
            if isinstance(old_chain_id, tuple):
                old_chain_id = old_chain_id[0]
                
            try:
                chain = model[old_chain_id]
            except KeyError:
                logger.warning(f"Chain {old_chain_id} not found in mmCIF file")
                continue
                
            # Determine if this is RNA or protein
            molecule_type = determine_molecule_type(chain)
            molecule_name = info['molecule_name']
            
            # Store chain mapping information
            chain_mapping[old_chain_id] = {
                'molecule_name': molecule_name,
                'molecule_type': molecule_type,
                'ranges': extract_residue_ranges(chain)
            }
            
            # Add to appropriate list
            if molecule_type == "rna":
                rna_chains.append(old_chain_id)
            else:
                protein_chains.append(old_chain_id)
                
        except Exception as e:
            logger.warning(f"Error processing chain {old_chain_id}: {e}")
    
    logger.info(f"Found {len(protein_chains)} protein chains and {len(rna_chains)} RNA chains")
    
    # Group chains by atom count and chain count
    logger.info("Calculating atom counts for chain distribution...")
    protein_chain_groups = distribute_chains_to_files(structure, protein_chains)
    logger.info(f"Split protein chains into {len(protein_chain_groups)} files")    

    # Process RNA chains
    rna_pdb_file = os.path.join(output_dir, "rna.pdb")
    if rna_chains:
        from copy import deepcopy
        rna_structure = deepcopy(structure)
        chain_remapper = PDBChainRemapper()  # Fresh remapper for RNA
        
        rna_chain_map = rename_chains_in_structure(rna_structure, set(rna_chains), chain_remapper)
        
        io = PDBIO()
        io.set_structure(rna_structure)
        io.save(rna_pdb_file, ChainSelector(rna_chain_map.values()))
        
        # Update chain mapping
        for old_id, new_id in rna_chain_map.items():
            if old_id in chain_mapping:
                chain_mapping[old_id]['new_id'] = new_id
                chain_mapping[old_id]['pdb_file'] = "rna.pdb"
        
        logger.info(f"Saved RNA chains to {rna_pdb_file} with PDB-compatible IDs")
    
    # Process protein chains in groups
    protein_pdb_files = []
    for group_idx, chain_group in enumerate(protein_chain_groups):
        protein_pdb_file = os.path.join(output_dir, f"proteins_{group_idx+1}.pdb")
        protein_pdb_files.append(os.path.basename(protein_pdb_file))
        
        from copy import deepcopy
        protein_structure = deepcopy(structure)
        chain_remapper = PDBChainRemapper()  # Fresh remapper for each group
        
        protein_chain_map = rename_chains_in_structure(protein_structure, set(chain_group), chain_remapper)
        
        io = PDBIO()
        io.set_structure(protein_structure)
        combined_writer = CombinedPDBWriter(chain_ids=protein_chain_map.values(), max_atoms=75000)
        io.save(protein_pdb_file, combined_writer)

        if combined_writer.skipped_atoms > 0:
            logger.warning(f"WARNING: Skipped {combined_writer.skipped_atoms} atoms due to PDB format limit!")

        logger.info(f"Saved {combined_writer.atom_count} atoms to {protein_pdb_file}")

        # Update chain mapping with only processed chains
        for old_id, new_id in protein_chain_map.items():
            if old_id in chain_mapping and new_id in combined_writer.processed_chains:
                chain_mapping[old_id]['new_id'] = new_id
                chain_mapping[old_id]['pdb_file'] = os.path.basename(protein_pdb_file)
    
    # Log the chain ID mapping
    logger.info("Chain ID mapping for PDB compatibility:")
    for old_id, new_data in chain_mapping.items():
        if 'new_id' in new_data:
            pdb_file = new_data.get('pdb_file', 'unknown')
            logger.info(f"  {old_id} → {new_data['new_id']} in {pdb_file} ({new_data['molecule_name']})")
    
    # Create topology file
    topology_file = os.path.join(output_dir, "topology.dat")
    with open(topology_file, "w") as f:
        # Write header
        f.write("|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|\n")
        f.write("|\n")
        
        # Process RNA chains first
        for old_chain_id in rna_chains:
            info = chain_mapping[old_chain_id]
            if 'new_id' not in info:
                continue  # Skip if no new ID was assigned
                
            new_chain_id = info['new_id']
            molecule_name = info['molecule_name']
            
            # For RNA entries, use RNA1 for ribosomal proteins
            rna_name = "RNA1" if molecule_name.startswith("RP") else molecule_name
            
            for start, end in info['ranges']:
                entry = generate_topology_entry(
                    rna_name, "rna", new_chain_id, 
                    (start, end), "rna.pdb", rna_fasta, "18SRNA,RNA"
                )
                f.write(f"{entry}\n\n")
        
        # Process protein chains
        for old_chain_id in protein_chains:
            info = chain_mapping[old_chain_id]
            if 'new_id' not in info or 'pdb_file' not in info:
                continue  # Skip if no new ID was assigned
                
            new_chain_id = info['new_id']
            molecule_name = info['molecule_name']
            pdb_file = info['pdb_file']
            
            for start, end in info['ranges']:
                entry = generate_topology_entry(
                    molecule_name, "protein", new_chain_id, 
                    (start, end), pdb_file, protein_fasta, molecule_name
                )
                f.write(f"{entry}\n")
    
    logger.info(f"Created topology file: {topology_file}")
    return chain_mapping

def main():
    """Main function to process files and generate output."""
    args = parse_arguments()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    try:
        # Resolve file paths
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir)
        
        mmcif_path = input_dir / args.mmcif if not os.path.isabs(args.mmcif) else args.mmcif
        fasta_path = input_dir / args.fasta if not os.path.isabs(args.fasta) else args.fasta
        
        # Check if input files exist
        if not os.path.exists(mmcif_path):
            logger.error(f"Input mmCIF file not found: {mmcif_path}")
            sys.exit(1)
        
        if not os.path.exists(fasta_path):
            logger.error(f"Input FASTA file not found: {fasta_path}")
            sys.exit(1)
            
        # Extract structures and create topology file
        chain_mapping = extract_protein_structures(
            mmcif_path, fasta_path, output_dir,
            args.rna_fasta, args.protein_fasta
        )
        
        logger.info("Processing completed successfully")
        
    except Exception as e:
        logger.error(f"Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()