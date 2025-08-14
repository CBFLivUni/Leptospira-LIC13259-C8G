import os
import csv
import argparse
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import itertools
import numpy as np

# https://www.blopig.com/blog/2023/05/checking-your-pdb-file-for-clashing-atoms/ 
# Partially based off above script + lots of help from chatGPT
# Default threshold is 2.14 - 
#    1.70 nm atomic radius for carbon. 1.7 + 1.7 = 3.4. 3.4 * 0.63 clash cut-off from above script

def parse_pdb(filename):
    """Parse a PDB file and extract chain information."""
    parser = PDBParser()
    structure = parser.get_structure('protein', filename)
    chains = list(structure.get_chains())
    return chains

def calculate_distance(atom1, atom2):
    """Calculate distance between two atoms."""
    coord1 = atom1.get_coord()
    coord2 = atom2.get_coord()
    distance = np.linalg.norm(coord1 - coord2)
    return distance

def get_chain_length(chain):
    """Calculate the length of the protein chain."""
    residues = list(chain.get_residues())
    return len(residues)

def check_clashes(chains, threshold=2.14):
    """Check for clashes between c-alpha/backbone atoms of all chains."""
    clashes = 0
    total_length = 0
    for chain in chains:
        total_length += get_chain_length(chain)

    for chain1, chain2 in itertools.combinations(chains, 2):
        for residue1 in chain1:
            for residue2 in chain2:
                for atom1 in residue1:
                    for atom2 in residue2:
                        if atom1.name == 'CA' and atom2.name == 'CA':
                            distance = calculate_distance(atom1, atom2)
                            if distance < threshold:
                                clashes += 1

    clashes_per_100_residues = (clashes / total_length) * 100
    return clashes, clashes_per_100_residues

def main(input_directory, output_csv, threshold):
    # Suppress PDB parser warnings
    warnings.simplefilter("ignore", PDBConstructionWarning)

    # Write CSV header
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['PDB_File', 'Clashes_Count', 'Clashes_per_100_Residues'])

    # Iterate over PDB files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".pdb"):
            pdb_file = os.path.join(input_directory, filename)
            chains = parse_pdb(pdb_file)
            clashes, clashes_per_100_residues = check_clashes(chains, threshold)
            with open(output_csv, 'a', newline='') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow([filename, clashes, clashes_per_100_residues])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check clashes between protein chains in PDB files.')
    parser.add_argument('input_directory', type=str, help='Path to the directory containing the PDB files')
    parser.add_argument('--threshold', type=float, default=2.14, help='Clash distance threshold in angstroms (default: 2.0)')
    parser.add_argument('--output_csv', type=str, default='clashes_summary.csv', help='Output CSV file path (default: clashes_summary.csv)')
    args = parser.parse_args()

    main(args.input_directory, args.output_csv, args.threshold)