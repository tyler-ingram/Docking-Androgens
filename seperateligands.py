from Bio.PDB import PDBParser, PDBIO, Select
import os
"""
This script extracts individual ligands from a PDB file and saves them as separate PDB files.
It assumes that the input PDB file contains multiple ligands and that each ligand is represented by a unique residue name, chain ID, and residue sequence number.
Run this script after running the separate.py script to separate the ligands from the PDB file."""
class SingleLigandSelect(Select):
    def __init__(self, target_residue):
        self.target_residue = target_residue

    def accept_residue(self, residue):
        # Check if the residue is the target ligand
        return residue == self.target_residue

# Load structure
parser = PDBParser(QUIET=True)


# Prepare output directory
input_dir = "Dataset"
output_dir = "ligand_pdbs"
os.makedirs(output_dir, exist_ok=True)

# Want unique ligands
ligands_seen = set()
io = PDBIO()
# Loop through all PDB files in the input directory
for filename in os.listdir(input_dir):
    filepath = os.path.join(input_dir, filename)
    try:
        structure = parser.get_structure(filename, filepath)
    except Exception as e:
        print(f"Failed to parse {filename}: {e}")
        continue
    #Try to extract ligands from the structure
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is a ligand (not water or standard residue)
                hetfield, resseq, icode = residue.id
                resname = residue.resname.strip()
                if hetfield != " " and resname != "HOH": 
                    #Make sure no duplicates are saved
                    ligand_key = (resname)
                    if ligand_key in ligands_seen:
                        continue
                    ligands_seen.add(ligand_key)
                    # Build filename and write file
                    outfile = f"{resname}" + ".pdb"
                    filepath = os.path.join(output_dir, outfile)
                    # Save the ligand
                    io.set_structure(structure)
                    io.save(filepath, SingleLigandSelect(residue))

print(f"Saved {len(ligands_seen)} ligands to folder '{output_dir}'")
