import os
import subprocess
from pocket_definition import get_pocket_coordinates
"""
This script runs the docking process using QVina2.
Run this script after running the separate.py, seperateligands.py, and pdbtopdbqt.py scripts.
"""

input_dir = "ligand_pdbqt"
output_dir = "results"
receptor = "Project1_Data/protein.pdbqt"
pdb_file = "Project1_Data/protein.pdb"  
#Got from coach for their rank 1 pocket for protein.pdb
#See the metadata folder to see all pockets
residue_ids = [33, 36, 37, 40, 43, 74, 77, 81, 84, 96, 112, 115, 116, 119, 205, 208, 209, 212, 227]
residue_ids = [r + 668 for r in residue_ids]  # Adjust for the offset coach starts from 1 while our pdb starts from 669
#Get the box of the pocket for docking
center, size = get_pocket_coordinates(pdb_file, residue_ids)
os.makedirs(output_dir, exist_ok=True)
# Loop through all the PDBQT files in the input directory and dock them into the receptor
for filename in os.listdir(input_dir):
    if filename.lower().endswith(".pdbqt"):
        input_path = os.path.join(input_dir, filename)
        output_filename = os.path.splitext(filename)[0] + "_out.pdbqt"
        log_filename = os.path.splitext(filename)[0] + "_log.txt"
        output_path = os.path.join(output_dir, output_filename)
        log_path = os.path.join(output_dir, log_filename)
        # Dock the ligand using QVina2
        # The output file will be saved in the results folder
        # The log file will be saved in the results folder
        try:
            result = subprocess.run(
                ["qvina2", "--receptor", receptor, "--ligand", input_path, "--out", output_path,
                 "--log", log_path,
                 "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
                 "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2])],
                capture_output=True, text=True
            )
            if result.returncode != 0:
                raise RuntimeError(result.stderr.strip())
        except Exception as e:
            print(f"Failed to dock {filename}: {e}")