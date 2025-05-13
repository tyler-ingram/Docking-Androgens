import os
import subprocess
from pocket_definition import get_pocket_coordinates

input_dir = "ligand_pdbqt"
output_dir = "results"
receptor = "other_pdbs/protein.pdbqt"
pdb_file = "./other_pdbs/protein.pdb"  
#Got from coach for their rank 1 pocket for protein.pdb
residue_ids = [33, 36, 37, 40, 43, 74, 77, 81, 84, 96, 112, 115, 116, 119, 205, 208, 209, 212, 227]
residue_ids = [r + 668 for r in residue_ids]  # Adjust for the offset coach starts from 1 while our pdb starts from 669
#Get the coordinates of the pocket
center, size = get_pocket_coordinates(pdb_file, residue_ids)
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    if filename.lower().endswith(".pdbqt"):
        input_path = os.path.join(input_dir, filename)
        output_filename = os.path.splitext(filename)[0] + "_out.pdbqt"
        log_filename = os.path.splitext(filename)[0] + "_log.txt"
        output_path = os.path.join(output_dir, output_filename)
        log_path = os.path.join(output_dir, log_filename)
        #print(f"Docking {input_path} into {output_filename}")
        #print(receptor)
        #print(input_path)
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

            #print(f"Docked {filename} into {output_filename}")
        except Exception as e:
            print(f"Failed to dock {filename}: {e}")