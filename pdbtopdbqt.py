import os
import subprocess
"""This script converts PDB files to pdbqt files using Open Babel.
PDBQT files are needed for docking with Vina.
This script should be run after running the separate.py and seperateligands.py scripts.
"""
#Ligands are in the ligand_pdbs folder
#PDBQT files are saved in the ligand_pdbqt folder
input_dir = "ligand_pdbs"
output_dir = "ligand_pdbqt"
os.makedirs(output_dir, exist_ok=True)
#Loop through all the PDB files in the input directory
for filename in os.listdir(input_dir):
    if filename.lower().endswith(".pdb"):
        input_path = os.path.join(input_dir, filename)
        #Want to save the PDBQT files in the ligand_pdbqt folder
        output_filename = os.path.splitext(filename)[0] + ".pdbqt"
        output_path = os.path.join(output_dir, output_filename)
        #Convert the PDB files to PDBQT files using Open Babel
        try:
            result = subprocess.run(
                ["obabel", input_path, "-O", output_path],
                capture_output=True, text=True
            )
            if result.returncode != 0:
                raise RuntimeError(result.stderr.strip())

            print(f"Converted {filename} → {output_filename}")
        except Exception as e:
            print(f"Failed to convert {filename}: {e}")