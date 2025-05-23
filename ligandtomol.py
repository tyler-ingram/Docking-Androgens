import os
import subprocess
"""This script converts PDB files to MOL files using Open Babel.
MOL files are needed for calculating Tanimoto similarity using RDKit.
This script should be run after running the separate.py and seperateligands.py scripts.
"""
input_dir = "ligand_pdbs"
output_dir = "ligand_mols"
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    if filename.lower().endswith(".pdb"):
        input_path = os.path.join(input_dir, filename)
        output_filename = os.path.splitext(filename)[0] + ".mol"
        output_path = os.path.join(output_dir, output_filename)

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
