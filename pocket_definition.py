from Bio.PDB import PDBParser


pdb_file = "./other_pdbs/protein.pdb"  
#Got from coach for their rank 1 pocket
residue_ids = [33, 36, 37, 40, 43, 74, 77, 81, 84, 96, 112, 115, 116, 119, 205, 208, 209, 212, 227]
residue_ids = [r + 668 for r in residue_ids]  # Adjust for the offset coach starts from 1 while our pdb starts from 669

# === Load structure ===
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

# === Collect atom coordinates from target residues ===
coords = []
for model in structure:
    for chain in model:
        for res in chain:
            if res.get_id()[1] in residue_ids:
                for atom in res:
                    coords.append(atom.coord)

# === Calculate bounding box ===
import numpy as np
coords = np.array(coords)
min_coords = coords.min(axis=0)
max_coords = coords.max(axis=0)
center = (min_coords + max_coords) / 2
size = max_coords - min_coords

print(f"Bounding Box Center: {center}")
print(f"Box Size (X,Y,Z): {size}")