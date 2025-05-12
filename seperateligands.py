from Bio.PDB import PDBParser, PDBIO, Select
import os

class SingleLigandSelect(Select):
    def __init__(self, target_residue):
        self.target_residue = target_residue

    def accept_residue(self, residue):
        return residue == self.target_residue

# Load structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", "ligands.pdb")  # <-- your file

# Prepare output directory
output_dir = "separated_ligands"
os.makedirs(output_dir, exist_ok=True)

# Loop through and extract individual ligands
ligands_seen = set()
io = PDBIO()

for model in structure:
    for chain in model:
        for residue in chain:
            hetfield, resseq, icode = residue.id
            resname = residue.resname.strip()
            if hetfield != " " and resname != "HOH":  # It's a ligand (not protein or water)
                ligand_key = (resname, chain.id, resseq, icode)
                if ligand_key in ligands_seen:
                    continue
                ligands_seen.add(ligand_key)

                # Build filename and write file
                filename = f"{resname}_{chain.id}_{resseq}{icode or ''}.pdb"
                filepath = os.path.join(output_dir, filename)

                io.set_structure(structure)
                io.save(filepath, SingleLigandSelect(residue))

print(f"Saved {len(ligands_seen)} ligands to folder '{output_dir}'")
