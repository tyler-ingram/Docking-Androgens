from Bio.PDB import PDBParser, PDBIO, Select
"""
This script extracts the protein and ligands from a PDB file and saves them as separate PDB files.
It assumes that the input PDB file contains a protein and ligands, and that the ligands are represented by unique residue names.
Run this script first to separate the protein and ligands from the PDB file.
"""
# Get the protein
class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "  # standard residues only
#Get the ligands, which are not water and not standard residues
class LigandSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] != " " and residue.resname != "HOH"

#Get the structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", "1z95.pdb")  # <-- change filename

# Save protein
protein_io = PDBIO()
protein_io.set_structure(structure)
protein_io.save("protein-2.pdb", ProteinSelect())

# Save all ligands
ligand_io = PDBIO()
ligand_io.set_structure(structure)
ligand_io.save("ligands-2.pdb", LigandSelect())

print("Saved: protein.pdb and ligands-2.pdb")
