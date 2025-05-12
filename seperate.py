from Bio.PDB import PDBParser, PDBIO, Select

# Define Select classes
class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "  # standard residues only

class LigandSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] != " " and residue.resname != "HOH"  # HETATM and not water

# Parse PDB
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
