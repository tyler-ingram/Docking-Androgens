from Bio.PDB import PDBParser
import numpy as np
def get_pocket_coordinates(pdb_file, residue_ids):
    """
    Get the coordinates of the specified residues from a PDB file.
    
    Args:
        pdb_file (str): Path to the PDB file.
        residue_ids (list): List of residue IDs to extract coordinates for.
        
    Returns:
        Pair of numpy arrays: The center and size of the pocket.
    """
    # Get the structure from the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # Get the coordinates of the specified residues
    coords = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_id()[1] in residue_ids:
                    for atom in res:
                        coords.append(atom.coord)
    #Get the min and max coordinates for x y and z to create the box
    coords = np.array(coords)
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    # Calculate the center and size of the pocket
    center = (min_coords + max_coords) / 2
    size = max_coords - min_coords

    return center, size