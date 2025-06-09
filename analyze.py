import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
import os


def get_center(id):
    pdbqt_file = f'./results/{id}_out.pdbqt'
    coords = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    coords = np.array(coords)
    return np.mean(coords, axis=0)

def calculate_centroid_distance(center1, center2):
    """
    Calculate the Euclidean distance between two centers.
    :param center1: The first center as a numpy array.
    :param center2: The second center as a numpy array.
    :return: The Euclidean distance between the two centers.
    """
    return np.linalg.norm(center1 - center2)

def tanimoto_similarity(id1, id2):
    """
    Calculate the Tanimoto similarity between two molecules.
    """
    mol1 = Chem.MolFromMolFile(f'ligand_mols/{id1}.mol')
    mol2 = Chem.MolFromMolFile(f'ligand_mols/DHT.mol')
    fp1 = AllChem.GetMorganFingerprint(mol1, radius=2)
    fp2 = AllChem.GetMorganFingerprint(mol2,radius=2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def get_binding_affinity(id):
    """
    Extract the binding affinity from a PDBQT file.
    :param id: The id of the docked ligand, which corresponds to a PDBQT file.
    :return: The binding affinity as a float.
    """
    pdbqt_file = f'./results/{id}_log.txt'
    if not os.path.exists(pdbqt_file):
        raise FileNotFoundError(f"PDBQT file for {id} does not exist.")
    with open(pdbqt_file) as f:
        for line in f:
            if line.strip().startswith("1 "):
                affinity = float(line.split()[0])
                break
    return affinity