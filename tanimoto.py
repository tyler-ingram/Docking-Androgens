from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
import os

"""
This script calculates the Tanimoto similarity between a reference ligand and a set of ligands in a specified folder.
It uses the RDKit library to generate fingerprints and compute the similarity. This is meant to be used after running run.py and ligandtomol.py.
"""

folder_path = './ligand_mols'
output_path = './results'
os.makedirs(output_path, exist_ok=True)
#DHT is the reference ligand
mol1 = Chem.MolFromMolFile('ligand_mols/DHT_A_931.mol')
#Loop through all the ligands in the folder
for file in os.listdir(folder_path):
    mol2 = Chem.MolFromMolFile(os.path.join(folder_path, file))
    if mol2 is not None:
        # Generate fingerprints
        fp1 = AllChem.GetMorganFingerprint(mol1, 2)
        fp2 = AllChem.GetMorganFingerprint(mol2, 2)
        # Calculate Tanimoto similarity
        tanimoto_similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
        print(f'Tanimoto similarity between {file} and DHT_A_931: {tanimoto_similarity:.4f}')
        prefix = file.split('.')[0]
        #Append the results to a log file after doing the docking
        open(f'./results/{prefix}_log.txt', 'a').write(f'\nTanimoto index: {tanimoto_similarity:.4f}\n')
    else:
        print(f'Failed to read {file}')

def tanimoto_similarity(id1, id2):
    """
    Calculate the Tanimoto similarity between two molecules.
    """
    mol1 = Chem.MolFromMolFile(f'ligand_mols/{id1}.mol')
    mol2 = Chem.MolFromMolFile(f'ligand_mols/{id2}.mol')
    fp1 = AllChem.GetMorganFingerprint(mol1, radius=2)
    fp2 = AllChem.GetMorganFingerprint(mol2,radius=2)
    return DataStructs.TanimotoSimilarity(fp1, fp2)