from rdkit import Chem
from rdkit.Chem import DataStructs, AllChem
import os

folder_path = './ligand_mols'
output_path = './results'
os.makedirs(output_path, exist_ok=True)
mol1 = Chem.MolFromMolFile('ligand_mols/DHT_A_931.mol')
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
        open(f'./results/{prefix}_log.txt', 'a').write(f'Tanimoto index: {tanimoto_similarity:.4f}\n')
    else:
        print(f'Failed to read {file}')
