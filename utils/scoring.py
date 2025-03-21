import numpy as np
from rdkit.Chem import AllChem
from tqdm import tqdm

def best_o3a_match(query_confs: list, target_mol, verbose: bool = False):
    """
    Compare each query conformer against the target molecule using O3A scoring,
    and return the best matching conformer and its score.
    """
    scores = []
    for i in tqdm(range(len(query_confs)), disable=not verbose):
        o3a = AllChem.GetO3A(query_confs[i], target_mol)
        scores.append(o3a.Score())
    best_index = int(np.argmax(scores))
    return query_confs[best_index], scores[best_index]