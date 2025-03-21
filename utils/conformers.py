from rdkit import Chem
from rdkit.Chem import AllChem

def generate_conformers(smiles: str, nConfs: int = 25, maxIters: int = 500, verbose: bool = False) -> Chem.Mol:
    """
    Generate conformers for a molecule specified by its SMILES string.
    Returns a molecule with multiple conformers.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = AllChem.AddHs(mol)
    conf_ids = AllChem.EmbedMultipleConfs(
        mol, nConfs, randomSeed=42, 
        useRandomCoords=False, enforceChirality=True,
        useExpTorsionAnglePrefs=True, useBasicKnowledge=True,
        useSmallRingTorsions=True, useMacrocycleTorsions=True
    )
    failed = 0
    for cid in conf_ids:
        opt = AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=maxIters)
        if opt == -1 and verbose:
            print(f'Force field could not be set up for conformer {cid}!')
        failed += opt
    if verbose:
        print(f'{failed} conformer minimizations failed to converge')
    return mol

def align_conformers(multi_conf_mol: Chem.Mol) -> list:
    """
    Align conformers within a molecule and return them as a list of aligned molecules.
    """
    AllChem.AlignMolConformers(multi_conf_mol)
    aligned_mols = []
    for conf in multi_conf_mol.GetConformers():
        block = Chem.MolToMolBlock(multi_conf_mol, confId=conf.GetId())
        aligned_mol = Chem.MolFromMolBlock(block, removeHs=False)
        aligned_mols.append(aligned_mol)
    return aligned_mols
