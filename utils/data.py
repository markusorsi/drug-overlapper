import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def load_data(sdf_path: str) -> pd.DataFrame:
    """Load molecules from an SDF file and return a DataFrame with SMILES and heavy atom count."""
    supplier = Chem.SDMolSupplier(sdf_path)
    molecules = [Chem.AddHs(mol) for mol in supplier if mol is not None]
    smiles = [Chem.MolToSmiles(mol) for mol in molecules]
    hac = [rdMolDescriptors.CalcNumHeavyAtoms(mol) for mol in molecules]

    df = pd.DataFrame({'smiles': smiles, 'mol': molecules, 'hac': hac})
    # Filter out molecules with less than 5 heavy atoms
    df = df[df['hac'] >= 5]
    return df
