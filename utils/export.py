import os
import base64
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import HTML
import pandas as pd

def smiles_to_img_base64(smiles: str, size=(400, 400)) -> str:
    """
    Convert an RDKit molecule from a SMILES string to a base64-encoded PNG image embedded in an HTML img tag.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.RemoveHs(mol)
    img = Draw.MolToImage(mol, size=size)
    buffer = BytesIO()
    img.save(buffer, format="PNG")
    img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
    return f'<img src="data:image/png;base64,{img_str}" width="{size[0]}" height="{size[1]}"/>'

def export_to_html(df, run_id: str, filename: str = 'overlaps.html'):
    """
    Export the DataFrame as an HTML file with embedded molecule images to the directory outputs/{run_id}/filename.
    """
    output_dir = os.path.join("outputs", run_id)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    html_content = HTML(df.to_html(escape=False)).data
    with open(output_path, 'w') as f:
        f.write(html_content)
    print(f"HTML file exported to {output_path}")

def export_to_pickle(df, run_id: str, filename: str = 'dataframe.pkl'):
    """
    Export the DataFrame as a pickle file to the directory outputs/{run_id}/filename.
    """
    output_dir = os.path.join("outputs", run_id)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    df.to_pickle(output_path)
    print(f"DataFrame exported to {output_path}")

def load_pickle(filepath: str) -> pd.DataFrame:
    """
    Load and return a DataFrame from a pickle file.
    """
    return pd.read_pickle(filepath)
