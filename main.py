import argparse
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
tqdm.pandas()

import utils.data as data
import utils.conformers as conformers
import utils.scoring as scoring
import utils.export as export

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Molecule Analysis Pipeline")
    parser.add_argument(
        "--query",
        required=True,
        help="SMILES string of the query molecule (required)."
    )
    parser.add_argument(
        "--conformers",
        type=int,
        default=25,
        help="Number of conformers to generate (default: 25)."
    )
    parser.add_argument(
        "--run_name",
        default="out",
        help="Run name for the output files (default: 'out')."
    )
    parser.add_argument(
        "--top",
        type=int,
        default=250,
        help="Number of top scoring molecules to export (default: 250)."
    )
    parser.add_argument(
        "--sdf",
        default="data/drugbank.sdf",
        help="Path to the input SDF file (default: data/drugbank.sdf)."
    )
    args = parser.parse_args()
    return args.query, args.conformers, args.run_name, args.top, args.sdf

def generate_query_conformers(query_smiles, n_conformers):
    """
    Generate and align conformers for the query molecule.
    
    Parameters:
      - query_smiles (str): SMILES string for the query molecule.
      - n_conformers (int): Number of conformers to generate.
    
    Returns:
      - list: A list of aligned query conformers.
    """
    query_mol = conformers.generate_conformers(query_smiles, nConfs=n_conformers)
    return conformers.align_conformers(query_mol)

def score_molecules(df, query_confs):
    """
    Compute the best overlap (O3A score) for each molecule in the DataFrame.
    
    Parameters:
      - df (DataFrame): DataFrame containing molecules.
      - query_confs (list): List of query molecule conformers.
    
    Returns:
      - DataFrame: Sorted DataFrame with additional 'conformer' and 'score' columns.
    """
    def apply_scoring(row):
        best_conf, score = scoring.best_o3a_match(query_confs, row['mol'], verbose=False)
        return pd.Series({'conformer': best_conf, 'score': score})
    
    df[['conformer', 'score']] = df.progress_apply(apply_scoring, axis=1)
    return df.sort_values('score', ascending=False)

def prepare_export(df, top_n):
    """
    Prepare the DataFrame for export by processing the top scoring molecules.
    
    Parameters:
      - df (DataFrame): DataFrame with scored molecules.
      - top_n (int): Number of top molecules to export.
    
    Returns:
      - DataFrame: Processed DataFrame ready for export.
    """
    df_export = df.head(top_n).copy()
    df_export['mol'] = df_export['mol'].apply(Chem.RemoveHs)
    df_export['smiles'] = df_export['mol'].apply(Chem.MolToSmiles)
    df_export['image'] = df_export['smiles'].apply(export.smiles_to_img_base64)
    return df_export

def main():
    query_smiles, n_conformers, run_name, top_n, sdf_path = parse_arguments()
    df = data.load_data(sdf_path)
    query_confs = generate_query_conformers(query_smiles, n_conformers)
    scored_df = score_molecules(df, query_confs)
    df_export = prepare_export(scored_df, top_n)
    
    export.export_to_html(df_export, run_id=run_name, filename="overlaps.html")
    export.export_to_pickle(df_export, run_id=run_name, filename="dataframe.pkl")

if __name__ == "__main__":
    main()
