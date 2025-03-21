# Drug Overlapper

Drug Overlapper is a modular pipeline for analyzing molecules from an SDF file (e.g., DrugBank). It generates conformers for a query molecule (via a provided SMILES string), computes overlap (O3A) scores against a database of molecules, and exports the best overlapping results as both an HTML file (with embedded molecule images) and a pickle file containing the processed DataFrame.

## Features

- **Data Loading:** Reads molecules from an SDF file.
- **Conformer Generation:** Generates and aligns multiple conformers for a query molecule.
- **Scoring:** Computes the best overlap (O3A score) between the query molecule conformers and database molecules.
- **Export:** Saves the results as an HTML file (with images) and as a pickle file for later use.
- **Modular Design:** Organized into separate modules for data handling, conformers, scoring, and export for ease of maintenance and extension.

## Directory Structure

```bash
├── data
├── main.py
├── readme.md
└── utils
    ├── conformers.py
    ├── data.py
    ├── export.py
    └── scoring.py
```
## Installation

It is recommended to use a virtual environment. For example:

```bash
python -m venv venv
source venv/bin/activate  # On Windows, use venv\Scripts\activate
pip install requirements.txt
```

## Usage

Before running the pipeline, ensure that your SDF file (e.g., drugbank.sdf) is placed in the data folder.

### Command-Line Arguments

The script accepts the following command-line arguments:

### Command-Line Arguments

- `--query`: (Required) The SMILES string of the query molecule.
- `--sdf`: (Optional) Path to the input SDF file (default: data/drugbank.sdf).
- `--conformers`: (Optional) Number of conformers to generate (default: 25).
- `--run_name`: (Optional) Run name for the output files. Results will be saved under outputs/{run_name} (default: out).
- `--top`: (Optional) Number of top scoring molecules to export (default: 250).

### Example Command

```bash
python main.py --query "CC(=O)OC1=CC=CC=C1C(=O)O" --sdf "data/myfile.sdf" --conformers 10 --run_name "aspirin" --top 100
```

This command will:

- Generate 10 conformers for the query molecule.
- Load molecules from the SDF file at "data/myfile.sdf".
- Compute and sort molecules based on the best O3A score.
- Export the top 100 molecules as:
    - An HTML file ("overlaps.html") with embedded molecule images.
    - A pickle file ("dataframe.pkl") containing the processed DataFrame.

Both files will be saved in the folder: outputs/my_run/.

### Loading the Exported Pickle 

You can load the exported DataFrame later using the helper function in `utils/data.py`:

```python
from export.data import load_pickle 

df = load_pickle("outputs/my_run/results.pkl")
print(df.head())
``` 

## Contributing 

Contributions and suggestions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.