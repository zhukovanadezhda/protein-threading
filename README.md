# Protein Structure Prediction Using Threading

This project focuses on protein structure prediction using the [threading (fold recognition) method](https://en.wikipedia.org/wiki/Threading_(protein_sequence)). Threading is a computational approach that aligns a target protein sequence with known structures to predict its 3D structure, even when sequence similarity is low. This method identifies candidate templates by considering both sequence and structural similarities, such as predicted secondary structures and solvent accessibility. By mapping the target sequence onto structurally similar templates, threading can accurately predict protein folds, making it a valuable tool in bioinformatics, medicine, and biotechnology.

<p align="center">
  <img align="center" width="800px" 
    src="doc/assets/fold-recognition.png"
    alt="protein">
</p>
<p align="center">
  <i>
    Fig. 1: Protein fold recognition summary.
  </i>
</p>


## Installation

### Clone the repository

```bash
git clone git@github.com:zhukovanadezhda/protein-threading.git
cd protein-threading
```
### Setup the conda environment

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html). Create the `protein-threading` conda environment:

```bash
conda env create -f environment.yml
```

### Load the environment

```bash
conda activate protein-threading
```

> **Note:** to deactivate an active environment, use:
> ```bash
> conda deactivate
> ```

## Usage

> **Note:** Before running the program, ensure the `src/config.py` file is properly configured to set up your working directories.  
> - **To run an example:** Select one of the provided example directories (e.g., `data/example1`).  
> - **To use custom data:** Create your own directory (e.g., `data/your_dir`) with subdirectories for `sequences` and `templates`. Place your data in these folders and update the paths in `src/config.py`.
>   
> Example `src/config.py` modification:
> ```python
> # Directory paths
> TEMPLATES_DIR = 'data/your_dir/templates/'
> SEQUENCES_DIR = 'data/your_dir/sequences/'
> ```

To run the program, use the following command:

```python
python src/main.py [-h] [--sequences SEQUENCES] [--templates TEMPLATES] [--output_file OUTPUT_FILE] \
                   [--jobs JOBS] [--dry_run] [--verbose]
```

## Arguments

| Argument                  | Description                                                   | Default           |
|:-------------------------:|---------------------------------------------------------------|-------------------|
| `-h`                      | Show a help message and exit.                                 |                   |
| `--sequences`             | Comma-separated list of sequence filenames (`.fasta` format). | All files from `SEQUENCES_DIR` from `src/config.py`. |
| `--templates`             | Comma-separated list of template filenames (`.pdb` format).   | All files from `TEMPLATES_DIR` from `src/config.py`. |
| `--output_file`           | Name of the output CSV file.                                  | `results/energy_scores.csv`|
| `--jobs`                  | Number of parallel jobs to run.                               | All cores         |
| `--dry_run`               | If set, only log actions without processing.                  | False (not set)   |
| `--verbose`               | If set, verbose output enabled.                               | False (not set)   |

## Example

#### Input

> `src/config.py` :
> ```python
> # Directory paths
> TEMPLATES_DIR = 'data/example1/templates/'
> SEQUENCES_DIR = 'data/example1/sequences/'
> ```

```python
python src/main.py --sequences 5AWL.fasta --output_file results/example1_result.csv
```

#### Output

```python
2024-09-06 13:30:02,290 - INFO - Loading DOPE score data...
2024-09-06 13:30:02,732 - INFO - Processing sequences and templates...
2024-09-06 13:30:02,732 - INFO - Processing sequence 1L2Y.fasta, length: 20
2024-09-06 13:30:03,455 - INFO - Processing template 1le0.pdb with 12 residues.
2024-09-06 13:30:03,456 - INFO - Estimated time: 1-12 sec (0-0 min)
2024-09-06 13:30:03,457 - INFO - Processing template 1le1.pdb with 12 residues.
2024-09-06 13:30:03,457 - INFO - Estimated time: 1-12 sec (0-0 min)
2024-09-06 13:30:03,466 - INFO - Processing template 1vii.pdb with 36 residues.
2024-09-06 13:30:03,466 - INFO - Estimated time: 9-104 sec (0-2 min)
2024-09-06 13:30:03,470 - INFO - Processing template 5awl.pdb with 10 residues.
2024-09-06 13:30:03,470 - INFO - Estimated time: 1-8 sec (0-0 min)
2024-09-06 13:30:03,472 - INFO - Processing template 1l2y.pdb with 20 residues.
2024-09-06 13:30:03,472 - INFO - Estimated time: 3-32 sec (0-0 min)
2024-09-06 13:30:03,513 - INFO - Processing template 1crn.pdb with 46 residues.
2024-09-06 13:30:03,513 - INFO - Estimated time: 14-169 sec (0-3 min)
2024-09-06 13:30:14,150 - INFO - Processed template 5awl.pdb. Energy score: -525.06
2024-09-06 13:30:18,791 - INFO - Processed template 1le1.pdb. Energy score: -419.29
2024-09-06 13:30:19,292 - INFO - Processed template 1le0.pdb. Energy score: -394.24
2024-09-06 13:30:44,532 - INFO - Processed template 1l2y.pdb. Energy score: -469.83
2024-09-06 13:31:49,776 - INFO - Processed template 1vii.pdb. Energy score: -1051.97
2024-09-06 13:32:42,218 - INFO - Processed template 1crn.pdb. Energy score: -1315.65
2024-09-06 13:32:42,231 - INFO - Energy scores saved to 'results/example1_energy_scores.csv'.
```

#### Result

> `results/example1_result.csv`:
> |              | 1crn.pdb                        | 1l2y.pdb                        | 1le0.pdb                        | 1le1.pdb                        | 1vii.pdb                        | 5awl.pdb                        |
> |:--------------:|:----------------------------------:|:----------------------------------:|:----------------------------------:|:----------------------------------:|:----------------------------------:|:----------------------------------:|
> | 1L2Y.fasta    | -1315.65                          | -469.83                         | -394.24                          | -419.29                           | -1051.97                          | -525.06                          |
> | ![1l2y](doc/assets/1l2y.png) | ![1crn](doc/assets/1crn.png)   | ![1l2y](doc/assets/1l2y.png)    | ![1le0](doc/assets/1le0.png)    | ![1le1](doc/assets/1le1.png)    | ![1vii](doc/assets/1vii.png)    | ![5awl](doc/assets/5awl.png)    |





## Contact

For questions or issues, please open an issue on GitHub or contact [nadiajuckova@gmail.com](mailto:nadiajuckova@gmail.com).
