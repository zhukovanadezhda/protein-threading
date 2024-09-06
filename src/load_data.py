"""
Sequence and Structure Data Processing Module

This module provides functionalities to read FASTA sequences, extract
C-alpha coordinates from PDB files, compute distance matrices, and
load DOPE score data. These functionalities are essential for
threading algorithms used in sequence-structure prediction.

Functions:
- read_fasta: Reads a FASTA file and returns the sequence as a string.
- pdb_to_c_alpha_coordinates: Converts a PDB file to a numpy array of 
  C-alpha coordinates.
- coordinates_to_distance_matrix: Converts a numpy array of coordinates 
  to a distance matrix.
- load_dope: Loads and prepares the DOPE score data from the given URL.

Usage:
    This module can be imported and used to preprocess data for threading
    algorithms in sequence-structure prediction. 

Example:
    from load_data import read_fasta, pdb_to_c_alpha_coordinates, \
                             coordinates_to_distance_matrix, load_dope
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadiajuckova@gmail.com"
__copyright__ = "MIT"
__date__ = "2024-09-06"
__version__ = "1.0.0"


import logging
import os
from Bio.SeqUtils import IUPACData
import numpy as np
import pandas as pd

# Setup logging configuration if not already configured
if not logging.getLogger().hasHandlers():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')


def read_fasta(fasta_file: str) -> str:
    """
    Read a FASTA file and return the sequence as a string.

    Args:
        fasta_file (str): The path to the FASTA file.

    Returns:
        str: The sequence string.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
        ValueError: If the sequence contains invalid residues.
    """
    logging.debug(f"Reading FASTA file from '{fasta_file}'")
    try:
        with open(fasta_file, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found.")
    
    sequence = ''.join(line.strip().upper() 
                       for line in lines 
                       if not line.startswith('>'))
    
    # Validate the sequence to ensure it only contains known residues
    valid_residues = set(IUPACData.protein_letters)
    if not set(sequence).issubset(valid_residues):
        invalid_residues = set(sequence) - valid_residues
        raise ValueError(
            f"Sequence contains invalid residues: "
            f"{', '.join(invalid_residues)}")
    
    logging.debug(f"FASTA sequence length: {len(sequence)}")
    return sequence


def pdb_to_c_alpha_coordinates(pdb_file: str) -> np.ndarray:
    """Convert a PDB file to a numpy array of C-alpha coordinates.

    Args:
        pdb_file (str): The path to the PDB file.

    Returns:
        np.ndarray: A numpy array of C-alpha coordinates with shape (n, 3),
                    where n is the number of residues in the protein.
    """
    logging.debug(f"Extracting C-alpha coordinates from PDB file '{pdb_file}'")
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    coords = []
    for line in lines:
        if line.startswith('ATOM') and 'CA' in line:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
    
    coords = np.array(coords)
    logging.debug(f"Extracted {coords.shape[0]} C-alpha coordinates")
    return coords


def coordinates_to_distance_matrix(coords: np.ndarray) -> np.ndarray:
    """Convert a numpy array of coordinates to a distance matrix.

    Args:
        coords (np.ndarray): A numpy array of coordinates with shape (n, 3),
                             where n is the number of residues.

    Returns:
        np.ndarray: A distance matrix with shape (n, n).
    """
    logging.debug("Calculating distance matrix")
    n = coords.shape[0]
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            dist_matrix[i, j] = np.linalg.norm(coords[i] - coords[j])
            dist_matrix[j, i] = dist_matrix[i, j]
    
    logging.debug(f"Distance matrix shape: {dist_matrix.shape}")
    return dist_matrix


def load_dope(url: str) -> pd.DataFrame:
    """
    Load and prepare the DOPE score data from the given URL.
    
    Args:
        url (str): The URL to the DOPE score data file.
    
    Returns:
        pd.DataFrame: Processed DataFrame with only CA-CA interactions and 
                      converted residue names to single-letter format.
    """
    logging.debug(f"Loading DOPE score data from '{url}'")
    # Define column names (res1, atom1, res2, atom2, and distance columns)
    residue_columns = ["res1", "atom1", "res2", "atom2"]
    distance_columns = [str(i) for i in np.arange(0.25, 15.25, 0.5)]
    columns = residue_columns + distance_columns
    
    # Load the DOPE data from the URL
    df_dope = pd.read_csv(url, sep=r'\s+', header=None, names=columns)
    
    # Filter for CA-CA interactions and drop the atom columns
    df_dope = df_dope[(df_dope["atom1"] == "CA") & (df_dope["atom2"] == "CA")]
    df_dope = df_dope.drop(columns=["atom1", "atom2"]).reset_index(drop=True)
    
    # Convert residue names from 3-letter to 1-letter format
    df_dope["res1"] = df_dope["res1"].apply(
        lambda x: IUPACData.protein_letters_3to1[x.lower().capitalize()]
    )
    df_dope["res2"] = df_dope["res2"].apply(
        lambda x: IUPACData.protein_letters_3to1[x.lower().capitalize()]
    )
    
    logging.debug(f"DOPE score data loaded with shape: {df_dope.shape}")
    return df_dope
