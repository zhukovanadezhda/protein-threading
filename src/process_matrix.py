"""
Matrix Initialization and Filling Module

This module provides functionalities for initializing and filling 
low-level matrices used in sequence-structure alignment and scoring. 
The matrices are filled using distance matrices and DOPE score data 
to facilitate threading algorithms for sequence-structure prediction.

Functions:
- initialize_low_level_matrices: Initializes a 4D matrix filled with NaN.
- set_boundary_conditions: Sets boundaries for a given low-level matrix.
- fill_matrix_region: Fills a region of a low-level matrix with DOPE scores.
- fill_low_level_matrices: Fills multiple low-level matrices by setting 
  boundary conditions and filling matrix regions.
- fill_region: Fills a specified region of a low-level matrix with DOPE scores.
- fill_high_level_matrix: Fills a high-lvl matrix using the low-lvl matrices.

Usage:
    Import this module and use its functions to preprocess data for 
    threading algorithms in sequence-structure prediction.

Example:
    from matrix_operations import initialize_low_level_matrices, \
                                 set_boundary_conditions, fill_matrix_region, \
                                 fill_low_level_matrices, fill_region, \
                                 fill_high_level_matrix
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadiajuckova@gmail.com"
__copyright__ = "MIT"
__date__ = "2024-09-06"
__version__ = "1.0.0"


import logging
import numpy as np
import pandas as pd
from process_dope import find_dope_score

# Setup logging configuration if not already configured
if not logging.getLogger().hasHandlers():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')


def initialize_low_level_matrices(n: int, m: int) -> np.ndarray:
    """
    Initialize a low-level matrix filled with NaN values.
    
    Args:
        n (int): Number of rows for the matrix.
        m (int): Number of columns for the matrix.
    
    Returns:
        np.ndarray: 4D matrix initialized with NaN values.
    """
    logging.debug(f"Initializing low-level matrices with dimensions ({n}, {m})")
    return np.full((n, m, n, m), np.nan, dtype=float)


def set_boundary_conditions(
    low_level_matrix: np.ndarray, 
    i: int, 
    j: int
    ) -> None:
    """
    Set the boundary conditions for a low-level matrix.
    
    Args:
        low_level_matrix (np.ndarray): The matrix to set boundary conditions.
        i (int): Index for the row boundary.
        j (int): Index for the column boundary.
    """
    logging.debug(f"Setting boundary conditions for matrix "
                  f"with boundaries ({i}, {j})")
    low_level_matrix[0, :j+1] = 0  # Set first row up to j
    low_level_matrix[:i+1, 0] = 0  # Set first column up to i


def fill_matrix_region(
    low_level_matrix: np.ndarray, 
    k_range: range, 
    l_range: range, 
    sequence_j: str, 
    dist_matrix_i: np.ndarray, 
    df_dope: pd.DataFrame, 
    sequence: list, 
    skip: set = None
) -> None:
    """
    Fill a region of a low-level matrix dynamically calculating the dope score.
    
    Args:
        low_level_matrix (np.ndarray): The matrix to fill.
        k_range (range): The range for the k indices.
        l_range (range): The range for the l indices.
        sequence_j (str): Residue at position j in the sequence.
        dist_matrix_i (np.ndarray): The i-th row of the distance matrix.
        df_dope (pd.DataFrame): DataFrame containing DOPE scores.
        sequence (list): Sequence of residues.
        skip (set, optional): Set of (k, l) tuples to skip. Defaults to None.
    """
    logging.debug(f"Filling matrix region for sequence residue '{sequence_j}'")
    for k in k_range:
        for l in l_range:
            if skip and (k, l) in skip:
                continue
            # Calculate DOPE score and fill the matrix
            dope = find_dope_score(
                res1=sequence_j,
                res2=sequence[l],
                distance=dist_matrix_i[k],
                dope_df=df_dope
            )
            # Find the minimum value from neighboring cells
            min_choice = np.nanmin([
                low_level_matrix[k-1, l],  # From above
                low_level_matrix[k, l-1],  # From left
                low_level_matrix[k-1, l-1] # From diagonal
            ])
            low_level_matrix[k, l] = min_choice + dope
            logging.debug(f"Set low_level_matrix[{k}, {l}] "
                          f"= {low_level_matrix[k, l]}")


def fill_low_level_matrices(
    n: int, 
    m: int, 
    sequence: list, 
    dist_matrix: np.ndarray, 
    df_dope: pd.DataFrame
) -> np.ndarray:
    """
    Fill low-level matrices by setting boundary conditions and filling regions.
    
    Args:
        n (int): Number of rows for the matrix.
        m (int): Number of columns for the matrix.
        sequence (list): List of sequence residues.
        dist_matrix (np.ndarray): Distance matrix for the sequence.
        df_dope (pd.DataFrame): DataFrame containing DOPE scores.
    
    Returns:
        np.ndarray: 4D low-level matrices filled with computed values.
    """
    logging.debug(f"Filling low-level matrices with dimensions ({n}, {m})")
    low_level_matrices = initialize_low_level_matrices(n, m)

    for i in range(n):
        for j in range(m):
            low_level_matrix = low_level_matrices[i, j]
            
            # Set initial boundary conditions
            set_boundary_conditions(low_level_matrix, i, j)

            # Fill the left-top region from (1,1) to (i,j)
            fill_matrix_region(
                low_level_matrix=low_level_matrix,
                k_range=range(1, i+1),
                l_range=range(1, j+1),
                sequence_j=sequence[j],
                dist_matrix_i=dist_matrix[i],
                df_dope=df_dope,
                sequence=sequence
            )

            # Fill the right-bottom region from (i,j) to (n-1,m-1), skip (i,j)
            fill_matrix_region(
                low_level_matrix=low_level_matrix,
                k_range=range(i, n),
                l_range=range(j, m),
                sequence_j=sequence[j],
                dist_matrix_i=dist_matrix[i],
                df_dope=df_dope,
                sequence=sequence,
                skip={(i, j)}
            )

    return low_level_matrices


def fill_region(
    low_level_matrix: np.ndarray, 
    k_range: range, 
    l_range: range, 
    sequence_j: str, 
    dist_matrix_i: np.ndarray, 
    df_dope: pd.DataFrame, 
    sequence: list, 
    skip: set = None
    ) -> None:
    """
    Fill a specified region of the low_level_matrix with DOPE scores.

    Args:
        low_level_matrix (np.ndarray): The matrix to fill.
        k_range (range): Range for the k indices.
        l_range (range): Range for the l indices.
        sequence_j (str): The residue at position j in the sequence.
        dist_matrix_i (np.ndarray): The i-th row of the distance matrix.
        df_dope (pd.DataFrame): DataFrame containing DOPE scores.
        sequence (list): The sequence of residues.
        skip (set, optional): Set of (k, l) tuples to skip. Defaults to None.
    """
    logging.debug(f"Filling region in low-level matrix "
                  f"for sequence residue '{sequence_j}'")
    for k in k_range:
        for l in l_range:
            if skip and (k, l) in skip:
                continue
            
            dope = find_dope_score(
                res1=sequence_j,
                res2=sequence[l],
                distance=dist_matrix_i[k],
                dope_df=df_dope
            )
            
            # Get the minimum choice from neighboring cells
            min_choice = np.nanmin([
                low_level_matrix[k-1, l],    # Top
                low_level_matrix[k, l-1],    # Left
                low_level_matrix[k-1, l-1]   # Diagonal
            ])
            
            # Fill the matrix at position (k, l)
            low_level_matrix[k, l] = min_choice + dope
            logging.debug(f"Set low_level_matrix[{k}, {l}] "
                          f"= {low_level_matrix[k, l]}")


def fill_high_level_matrix(low_level_matrices: np.ndarray) -> np.ndarray:
    """
    Fill the high-level matrix using the low-level matrices.

    Args:
        low_level_matrices (np.ndarray): The low-level matrices.

    Returns:
        np.ndarray: The filled high-level matrix.
    """
    logging.debug("Filling high-level matrix using low-level matrices")
    n, m = low_level_matrices.shape[:2]
    
    # Initialize high_level_matrix with 0 values
    high_level_matrix = np.full((n, m), 0, dtype=float)
    
    for i in range(1, n):
        for j in range(1, m):
            min_choice = min(high_level_matrix[i-1, j], 
                             high_level_matrix[i, j-1], 
                             high_level_matrix[i-1, j-1])
            high_level_matrix[i, j] = min_choice + \
                                      low_level_matrices[i, j][-1, -1]
            logging.debug(f"Set high_level_matrix[{i}, {j}] "
                          f"= {high_level_matrix[i, j]}")
    
    return high_level_matrix
