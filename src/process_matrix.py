"""
Matrix Initialization and Filling Module

This module provides functionalities for initializing and filling 
low-level matrices used in sequence-structure alignment and scoring. 
The matrices are filled using distance matrices and DOPE score data 
to enable threading algorithms for sequence-structure prediction.

Functions:
- initialize_low_level_matrices: Initializes a 4D matrix filled with NaN.
- set_boundary_conditions: Sets boundaries for a given low-level matrix.
- fill_matrix_region: Fills a region of a low-level matrix with DOPE scores.
- fill_low_level_matrices: Fills multiple low-level matrices by setting 
  boundary conditions and filling matrix regions.
- fill_high_level_matrix: Fills a high-lvl matrix using the low-lvl matrices.

Usage:
    Import this module and use its functions to preprocess data for 
    threading algorithms in sequence-structure prediction.

Example:
    from process_matrix import initialize_low_level_matrices, \
                               set_boundary_conditions, fill_matrix_region, \
                               fill_low_level_matrices, \
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
from joblib import Parallel, delayed
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
    logging.debug(f"Initializing low-lvl matrices with dimensions ({n}, {m})")
    return np.full((n, m, n, m), np.nan, dtype=float)


def set_boundary_conditions(low_level_matrix: np.ndarray, i: int, j: int,
                            sequence: str, dist_matrix: np.ndarray,
                            df_dope: pd.DataFrame, gap_score: float) -> None:

    if i == 0 and j == 0:
        low_level_matrix[0, 0] = 0
    else:
        dope_score = find_dope_score(
            res1=sequence[0], 
            res2=sequence[j],  
            distance=dist_matrix[0, i],
            dope_df=df_dope
        )
        low_level_matrix[0, 0] = round(dope_score, 2)
    
    # Initialize first row (i=0) using values from the left
    for j in range(1, j+1):
        low_level_matrix[0, j] = round(low_level_matrix[0, j-1] + gap_score, 2)

    # Initialize first column (j=0) using values from above
    for i in range(1, i+1):
        low_level_matrix[i, 0] = round(low_level_matrix[i-1, 0] + gap_score, 2)



def fill_matrix_region(low_level_matrix: np.ndarray, k_range: range, 
                       l_range: range, sequence_j: str, 
                       dist_matrix_i: np.ndarray, df_dope: pd.DataFrame, 
                       sequence: list, gap_score: float,
                       skip: set = None) -> None:
    """
    Fill a region of a low-level matrix dynamically calculating the DOPE score.
    
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
            
            # Skip if the cell is in the skip set
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
            low_level_matrix[k, l] = np.nanmin([
                round(low_level_matrix[k-1, l] + gap_score, 2), # From above
                round(low_level_matrix[k, l-1] + gap_score, 2), # From left
                round(low_level_matrix[k-1, l-1] + dope, 2)     # From diagonal
            ])
            logging.debug(f"Set low_level_matrix[{k}, {l}] "
                          f"= {low_level_matrix[k, l]}")


def fill_low_level_matrices(n: int, m: int, sequence: list, 
                            dist_matrix: np.ndarray, gap_score: float,
                            df_dope: pd.DataFrame) -> np.ndarray:
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
            set_boundary_conditions(low_level_matrix=low_level_matrix,
                                    i=i, j=j,
                                    sequence=sequence,
                                    dist_matrix=dist_matrix,
                                    df_dope=df_dope,
                                    gap_score=gap_score)

            # Fill the left-top region from (1,1) to (i,j)
            fill_matrix_region(
                low_level_matrix=low_level_matrix,
                k_range=range(1, i+1),
                l_range=range(1, j+1),
                sequence_j=sequence[j],
                dist_matrix_i=dist_matrix[i],
                df_dope=df_dope,
                gap_score=gap_score,
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
                gap_score=gap_score,
                sequence=sequence,
                skip={(i, j)}
            )

    return low_level_matrices


def fill_high_level_matrix(low_level_matrices: np.ndarray, 
                           gap_score: float, sequence: str, 
                           print_alignments: bool) -> np.ndarray:
    """
    Fill the high-level matrix using the low-level matrices.

    Args:
        low_level_matrices (np.ndarray): The low-level matrices.
        gap_score (float): The gap score to be used.

    Returns:
        np.ndarray: The filled high-level matrix.
    """
    logging.debug("Filling high-level matrix using low-level matrices")
    n, m = low_level_matrices.shape[:2]
    
    # Initialize high_level_matrix
    high_level_matrix = np.full((n, m), np.nan, dtype=float)
    
    high_level_matrix[0, 0] = low_level_matrices[0, 0][-1, -1]
    
    for i in range(1, n):
        high_level_matrix[i, 0] = round(
            high_level_matrix[i-1, 0] + gap_score, 2
            )
    for j in range(1, m):
        high_level_matrix[0, j] = round(
            high_level_matrix[0, j-1] + gap_score, 2
            )

    for i in range(1, n):
        for j in range(1, m):
            high_level_matrix[i, j] = np.nanmin([
                round(high_level_matrix[i-1, j] + gap_score, 2), 
                round(high_level_matrix[i, j-1] + gap_score, 2), 
                round(
                    high_level_matrix[i-1, j-1] + low_level_matrices[i, j][-1, -1], 2
                    )
                ])
            
    # Print alignment
    if print_alignments:
        traceback_alignment(low_level_matrices, high_level_matrix, 
                            gap_score, sequence, n, m)
    
    return high_level_matrix


def traceback_alignment(low_level_matrices: np.ndarray,
                        high_level_matrix: np.ndarray, gap_score: float, 
                        sequence: str, n: int, m: int, 
                        max_line_length: int = 60) -> tuple:
    """Traceback the alignment and print the aligned sequences.

    Args:
        low_level_matrices (np.ndarray): Matrix of low-level matrices.
        high_level_matrix (np.ndarray): High-level matrix.
        gap_score (float): Gap score used in the alignment.
        sequence (str): Sequence to align.
        n (int): Length of the structure.
        m (int): Length of the sequence.
        max_line_length (int): Maximum character length per line (default: 60).

    Returns:
        tuple: Strings with the aligned indices, residues, and connectors.
    """
    aligned_indices = []
    aligned_residues = []
    alignment_connectors = []
    
    # Start from the bottom-right corner of the high-level matrix
    i, j = n - 1, m - 1
    
    while (i > 0) and (j > 0):
        current_value = high_level_matrix[i, j]
        low_level_value = low_level_matrices[i, j][-1, -1]
        
        vertical_move = high_level_matrix[i-1, j] + gap_score
        horizontal_move = high_level_matrix[i, j-1] + gap_score
        diagonal_move = high_level_matrix[i-1, j-1] + low_level_value
        
        if np.isclose(current_value, round(vertical_move, 2)):
            # Move up (gap in sequence[j])
            nb_spaces = 3 - len(str(i))
            aligned_indices.append(f'{" " * nb_spaces}{i}')
            aligned_residues.append('  -')
            alignment_connectors.append('   ')
            i -= 1
        
        elif np.isclose(current_value, round(horizontal_move, 2)):
            # Move left (gap in sequence[i])
            aligned_indices.append('  -')
            aligned_residues.append(f'  {sequence[j-1]}')
            alignment_connectors.append('   ')
            j -= 1
        
        elif np.isclose(current_value, round(diagonal_move, 2)):
            # Move diagonally (match/mismatch)
            nb_spaces = 3 - len(str(i))
            aligned_indices.append(f'{" " * nb_spaces}{i}')
            aligned_residues.append(f'  {sequence[j-1]}')
            alignment_connectors.append('  |')
            i -= 1
            j -= 1

    # Add the starting point (i=0, j=0)
    aligned_indices.append('  0')
    aligned_residues.append(f'  {sequence[0]}')
    alignment_connectors.append('  |')
    
    # Reverse the lists since we built them backward
    aligned_indices.reverse()
    aligned_residues.reverse()
    alignment_connectors.reverse()
    
    # Convert lists to strings
    aligned_indices_str = ''.join(aligned_indices)
    alignment_connectors_str = ''.join(alignment_connectors)
    aligned_residues_str = ''.join(aligned_residues)
    
    # Add line breaks every 60 characters
    aligned_indices_lines = [aligned_indices_str[i:i+max_line_length] 
                            for i in range(
                                0, 
                                len(aligned_indices_str), 
                                max_line_length
                                )
                            ]
    alignment_connectors_lines = [alignment_connectors_str[i:i+max_line_length] 
                                          for i in range(
                                              0, 
                                              len(alignment_connectors_str), 
                                              max_line_length
                                              )
                                 ]
    aligned_residues_lines = [aligned_residues_str[i:i+max_line_length] 
                                      for i in range(
                                          0, 
                                          len(aligned_residues_str), 
                                          max_line_length
                                          )
                             ]
    
    alignment_str = ''
    
    # Combine them in alternating order (indices -> connectors -> residues)
    for idx_line, con_line, res_line in zip(
                                            aligned_indices_lines, 
                                            alignment_connectors_lines, 
                                            aligned_residues_lines
                                            ):
        alignment_str += idx_line + '\n' + con_line + '\n' + res_line + '\n'
        
    print(alignment_str[:-1])
    
    return aligned_indices_str, alignment_connectors_str, aligned_residues_str

