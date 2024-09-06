"""
DOPE Score Processing Module

This module provides functionality to find and interpolate DOPE scores 
for given residue pairs and distances based on a DOPE score data frame.

DOPE (Discrete Optimized Protein Energy) score is a statistical potential
used to evaluate the energy of protein structures. This module allows the 
user to find the score between two residues at a specific distance, using 
linear interpolation if the exact distance is not present in the DOPE data.

Usage:
    This module can be imported and used as part of a larger pipeline or for 
    standalone DOPE score calculations.

Example:
    from process_dope import find_dope_score
    score = find_dope_score('A', 'W', 3.5, dope_df, verbose=True)
    print(f"DOPE Score: {score}")
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadiajuckova@gmail.com"
__copyright__ = "MIT"
__date__ = "2024-09-06"
__version__ = "1.0.0"


import logging
import numpy as np
import pandas as pd

# Setup logging configuration if not already configured
if not logging.getLogger().hasHandlers():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')


def find_dope_score(res1, res2, distance, dope_df, verbose=False):
    """
    Find the DOPE score for a given pair of residues and distance.
    If the distance falls between two known distances, the score is linearly 
    interpolated between the two nearest values. If the distance is out of the 
    range, the closest available score is returned.

    Args:
        res1 (str): Residue 1 name (e.g., 'A').
        res2 (str): Residue 2 name (e.g., 'W').
        distance (float): The distance between residue1 and residue2.
        dope_df (pd.DataFrame): A DataFrame containing the DOPE scores, with 
                                columns 'res1', 'res2', and subsequent distance 
                                columns holding the DOPE scores.
        verbose (bool): If True, log detailed debug information.

    Returns:
        float: The interpolated or closest DOPE score for the given residue 
               pair and distance.
    
    Raises:
        ValueError: If no matching residue pair is found in the DataFrame.
        IndexError: If the provided distance is out of range of distances.
        Exception: For any other unexpected errors.
    """
    try:
        if verbose:
            logging.getLogger().setLevel(logging.DEBUG)
            logging.debug(f"Looking for DOPE score between {res1} and {res2} "
                          f"at distance {distance}")

        # Filter the DataFrame for the specific residue pair
        filtered_df = dope_df[
            (dope_df['res1'] == res1) & (dope_df['res2'] == res2)
        ]
        if filtered_df.empty:
            raise ValueError(f"No matching DOPE for residues {res1}-{res2}.")

        # Extract distance columns from the DataFrame
        distance_columns = dope_df.columns[2:].astype(float)

        # Check if the distance is outside the range of available columns
        if distance <= distance_columns.min():
            if verbose:
                logging.debug(f"Distance {distance} is below range, "
                              f"using minimum score")
            return filtered_df.iloc[0, 2]  # First column's score
        elif distance >= distance_columns.max():
            if verbose:
                logging.debug(f"Distance {distance} is above range, "
                              f"using maximum score")
            return filtered_df.iloc[0, -1]  # Last column's score

        # Find the two closest distances to the given distance
        lower_idx = np.searchsorted(distance_columns, distance) - 1
        upper_idx = lower_idx + 1

        lower_dist = distance_columns[lower_idx]
        upper_dist = distance_columns[upper_idx]

        if verbose:
            logging.debug(f"Interpolating between distances {lower_dist} "
                          f"and {upper_dist} for score")

        # DOPE scores at the two closest distances
        lower_score = filtered_df.iloc[0, 2 + lower_idx]
        upper_score = filtered_df.iloc[0, 2 + upper_idx]

        # Linear interpolation
        dope_score = lower_score + (upper_score - lower_score) * \
                     (distance - lower_dist) / (upper_dist - lower_dist)

        if verbose:
            logging.debug(f"Calculated DOPE score: {dope_score}")
            
        return dope_score

    except ValueError as ve:
        logging.error(f"ValueError: {ve}")
        raise
    except IndexError as ie:
        logging.error(f"IndexError: {ie} - Distance out of bounds.")
        raise
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise