"""
Sequence Shuffling and Z-Score Calculation Script

This script processes protein sequences by shuffling them multiple times to 
generate energy score distributions. It calculates z-scores for each 
sequence-template pair based on these distributions using DOPE scoring. 
The script also performs normality testing on the shuffled scores if the 
number of shuffles is below a threshold.

Features:
- Load a CSV file containing sequence and template scores.
- Generate shuffled versions of each sequence.
- Calculate energy scores for both original and shuffled sequences.
- Perform Shapiro-Wilk test on shuffled energy scores if the number of shuffles
  is less than 30, logging warnings if the distribution is not normal.
- Compute z-scores for the original sequences based on the energy score 
  distribution of shuffled sequences.
- Save the calculated z-scores to a CSV file.

Usage:
    python script.py --input_csv <path_to_input_csv> \
                     --output_file <path_to_output_csv> \
                     --gap_score <gap_score> \
                     --n_shuffles <number_of_shuffles>

Arguments:
    --input_csv : Path to the input CSV file with sequence and template scores.
    --output_file : Path where the output CSV file with z-scores will be saved.
    --gap_score : The gap score to use for energy calculation.
    --n_shuffles : The number of times to shuffle each sequence.
"""

import argparse
import pandas as pd
import numpy as np
import logging
from tqdm import tqdm
import os
from scipy.stats import shapiro

from config import DOPE_URL, TEMPLATES_DIR, SEQUENCES_DIR
from load_data import load_dope, read_fasta
from main import process_sequences_and_templates


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s')


def shuffle_sequence(sequence, n_shuffles):
    """Generate multiple shuffled versions of a sequence.

    Args:
        sequence (str): The sequence to shuffle.
        n_shuffles (int): The number of shuffled sequences to generate.

    Returns:
        list: A list of shuffled sequences.
    """
    sequence_list = list(sequence)
    shuffled_sequences = []
    for _ in range(n_shuffles):
        np.random.shuffle(sequence_list)
        shuffled_sequences.append(''.join(sequence_list))
    logging.debug(f"{n_shuffles} shuffled sequences generated.")
    return shuffled_sequences


def calculate_z_scores(original_scores, shuffled_scores):
    """Calculate z-scores based on the original and shuffled scores.

    Args:
        original_scores (numpy.ndarray): The original sequence energy scores.
        shuffled_scores (numpy.ndarray): The shuffled sequences energy scores.

    Returns:
        numpy.ndarray: The z-scores for the original sequence.
    """
    mean_shuffled = np.mean(shuffled_scores, axis=0)
    std_shuffled = np.std(shuffled_scores, axis=0)
    z_scores = (original_scores - mean_shuffled) / std_shuffled
    logging.debug("Z-scores calculated.")
    return np.round(z_scores, 2)


def perform_shapiro_test(shuffled_energy_scores, seq_file):
    """Perform Shapiro-Wilk test to check for normality of the shuffled scores.

    Args:
        shuffled_energy_scores (list): List of shuffled energy scores.
        seq_file (str): The sequence file name for logging purposes.

    Returns:
        None: Logs a warning if the distribution is not normal.
    """
    shuffled_energy_scores_flat = np.array(shuffled_energy_scores).flatten()
    stat, p_value = shapiro(shuffled_energy_scores_flat)
    if p_value < 0.05:
        logging.warning(f"The distribution of shuffled energy scores for "
                        f"sequence {seq_file} is not normally distributed "
                        f"(p-value = {p_value:.4f}).")


def process_sequence(df, seq_file, original_seq, templates, df_dope, 
                     templates_dir, gap_score, n_shuffles):
    """Process a single sequence: shuffle, calculate energies, and z-scores.

    Args:
        df (pandas.DataFrame): The DataFrame containing energy scores.
        seq_file (str): The sequence file name.
        original_seq (str): The original sequence.
        templates (list): The list of templates.
        df_dope (pandas.DataFrame): The DOPE scores DataFrame.
        templates_dir (str): The templates directory path.
        gap_score (float): The gap score to use for energy calculation.
        n_shuffles (int): The number of shuffled sequences to generate.

    Returns:
        numpy.ndarray: Z-scores for the original sequence.
    """
    logging.info(f"Processing sequence {seq_file}...")

    # Generate shuffled sequences
    shuffled_sequences = shuffle_sequence(original_seq, n_shuffles)

    # Calculate energy scores for the shuffled sequences
    shuffled_energy_scores = []
    for shuffled_seq in shuffled_sequences:
        shuffled_energy = process_sequences_and_templates(
            sequences=[(seq_file, shuffled_seq)],
            templates=templates,
            df_dope=df_dope,
            templates_dir=templates_dir,
            gap_score=gap_score
        )
        shuffled_energy_scores.append(shuffled_energy)

    # Perform Shapiro-Wilk test if n_shuffles < 30
    if n_shuffles < 30:
        perform_shapiro_test(shuffled_energy_scores, seq_file)

    # Retrieve original sequence energy scores
    original_energy_scores = df.loc[seq_file].values

    # Calculate z-scores
    shuffled_energy_scores_array = np.array(shuffled_energy_scores)
    z_scores = calculate_z_scores(original_energy_scores, 
                                  shuffled_energy_scores_array)

    return z_scores


def main(input_csv, output_file, gap_score, n_shuffles):
    """Main function to shuffle sequences and calculate z-scores.

    Args:
        input_csv (str): Path to the input CSV file.
        output_file (str): Path to save the output CSV file.
        gap_score (float): Gap score to use for energy calculation.
        n_shuffles (int): Number of times to shuffle each sequence.
    """
    logging.debug("Starting the z-score calculation process.")
    
    # Check if n_shuffles is > 1 
    if n_shuffles < 1:
        logging.error("Number of shuffles must be greater than 1.")
        return
    
    # Warn if n_shuffles is < 30
    if n_shuffles < 30:
        logging.warning("Number of shuffles is less than 30. "
                        "Shapiro-Wilk test will be performed to check for "
                        "normality of the shuffled scores.")

    # Read the input CSV
    try:
        df = pd.read_csv(input_csv, index_col=0)
        logging.debug(f"Input CSV file {input_csv} loaded successfully.")
    except FileNotFoundError:
        logging.error(f"Input CSV file {input_csv} not found.")
        return

    # Extract sequences and templates
    sequences = df.index.tolist()
    templates = df.columns.tolist()

    # Load DOPE scores
    df_dope = load_dope(DOPE_URL)

    # Initialize a DataFrame to hold z-scores
    z_scores_df = pd.DataFrame(index=sequences, columns=templates)

    # Process each sequence and calculate z-scores
    for seq_file in tqdm(sequences[::-1], desc="Processing sequences"):
        fasta_path = os.path.join(SEQUENCES_DIR, seq_file)
        if not os.path.exists(fasta_path):
            logging.warning(f"FASTA file {fasta_path} not found, skipping.")
            continue

        # Load the original sequence
        original_seq = read_fasta(fasta_path)

        # Process the sequence
        z_scores = process_sequence(
            df, seq_file, original_seq, templates, df_dope, 
            TEMPLATES_DIR, gap_score, n_shuffles
        )

        # Store z-scores in the DataFrame
        z_scores_df.loc[seq_file] = z_scores

    # Save the z-scores to an output file
    z_scores_df.to_csv(output_file)
    logging.info(f"Z-scores saved to {output_file}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Shuffle sequences and calculate energy scores.'
    )
    parser.add_argument(
        '--input_csv', 
        type=str, 
        required=True, 
        help='Input CSV file with sequence and structure scores.'
    )
    parser.add_argument(
        '--output_file',
        type=str, 
        required=True, 
        help='Output file to save shuffled scores.'
    )
    parser.add_argument(
        '--gap_score', 
        type=float, 
        required=True, 
        help='Gap score to use for energy calculation.'
    )
    parser.add_argument(
        '--n_shuffles', 
        type=int, 
        default=100,
        help='Number of times to shuffle each sequence.' 
    )

    args = parser.parse_args()

    # Run the main function with parsed arguments
    main(args.input_csv, args.output_file, args.gap_score, args.n_shuffles)
