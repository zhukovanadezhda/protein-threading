"""
Sequence and Template Processing Script

This script processes sequences and templates to calculate energy scores
using DOPE scoring. It loads sequence and template files, computes distance
matrices, and fills low-level and high-level matrices to derive energy scores.

Features:
- Load sequence and template files.
- Compute distance matrices from template coordinates.
- Fill low-level matrices with DOPE scores.
- Calculate energy scores for sequence-template pairs.
- Save the results to a CSV file.
"""

__authors__ = "Nadezhda Zhukova"
__contact__ = "nadiajuckova@gmail.com"
__copyright__ = "MIT"
__date__ = "2024-09-06"
__version__ = "1.2.0"

import argparse
from joblib import Parallel, delayed
import os
import pandas as pd
from multiprocessing import cpu_count
import logging
from config import DOPE_URL, SEQUENCES_DIR, TEMPLATES_DIR
from load_data import (load_dope, read_fasta, pdb_to_c_alpha_coordinates,
                       coordinates_to_distance_matrix)
from process_matrix import fill_low_level_matrices, fill_high_level_matrix

# Setup logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)


def load_templates(templates_dir: str) -> list:
    """
    Load template filenames from the specified directory.

    Args:
        templates_dir (str): Directory containing template files.

    Returns:
        list: List of template filenames.
    """
    return os.listdir(templates_dir)


def load_sequences(sequences_dir: str) -> list:
    """
    Load sequence filenames from the specified directory.

    Args:
        sequences_dir (str): Directory containing sequence files.

    Returns:
        list: List of sequence filenames.
    """
    return os.listdir(sequences_dir)


def process_template(
    template: str, templates_dir: str, sequence: str, df_dope: pd.DataFrame,
    jobs: int,
    verbose: bool = False
) -> float:
    """
    Process a single template by calculating the energy score for a sequence.

    Args:
        template (str): The template file name.
        templates_dir (str): Directory where templates are stored.
        sequence (str): Sequence of residues.
        df_dope (pd.DataFrame): DataFrame containing DOPE scores.
        verbose (bool): If True, enables verbose output.

    Returns:
        float: Computed energy score for the template.

    Raises:
        Exception: If processing the template fails.
    """
    try:
        pdb_file = os.path.join(templates_dir, template)
        coords = pdb_to_c_alpha_coordinates(pdb_file)
        dist_matrix = coordinates_to_distance_matrix(coords)
        n = coords.shape[0]
        m = len(sequence)
        complexity = n**2 * m**2
        logging.info(f"Processing template {template} with {n} residues.")
        logging.info(
            f"Estimated time: "
            f"{round(complexity * 2e-4 / jobs)}-{round(complexity * 2e-4)} sec "
            f"({round(complexity * 3e-6 / jobs)}-{round(complexity * 3e-6)} min)"
        )

        low_level_matrices = fill_low_level_matrices(
            n, m, sequence, dist_matrix, df_dope
        )
        energy_score = fill_high_level_matrix(low_level_matrices)[-1, -1]
        energy_score = round(energy_score, 2)
        logging.info(f"Processed template {template}. "
                     f"Energy score: {energy_score}")
        return energy_score
    except Exception as e:
        logging.error(f"Error processing template {template}: {e}")
        return float('nan')


def process_template_wrapper(
    template: str, sequence: str, templates_dir: str, df_dope: pd.DataFrame,
    jobs: int, verbose: bool
) -> float:
    """
    Wrapper function for process_template to use with multiprocessing.
    """
    return process_template(
        template, templates_dir, sequence, df_dope, jobs, verbose
    )


def process_sequences_and_templates(
    sequences: list, templates: list, df_dope: pd.DataFrame, 
    templates_dir: str, verbose: bool = False, 
    dry_run: bool = False, jobs: int = cpu_count()
) -> pd.DataFrame:
    """
    Process all sequences from the list and compare them
    with all templates, using parallel processing with joblib.

    Args:
        sequences (list): List of sequences to process.
        templates (list): List of templates to process.
        df_dope (pd.DataFrame): DataFrame containing DOPE scores.
        templates_dir (str): Directory containing template files.
        verbose (bool): If True, enables verbose output.
        dry_run (bool): If True, only log actions without processing.
        jobs (int): Number of parallel jobs to use for processing.

    Returns:
        pd.DataFrame: A DataFrame containing energy scores for all
                      sequence-template pairs.
    """
    energy_scores = {}

    for sequence_file, sequence in sequences:
        logging.info(
            f"Processing sequence {sequence_file}, length: {len(sequence)}"
        )

        energy_scores[sequence_file] = {}

        # Dry run - Log what would be processed
        if dry_run:
            for template in templates:
                logging.info(f"Dry run: would process {sequence_file} "
                             f"with template {template}")
            continue

        # Parallelize the template processing
        results = Parallel(n_jobs=jobs)(
            delayed(process_template_wrapper)(
                template, sequence, templates_dir, df_dope, jobs, verbose
            ) for template in templates
        )

        # Store the energy scores for each template
        for template, energy_score in zip(templates, results):
            energy_scores[sequence_file][template] = energy_score

    return pd.DataFrame(energy_scores).T


def save_energy_scores(
    energy_scores_df: pd.DataFrame, 
    filename: str = 'energy_scores.csv'
) -> None:
    """
    Save the energy scores DataFrame to a CSV file.

    Args:
        energy_scores_df (pd.DataFrame): DataFrame containing energy scores.
        filename (str): Path to the output CSV file.
    """
    energy_scores_df.to_csv(filename)
    logging.info(f"Energy scores saved to '{filename}'.")


def main() -> None:
    """
    Main function to execute the script. Loads DOPE score data, processes
    sequences and templates to compute energy scores, and saves the results
    to a CSV file.

    Command-line Arguments:
        --sequences (optional): A comma-separated list of sequence filenames.
        --templates (optional): A comma-separated list of template filenames.
        --output_file (optional): Path for the output CSV file.
        --jobs (optional): Number of parallel jobs to use.
        --dry_run (optional): If set, log actions without processing.
        --verbose (optional): Enables verbose output (debug-level logging).
    """
    parser = argparse.ArgumentParser(
        description="Process sequences and templates to calculate energy."
    )
    parser.add_argument(
        '--sequences', 
        type=str, 
        help='Comma-separated list of sequence filenames'
        )
    parser.add_argument(
        '--templates', 
        type=str, 
        help='Comma-separated list of template filenames'
        )
    parser.add_argument(
        '--output_file',
        type=str, 
        default='results/energy_scores.csv', 
        help='Name of the output CSV file'
        )
    parser.add_argument(
        '--jobs', 
        type=int, 
        default=cpu_count(), 
        help='Number of parallel jobs to run, default is all cores'
        )
    parser.add_argument(
        '--dry_run', 
        action='store_true', 
        help='If set, only log actions without processing'
        )
    parser.add_argument(
        '--verbose', 
        action='store_true', 
        help='Enable verbose output'
        )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logging.info("Loading DOPE score data...")
    df_dope = load_dope(DOPE_URL)

    # Load all sequences or only the specified ones
    if args.sequences:
        sequence_files = args.sequences.split(',')
        sequences = [
            (file, 
             read_fasta(os.path.join(SEQUENCES_DIR, file))) 
            for file in sequence_files
            ]
    else:
        sequence_files = load_sequences(SEQUENCES_DIR)
        sequences = [
            (file, read_fasta(os.path.join(SEQUENCES_DIR, file))) 
            for file in sequence_files
            ]

    # Load all templates or only the specified ones
    if args.templates:
        templates = args.templates.split(',')
    else:
        templates = load_templates(TEMPLATES_DIR)

    logging.info("Processing sequences and templates...")
    energy_scores_df = process_sequences_and_templates(
        sequences, 
        templates, 
        df_dope, 
        TEMPLATES_DIR, 
        verbose=args.verbose, 
        dry_run=args.dry_run, 
        jobs=args.jobs
    )

    if not args.dry_run:
        save_energy_scores(energy_scores_df, args.output_file)


if __name__ == "__main__":
    main()
