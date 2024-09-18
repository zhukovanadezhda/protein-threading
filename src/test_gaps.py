"""
Script for running sequence-structure matching tests with various gap scores, 
calculating performance, and generating result summaries. The script performs 
energy score calculations based on gap scores, checks for correctly predicted 
and homologous structures, and logs the process.

Functions:
    - calculate_performance: Calculates performance based on predictions.
    - run_tests: Runs the sequence-structure matching program for different 
      gap scores.
    - process_results: Processes the output files and calculates performance 
      results.
    - main: Main function to coordinate the tests and result processing.
"""

import argparse
import glob
import os
import subprocess
import pandas as pd
import logging

from config import gap_scores, homolog_pairs


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def calculate_performance(df, homolog_pairs):
    """
    Calculate performance metrics based on predicted structure matches.

    Args:
        df (pd.DataFrame): DataFrame containing energy scores for 
                           sequence-structure matching.
        homolog_pairs (dict): Dictionary mapping sequences to homologous 
                              structures.

    Returns:
        tuple: Performance score (float), count of correctly guessed 
               structures (int), count of similar homologous structures (int).
    """
    total_score = 0
    correctly_guessed_count = 0
    similar_structure_count = 0
    
    # Each sequence can contribute up to 3 points
    max_score = len(df) * 3  
    
    for sequence in df.index:
        correct_structure = sequence.replace('.fasta', '.pdb')
        homolog_structures = homolog_pairs.get(sequence, [])
        sorted_scores = df.loc[sequence].sort_values()

        # Check for correct structure match
        if (sorted_scores.index[0] == correct_structure or 
            sorted_scores.index[0] == correct_structure.lower()):
            total_score += 2
            correctly_guessed_count += 2
        
        # Check for homologous structures
        if any(homolog in sorted_scores.index[:2] for homolog in homolog_structures):
            total_score += 1
            similar_structure_count += 1

    performance = round(total_score / max_score, 2)
    correctly_guessed_count = int(correctly_guessed_count / 2)
    
    return performance, correctly_guessed_count, similar_structure_count


def run_tests(program_path, output_dir):
    """
    Run the sequence-structure matching program for each gap score and save 
    the results.

    Args:
        program_path (str): Path to the program to run (e.g., main.py).
        output_dir (str): Directory where output CSV files will be saved.
    """
    os.makedirs(output_dir, exist_ok=True)
    for gap_score in gap_scores:
        output_file = os.path.join(output_dir, f'energy_scores_{gap_score}.csv')
        command = ['python', program_path, '--gap_score', str(gap_score), 
                   '--output_file', output_file]
        
        logging.info(f'Running test with gap score {gap_score}...')
        try:
            subprocess.run(command, check=True)
            logging.info(f'Successfully ran test with gap score {gap_score}, '
                         f'output saved to {output_file}')
        except subprocess.CalledProcessError as e:
            logging.error(f'Error running test with gap score {gap_score}: {e}')


def process_results(output_dir, homolog_pairs, result_file):
    """
    Process the output CSV files and calculate the performance metrics for 
    each gap score.

    Args:
        output_dir (str): Directory containing the output CSV files.
        homolog_pairs (dict): Dictionary mapping sequences to homologous 
                              structures.
        result_file (str): Path to save the final performance results CSV.
    """
    csv_files = glob.glob(os.path.join(output_dir, 'energy_scores_*.csv'))
    performances = []

    for csv_file in csv_files:
        gap_score = float(csv_file.split('_')[-1].replace('.csv', ''))
        df = pd.read_csv(csv_file, index_col=0)
        performance, correctly_guessed, similar_structure = calculate_performance(
            df, homolog_pairs)

        performances.append({
            'gap_score': gap_score,
            'performance': performance,
            'correctly_guessed': correctly_guessed,
            'similar_structure': similar_structure
        })

    performance_df = pd.DataFrame(performances).sort_values(by='gap_score')\
                                               .reset_index(drop=True)
    performance_df.to_csv(result_file, index=False)

    logging.info(f'Performance results saved to {result_file}')


def main():
    """
    Main function to parse arguments, run the gap score tests, and process the 
    results.
    """
    parser = argparse.ArgumentParser(
        description="Run gap score tests and generate performance results.")
    parser.add_argument('--program_path', default='src/main.py', 
                        help="Path to the main program (e.g., main.py)")
    parser.add_argument('--output_dir', default='results/gaps_test', 
                        help="Directory to save the gap test results")
    parser.add_argument('--result_file', default='results/performance.csv', 
                        help="File to save the performance results")
    
    args = parser.parse_args()

    # Run tests with different gap scores
    run_tests(args.program_path, args.output_dir)
    
    # Process the results and generate performance CSV
    process_results(args.output_dir, homolog_pairs, args.result_file)


if __name__ == '__main__':
    main()
