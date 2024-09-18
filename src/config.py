"""
Configuration and constants for the protein threading algorithm.
"""

# DOPE score URL
DOPE_URL = "https://www.dsimb.inserm.fr/~gelly/doc/dope.par"

# Directory paths
TEMPLATES_DIR = 'data/example1/structures/'
SEQUENCES_DIR = 'data/example1/sequences/'

# Gap scores to test
gap_scores = [0, 0.1, 0.2, 0.3, 0.5, 1, 2, 5]

# Proteins with similar structures
homolog_pairs = {
    '5AWL.fasta': ['1l2y.pdb', '1vii.pdb', '1crn.pdb'],
    '1VII.fasta': ['1l2y.pdb', '1crn.pdb'],
    '1L2Y.fasta': ['1vii.pdb', '1crn.pdb'],
    '1CRN.fasta': ['1l2y.pdb','1vii.pdb', '1crn.pdb', '1le0.pdb', '1le1.pdb']
}
