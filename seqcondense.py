"""
    Sequence Consensus Condenser | RPINerd, 09/12/23

    Very early pseudocode!

    Take a user input of sequences and give back a minimal list of consensus
    sequences such that each are base-complete (i.e. no 'N's or gaps) and average
    out the conflicts between individuals in groups of similarity so that
    when designing targeting sequences the base mismatches are minimized but evenly
    distributed
"""

import os
import sys

import aln_tools
import config
import utils


def main(file):
    """
    Primary driving loop of the program
    Takes in the fasta file of sequences desired for consensus and directs
    all the sub functions to create said set of minimal consensus seqs
    """

    # Initial alignment for tree building and grouping of sequences
    fam_groups = aln_tools.init_align(file)

    # Associate each child in the families with their respective sequence
    family_dict = utils.associate_seqs(fam_groups, file)


if __name__ == "__main__":
    debug = True
    # Initialize logging
    utils.setup_logging(debug)

    # TODO should set up argparse for finer control over input and parameters
    input_fasta = sys.argv[1]

    # Set up folder for temporary alignment/tree files
    os.makedirs(config.TMP_DIR, exist_ok=True)

    main(input_fasta)

    # Clean up temporary folder, unless specified not to via debug
    if not debug:
        os.remove(config.TMP_DIR)
