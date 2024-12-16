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

from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo

import utils


def main(file):
    """
    Primary driving loop of the program
    Takes in the fasta file of sequences desired for consensus and directs
    all the sub functions to create said set of minimal consensus seqs
    """

    # Parse the provided file for sequences and store in list
    records = []
    try:
        for record in SeqIO.parse(file, "fasta"):
            records.append(record)
    except Exception as e:
        print(f"Error on record parsing: {e}")
        exit

    # Initial alignment for tree building and grouping of sequences
    utils.init_align(file)

    for sequence in records:
        reference_list = records.remove(sequence)
        for ref_sequence in reference_list:
            alignment = AlignIO.align()
            summary_align = AlignInfo.SummaryInfo(alignment)
            summary_align.dumb_consensus(float(sys.argv[2]))

            # if 'align better':
            #     pass
            #     # save to best alignment
            # else next fasta

        # save best match to new struct


if __name__ == "__main__":
    debug = True
    # Initialize logging
    utils.logging(debug)

    # TODO should set up argparse for finer control over input and parameters
    input_fasta = sys.argv[1]

    # Set up folder for temporary alignment/tree files
    os.makedirs(utils.TMP_DIR, exist_ok=True)

    main(input_fasta)

    # Clean up temporary folder, unless specified not to via debug
    if not debug:
        os.remove(utils.TMP_DIR)
