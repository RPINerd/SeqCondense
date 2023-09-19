"""
    Utility functions for the SeqCondense pipeline
"""

import logging

from Bio import SeqIO


def setup_logging(verbose) -> None:
    if verbose:
        logging.basicConfig(
            filename="scc.log",
            filemode="a",
            format="%(asctime)s - %(levelname)s - %(message)s",
            encoding="utf-8",
            datefmt="%H:%M:%S",
            level=logging.DEBUG,
        )
    else:
        logging.basicConfig(
            format="%(asctime)s - %(message)s",
            encoding="utf-8",
            datefmt="%M:%S",
            level=logging.INFO,
        )


def associate_seqs(fam_groups, file):
    """
    Associate each child in the families with their respective sequence
    """

    # Parse the provided file for sequences and store in a dictionary
    records = {}
    try:
        for record in SeqIO.parse(file, "fasta"):
            records[record.ID] = record.seq
    except Exception as e:
        print(f"Error on record parsing: {e}")
        exit

    # Associate each child in the families with their respective sequence
    associated_fam_groups = []
    for fam in fam_groups:
        new_fam = {}
        for child in fam:
            child.seq = next(record.seq for record in records if record.id == child.name)
