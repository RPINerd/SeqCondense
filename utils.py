"""
    Utility functions for the SeqCondense pipeline
"""

import logging
import subprocess

from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline as clustO_cline
from Bio.Phylo.TreeConstruction import DistanceCalculator

SEED_COUNT = 1
TMP_DIR = "./tmp/"


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


def align(seq_file) -> str:
    """A simple wrapper for performing a default Clustal Omega"""

    outfile_name = f"{TMP_DIR}{SEED_COUNT}.aln"
    cmd = clustO_cline(infile=seq_file, outfile=outfile_name, outfmt="phy", threads=8)
    logging.debug(f"clustal call: {cmd}")
    subprocess.check_call(cmd)
    return outfile_name


def init_align(ref_file) -> None:
    """
    Calls align to generate an MSA with all initial sequences then generates
    a pylogenetic tree to group sequences into most similar
    """

    alignment = align(ref_file)
    bioalign_obj = AlignIO.read(alignment, "phylip")
    distMatrix = DistanceCalculator("identity").get_distance(bioalign_obj)

    print(distMatrix)

    return
