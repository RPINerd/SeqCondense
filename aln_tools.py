"""
    Alignment focused tools for the pipeline
"""

import itertools
import logging
import shlex
import subprocess

from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline as clustO_cline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

import config


def align(seq_file) -> str:
    """A simple wrapper for performing a default Clustal Omega"""

    outfile_name = f"{config.TMP_DIR}{config.SEED_COUNT}.aln"
    cmd = clustO_cline(infile=seq_file, outfile=outfile_name, outfmt="phy", threads=8, force=True)
    call = shlex.split(str(cmd))
    logging.debug(f"clustal call: {call}")
    subprocess.check_call(call)
    config.SEED_COUNT += 1
    return outfile_name


def is_semipreterminal(clade):
    return any(child.is_terminal() for child in clade)


def generate_family_groups(tree) -> list:
    """
    Get the non-terminal clades, and for each one of them, check if any
    of the children are terminals. If so, add them to a family. The end
    result are a list of clades which are grouped in most-similar families.
    """

    nonTerminalClades = tree.get_nonterminals()
    groups = []
    for clade in nonTerminalClades:
        # If any of the children are terminals, add them to a family
        if is_semipreterminal(clade):
            family = []
            for child in clade:
                if child.is_terminal():
                    family.append(child)
            # Once all terminal children are gathered, add the family to the list
            # of groups if the family list is not empty
            groups.append(family) if family else None

    for fam in groups:
        logging.debug(f"Family: {[x.name for x in fam]}")

    return groups


def terminal_neighbor_dists(self):
    """Return a list of distances between adjacent terminals."""

    def generate_pairs(self):
        pairs = itertools.tee(self)
        next(pairs[1])
        return zip(pairs[0], pairs[1])

    return [self.distance(*i) for i in generate_pairs(self.find_clades(terminal=True))]


def init_align(ref_file):
    """
    Calls align to generate an MSA with all initial sequences then generates
    a pylogenetic tree to group sequences into most similar
    """

    # Create alignment and parse it into a AlignIO object
    alignment = align(ref_file)
    bioalign_obj = AlignIO.read(alignment, "phylip")

    # Create a distance matrix and tree from the MSA
    distMatrix = DistanceCalculator("identity").get_distance(bioalign_obj)
    tree = DistanceTreeConstructor().nj(distMatrix)
    logging.debug(f"Tree Full:\n{tree}")

    # Translate the tree into a list of similar-grouped families
    groups = generate_family_groups(tree)

    terminal_dists = terminal_neighbor_dists(tree)
    print(terminal_dists)

    return groups
