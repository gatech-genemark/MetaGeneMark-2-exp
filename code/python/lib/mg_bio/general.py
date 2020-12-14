# Author: Karl Gemayel
# Created: 2020-06-21, 5:19 p.m.

import logging
from typing import *
from Bio import SeqIO

log = logging.getLogger(__name__)


def compute_gc_from_sequences(sequences):
    # type: (Dict[str, SeqIO.SeqRecord]) -> float
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    for seq in sequences.values():
        for s in seq:
            if s.upper() in {"A", "C", "G", "T"}:
                counts[s] += 1

    total = sum(counts.values())
    count_gc = counts["G"] + counts["C"]

    if total == 0:
        return 0.0

    return 100 * count_gc / float(total)


def compute_single_gc_from_file(pf_sequences):
    # type: (str) -> float
    record_dict = SeqIO.to_dict(SeqIO.parse(pf_sequences, "fasta"))

    return compute_gc_from_sequences(record_dict)

