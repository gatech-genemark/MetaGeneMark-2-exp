# Author: Karl Gemayel
# Created: 2020-06-21, 6:35 p.m.

import logging
from typing import *

from Bio.Seq import Seq

from mg_io.general import write_to_file

log = logging.getLogger(__name__)


def write_sequence_list_to_fasta_file(sequences, pf_sequences):
    # type: (List[Seq], str) -> None

    data = ""
    for i in range(len(sequences)):
        data += ">{}\n{}\n".format(i, sequences[i])

    write_to_file(data, pf_sequences)
