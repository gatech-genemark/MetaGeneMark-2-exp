# Author: Karl Gemayel
# Created: 2020-06-21, 6:35 p.m.
import copy
import logging
from typing import *

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mg_general.labels import Labels
from mg_io.general import write_to_file

log = logging.getLogger(__name__)


def write_sequence_list_to_fasta_file(sequences, pf_sequences):
    # type: (List[Seq], str) -> None

    data = ""
    for i in range(len(sequences)):
        data += ">{}\n{}\n".format(i, sequences[i])

    write_to_file(data, pf_sequences)


def convert_multi_fasta_into_single_fasta(sequences, labels, new_seqname):
    # type: (Dict[str, Seq], Labels, str) -> [Dict[str, Seq], Labels]

    sorted_seqnames = sorted(sequences.keys())

    joined_seq = {new_seqname: ""}

    curr_offset = 0

    seqname_to_offset = dict()

    for seqname in sorted_seqnames:
        seqname_to_offset[seqname] = curr_offset
        seqname_sequence = sequences[seqname]
        joined_seq[new_seqname] += seqname_sequence
        curr_offset += len(seqname_sequence)

    # update all labels by corresponding offset
    labels = copy.deepcopy(labels)

    for l in labels:
        if l.seqname() in seqname_to_offset:
            l.coordinates().left += seqname_to_offset[l.seqname()]
            l.coordinates().right += seqname_to_offset[l.seqname()]
            l._fields["seqname"] = new_seqname

    joined_seq[new_seqname] = SeqRecord(seq=joined_seq[new_seqname].seq, id=new_seqname, name=new_seqname)

    return [joined_seq, labels]