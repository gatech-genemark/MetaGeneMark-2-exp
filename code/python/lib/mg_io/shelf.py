# Author: Karl Gemayel
# Created: 2020-06-21, 6:35 p.m.
import copy
import logging
from typing import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mg_general.general import os_join
from mg_general.labels import Labels
from mg_io.general import write_to_file
from mg_io.labels import read_labels_from_file

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


def read_sequences_for_gi(env, gi):
    # type: (Environment, GenomeInfo) -> Dict[str, SeqRecord]
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    return SeqIO.to_dict(SeqIO.parse(pf_sequence, "fasta"))


def read_labels_for_gi(env, gi, fn_labels="ncbi.gff"):
    # type: (Environment, GenomeInfo, str) -> Labels
    pf_labels = os_join(env["pd-data"], gi.name, fn_labels)
    return read_labels_from_file(pf_labels)