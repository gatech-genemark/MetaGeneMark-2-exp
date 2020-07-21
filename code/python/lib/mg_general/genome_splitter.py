# Author: Karl Gemayel
# Created: 7/20/20, 2:58 PM

import logging
from typing import *

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree, Interval

from mg_general import Environment
from mg_general.general import get_value, os_join
from mg_general.labels import Labels

log = logging.getLogger(__name__)


class GenomeSplitter:

    def __init__(self, sequences, fragment_size_nt, **kwargs):
        # type: (Dict[str, SeqRecord], int, Dict[str, Any]) -> None

        self._split_genome(sequences, fragment_size_nt, **kwargs)

    def _split_genome(self, sequences, fragment_size_nt, **kwargs):
        # type: (Dict[str, SeqRecord], int, Dict[str, Any]) -> None

        self.split_sequences_ = GenomeSplitter.split_fasta_into_chunks(sequences, fragment_size_nt, **kwargs)

    def write_to_file(self, pf_output, single=True):
        # type: (str, bool) -> None

        if single:
            with open(pf_output, "w") as f:
                SeqIO.write([x[0] for x in self.split_sequences_], f, "fasta")
        else:
            counter = 0
            for s in self.split_sequences_:
                pf_output_chunk = f"{pf_output}_{counter}"
                with open(pf_output_chunk, "w") as f:
                    SeqIO.write(s[0], f, "fasta")
                    counter += 1

    @staticmethod
    def split_fasta_into_chunks(sequences, chunk_size_nt, **kwargs):
        # type: (Dict[str, SeqRecord], int, Dict[str, Any]) -> List[[SeqRecord, int]]

        allow_split_in_cds = get_value(kwargs, "allow_splits_in_cds", True)
        labels = get_value(kwargs, "labels", required=not allow_split_in_cds)

        interval_labels = None
        if not allow_split_in_cds and labels:
            interval_labels = GenomeSplitter.split_labels_into_intervals(labels)

        list_chunk_info = list()
        counter = 0

        for seqname, seqrecord in sequences.items():

            offset = 0
            while offset < len(seqrecord):
                # compute chunk boundaries
                left = offset
                right_excluded = min(offset + chunk_size_nt, len(seqrecord))

                while interval_labels and left > 0 and interval_labels.overlaps_point(right_excluded-1) \
                        and right_excluded < len(seqrecord):
                    lab = interval_labels[right_excluded-1].pop().data
                    right_excluded = lab.right() + 1      # skip

                chunk = seqrecord[left:right_excluded]  # type: SeqRecord

                # add offset information to chunk (to reassemble later)
                chunk.id += f"_offset={offset}"

                list_chunk_info.append([chunk, offset])
                offset = right_excluded
                counter += 1

        return list_chunk_info
    #
    # @staticmethod
    # def merge_predictions_with_offsets(list_prediction_info, pf_output):
    #     # type: (List[str, str, int], str) -> None
    #
    #     list_label = list()  # type: List[Label]
    #
    #     for pi in list_prediction_info:
    #         pf_pred, seqname, offset = pi
    #
    #         chunk_labels = read_labels_from_file(pf_pred, shift=0)
    #
    #         # add offset
    #         for lab in chunk_labels:
    #             lab.coordinates().left += offset
    #             lab.coordinates().right += offset
    #             lab.set_seqname(seqname)
    #             list_label.append(lab)
    #
    #     labels = Labels(list_label)
    #     write_labels_to_file(labels, pf_output, shift_coordinates_by=0)

    @staticmethod
    def split_labels_into_intervals(labels):
        # type: (Labels) -> IntervalTree
        list_intervals = [
            Interval(l.left(), l.right(), l) for l in labels if l.right() > l.left()
        ]

        return IntervalTree(list_intervals)



