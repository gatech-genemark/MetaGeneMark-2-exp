from typing import *
import logging

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment, SeqRecord, Seq

from mg_general.general import get_value
from mg_general.shelf import list_find_first
from mg_io.general import write_to_file

logger = logging.getLogger(__name__)


class MSASinglePointMarker:

    def __init__(self, position, msa_length, **kwargs):
        # type: (Union[int, None], int,  Dict[str, Any]) -> None

        self.name = get_value(kwargs, "name", "mark")

        self.mark_position = int(position) if position is not None and position >= 0 else None
        self.msa_length = int(msa_length)
        self.mark = get_value(kwargs, "mark", "M")
        self.gap = get_value(kwargs, "gap", "-")

    def to_string(self, begin=None, end=None):
        # type: (Union[int, None], Union[int, None]) -> str

        begin = begin if begin is not None else 0
        end = end if end is not None else self.msa_length

        if self.mark_position is None:
            return self.gap * (end-begin)

        # my_str = self.gap * self.mark_position + self.mark + self.gap * (self.msa_length - self.mark_position - 1)
        my_str = MSASinglePointMarker.create_mark_line(self.mark_position, self.msa_length, mark_tag=self.mark)

        if begin != 0 or end != self.msa_length:
            return my_str[begin:end]
        else:
            return my_str

    def change_symbol(self, new_symbol, old_symbol=None):
        # type: (str, str) -> None

        self.mark = new_symbol


    @staticmethod
    def create_mark_line(mark_position, length, **kwargs):
        # type: (Union[int, None], int, Dict[str, Any]) -> str

        mark_tag = get_value(kwargs, "mark_tag", "M", invalid={None})

        if len(mark_tag) != 1:
            raise ValueError("Mark tag ({}) should have length of 1".format(mark_tag))

        if mark_position is None:
            mark_sequence = "-" * length
        else:
            mark_sequence = "-" * mark_position + mark_tag + "-" * (length - mark_position - 1)

        return mark_sequence


class MSAType:

    def __init__(self, alignments, **kwargs):
        # type: (MultipleSeqAlignment, Dict[str, Any]) -> None

        self.list_msa_markers = get_value(kwargs, "list_msa_markers", list())      # type: List[MSASinglePointMarker]

        self.list_alignment_sequences = [s for s in alignments]           # type: List[SeqRecord]

    def get_mark_position(self, name):
        # type: (str) -> Union[int, None]

        mark = list_find_first(self.list_msa_markers, lambda x: x.name == name)  # type: MSASinglePointMarker

        if mark is None:
            raise ValueError("Unknown mark name ({})".format(name))

        return mark.mark_position

    def get_marker(self, name):
        # type: (str) -> MSASinglePointMarker
        return list_find_first(self.list_msa_markers, lambda x: x.name == name)

    def number_of_sequences(self):
        return len(self.list_alignment_sequences)

    def alignment_length(self):
        return len(self.list_alignment_sequences[0])

    def remove_markers_with_name(self, name):
        # type: (str) -> None

        new_marker_list = [m for m in self.list_msa_markers if m.name != name]
        self.list_msa_markers = new_marker_list

    def remove_markers_with_names(self, list_names):
        # type: (Iterable[str]) -> None

        for m in list_names:
            self.remove_markers_with_name(m)

    def add_marker(self, marker, unique=False):
        # type: (MSASinglePointMarker, bool) -> None
        if unique:
            self.remove_markers_with_name(marker.name)

        self.list_msa_markers.append(marker)

    def get_marker_names(self):
        # type: () -> List[str]

        return [m.name for m in self.list_msa_markers]

    def sort_by_field(self, field_number):
        # type: (int) -> None

        self.list_alignment_sequences = sorted(self.list_alignment_sequences,
                                               key=lambda x: x.id.split("-")[field_number])


    @staticmethod
    def separate_marks_from_alignment(alignment):
        # type: (MultipleSeqAlignment) -> Dict[str, Union[List[MSASinglePointMarker], MultipleSeqAlignment]]

        def find_first_non_gap(my_str):
            # type: (str) -> Union[int, None]

            for i in range(len(my_str)):
                if my_str[i] != "-":
                    return i

            return None

        list_alignments_without_marks = list()
        marks = list()
        msa_length = alignment.get_alignment_length()

        for a in alignment:
            if a.id[0] == "#":          # it's a mark
                mark_position = find_first_non_gap(a.seq._data)
                mark_name = a.id[1:] if len(a.id) > 1 else ""
                marks.append(MSASinglePointMarker(mark_position, msa_length, name=mark_name))

            else:
                list_alignments_without_marks.append(a)

        return {"marks": marks, "alignment": MultipleSeqAlignment(list_alignments_without_marks)}

    # @classmethod
    # def init_from_file_DEPRECATED(cls, pf_msa):
    #     # type: (str) -> MSAType
    #     alignment_info = sbsp_io.msa.read_msa_and_extract_selected_and_ref_positions(pf_msa)
    #
    #
    #
    #     # alignment_info = MSAType.separate_marks_from_alignment(alignment)
    #     alignment = MultipleSeqAlignment(alignment_info["alignment"])
    #     list_msa_marks = [
    #         MSASinglePointMarker(alignment_info["pos-ref"], alignment.get_alignment_length(), name="ref"),
    #         MSASinglePointMarker(alignment_info["pos-selected"], alignment.get_alignment_length(), name="selected"),
    #     ]
    #
    #
    #     return cls(alignment_info["alignment"], list_msa_markers=list_msa_marks)

    @staticmethod
    def init_from_file(pf_msa):
        # type: (str) -> MSAType

        def read_as_standard_clustal(pf_msa):
            # type: (str) -> Dict

            return AlignIO.read(pf_msa, "clustal")

        alignments_from_file = read_as_standard_clustal(pf_msa)

        alignments_processed_list = list()
        list_markers_info = list()

        for a in alignments_from_file:
            if a.id[0] != "#":
                alignments_processed_list.append(a)
            else:
                marker_name = a.id[1:]
                non_gap_positions = [x for x in range(len(a.seq._data)) if a.seq._data[x] != "-"]
                position = None
                mark_tag = "M"
                if len(non_gap_positions) > 0:
                    position = non_gap_positions[0]
                    mark_tag = a.seq._data[position]

                list_markers_info.append((marker_name, position, mark_tag))

        msa_t = MSAType(MultipleSeqAlignment(alignments_processed_list))
        for m in list_markers_info:
            msa_t.add_marker(MSASinglePointMarker(m[1], msa_length=msa_t.alignment_length(), name=m[0], mark=m[2]))

        return msa_t


    def __getitem__(self, item):
        # type: (int) -> SeqRecord
        return self.list_alignment_sequences[item]

    @staticmethod
    def _format_headers_pretty(headers):
        # type: (List[str]) -> List[str]

        def pad_by(curr_string, pad_len, padding=" "):
            # type: (str, int, str) -> str
            return curr_string + (padding * pad_len)

        headers_split = [
            line.strip().split(";") for line in headers
        ]

        # fix upstream distance:
        # If -1, set to -
        # convert rest to int
        for line in headers_split:
            val = int(float(line[1]))

            if val == -1:
                line[1] = "-"
            else:
                line[1] = "{}".format(val)

        # remove dN column
        headers_split = [
            line_split[:min(len(line_split), 5)] for line_split in headers_split
        ]

        max_num_columns = max(len(h) for h in headers_split)

        # max_length_per_header_column = list()
        # for col in range(max_num_columns):
        #     max_num_columns.append(
        #         max(len(l[col]) for l in headers_split)
        #     )

        max_length_per_header_column = [
            max(len(l[col]) for l in headers_split if col < len(l)) for col in range(max_num_columns)
        ]

        headers_pretty = list()

        for line_split in headers_split:
            headers_pretty.append(
                "  ".join(pad_by(line_split[i], max_length_per_header_column[i] - len(line_split[i]))
                        for i in range(len(line_split)))
            )

        return headers_pretty

    def change_marker(self, marker_name, new_symbol, old_symbol=None):
        # type: (str, str, str) -> None
        """
        Changes the marker of "marker_name" from the old symbol to a new symbol
        :param marker_name:
        :param new_symbol:
        :param old_symbol:
        :return: None
        """

        try:
            marker = self.get_marker(marker_name)
            marker.change_symbol(new_symbol)
        except ValueError:
            logger.warning("Cannot change marker for unknown name: {}".format(marker_name))




    def _to_string_pretty(self, begin=None, end=None, **kwargs):
        # type: (Union[int, None], Union[int, None], Dict[str, Any]) -> str

        tag = get_value(kwargs, "tag", "", default_if_none=True)

        self.change_marker("q3prime", new_symbol="*")

        # add markers as sequence records
        seq_records = [
            SeqRecord(Seq(m.to_string(begin, end)), id="#{}".format(m.name)) for m in self.list_msa_markers
        ]

        if begin is not None or end is not None:
            begin = begin if begin is not None else 0
            end = end if end is not None else self.alignment_length()

        headers_old = [a.id for a in self.list_alignment_sequences]
        headers_new = MSAType._format_headers_pretty(headers_old)

        # add actual sequences
        for i, a in enumerate(self.list_alignment_sequences):
            a.id = headers_new[i]

            if begin is not None or end is not None:
                seq_records.append(a[begin:end])
            else:
                seq_records.append(a)

        # create alignment with markers
        alignment = MultipleSeqAlignment(seq_records)

        output_string = alignment.format("clustal")

        # Remove header
        output_string = output_string.replace("_", " ")

        output_string_array = output_string.split("\n")

        def get_summary_statistics_line_for_alignment():
            # type: () -> str

            ref_position = self.get_mark_position("ref")

            is_lorf = len(set(self[0][0:ref_position])) <= 1

            def count_lorf_targets_near_position(position):
                # type: (int) -> int

                count = 0

                for idx in range(1, self.number_of_sequences()):
                    j = 0
                    while True:
                        if position-j >= 0 and self[idx][position-j].isupper():
                            if len(set(self[idx][0:position-j])) <= 1:
                                count += 1
                            break
                        if position+j < self.alignment_length() and self[idx][position+j].isupper():
                            if len(set(self[idx][0:position+j])) <= 1:
                                count += 1
                            break

                        j += 1

                return count

            num_targets_that_are_lorf = count_lorf_targets_near_position(ref_position)

            return "{}: LORF={}\tTargetLORF={}".format(tag, str(is_lorf)[0], num_targets_that_are_lorf)

        output_string_array[0] = get_summary_statistics_line_for_alignment()

        output_string = "\n".join(output_string_array)

        return output_string



    def to_string(self, begin=None, end=None, **kwargs):
        # type: (Union[int, None], Union[int, None], Dict[str, Any]) -> str

        format = get_value(kwargs, "format", None)

        if format == "pretty":
            return self._to_string_pretty(begin, end, **kwargs)

        # add markers as sequence records
        seq_records = [
            SeqRecord(Seq(m.to_string(begin, end)), id="#{}".format(m.name)) for m in self.list_msa_markers
        ]

        if begin is not None or end is not None:
            begin = begin if begin is not None else 0
            end = end if end is not None else self.alignment_length()

        # add actual sequences
        for a in self.list_alignment_sequences:

            if begin is not None or end is not None:
                seq_records.append(a[begin:end])
            else:
                seq_records.append(a)

        # create alignment with markers
        alignment = MultipleSeqAlignment(seq_records)

        return alignment.format("clustal")

    def to_file(self, pf_out, begin=None, end=None, format=None, **kwargs):
        # type: (str, Union[int, None], Union[int, None]) -> None

        write_to_file(self.to_string(begin, end, format=format, **kwargs), pf_out)


