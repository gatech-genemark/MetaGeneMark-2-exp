import copy
from typing import *

from mg_general.general import create_gene_key, get_value


class Coordinates:

    def __init__(self, left=None, right=None, strand=None):
        self.left = left
        self.right = right
        self.strand = strand

    def to_string(self, field=None, shift_coordinates_by=0, delim="\t"):

        s = int(shift_coordinates_by)

        def stringify(element):
            if element is None:
                return ""
            return str(element)

        if field is None:
            return stringify(int(self.left) + s) + delim + stringify(int(self.right) + s) + delim + stringify(
                self.strand)

        if field == "left":
            return stringify(int(self.left) + s)

        if field == "right":
            return stringify(int(self.right) + s)

        if field == "strand":
            return stringify(self.strand)

        raise ValueError("Unrecognized field: " + stringify(field))

    def get_5prime(self):
        return self.left if self.strand == "+" else self.right

    def get_3prime(self):
        return self.right if self.strand == "+" else self.left

    @classmethod
    def from_fields(cls, fields):
        # type: (dict) -> Coordinates
        left = None
        right = None
        strand = None

        if "left" in fields:
            left = fields["left"]
        if "right" in fields:
            right = fields["right"]
        if "strand" in fields:
            strand = fields["strand"]

        return cls(left, right, strand)


class Label:

    def __init__(self, coordinates, seqname=None, **kwargs):
        # type: (Coordinates, str, Dict[str, Any]) -> None

        meta = get_value(kwargs, "meta", dict(), default_if_true=True)

        self._fields = dict()
        self._fields["seqname"] = seqname
        self._fields["coordinates"] = coordinates

        self.meta = meta

        # set default values for Nones
        for f in self._fields:
            if self._fields[f] is None:
                self._fields[f] = Label.default_for_key(f)

        self._attributes = get_value(kwargs, "attributes", dict(), defualt_if_none=True)

    def left(self):
        # type: () -> int
        return self.coordinates().left

    def right(self):
        # type: () -> int
        return self.coordinates().right

    def to_string(self, shift_coordinates_by=0):
        # FIXME: prints in GFF format only for now
        # seqname, source, feature, left, right, score, strand, frame, attributes

        out = self._fields["seqname"]
        out += "\t" "."
        out += "\t" + "CDS"
        out += "\t" + str(self._fields["coordinates"].to_string("left", shift_coordinates_by))
        out += "\t" + str(self._fields["coordinates"].to_string("right", shift_coordinates_by))
        out += "\t" "."
        out += "\t" + str(self._fields["coordinates"].to_string("strand"))
        out += "\t" + "."

        def attribute_to_strings():
            out = ""

            prefix = ""
            for k, v in self._attributes.items():
                out += "{}{}={}".format(prefix, k, v)
                prefix = ";"

            return out

        if len(self._attributes) > 0:
            out += "\t" + attribute_to_strings()
        else:
            out += "\t" + "."

        return out

    def get_attribute_value(self, attribute):
        # type: (str) -> Any

        return self._attributes[attribute] if attribute in self._attributes else None

    def set_attribute_value(self, attribute, value):
        # type: (str, Any) -> None
        self._attributes[attribute] = value

    def is_partial(self):
        # type: () -> bool

        if "partial" in self._attributes:
            if self._attributes["partial"].lower() == "true" or self._attributes["partial"] == True:
                return True

            if self._attributes["partial"] in {"01", "11", "10"}:
                return True

        return False

    def incomplete_at_5prime(self):
        # type: () -> bool
        if "partial" in self._attributes and len(self._attributes["partial"]) == 2:
            partial = self._attributes["partial"]
            if self.strand() == "+" and partial[0] == "1":
                return True
            if self.strand() == "-" and partial[1] == "1":
                return True

        return False

    def incomplete_at_3prime(self):
        # type: () -> bool
        if "partial" in self._attributes and len(self._attributes["partial"]) == 2:
            partial = self._attributes["partial"]
            if self.strand() == "+" and partial[1] == "1":
                return True
            if self.strand() == "-" and partial[0] == "1":
                return True

        return False


    def is_hypothetical(self):
        # type: () -> bool
        if self.get_attribute_value("product") is not None and "hypothetical" in self.get_attribute_value("product"):
            return True

        return False

    def is_frameshifted(self):
        # type: () -> bool

        if self.is_partial():
            return False

        length = self.coordinates().right - self.coordinates().left + 1
        if length % 3 != 0:
            return True

        return False

    @classmethod
    def minimum_set_of_field_names(cls):

        return {"seqname", "coordinates"}

    @classmethod
    def from_fields(cls, fields, **kwargs):
        #

        def get_keys_if_exist_or_none(keys, fields):
            new_fields = fields.copy()

            for k in keys:
                if k == "coordinates":
                    new_fields[k] = Coordinates.from_fields(fields)
                elif k not in new_fields:
                    new_fields[k] = None

            return new_fields

        minimum_keys = Label.minimum_set_of_field_names()
        new_fields = get_keys_if_exist_or_none(minimum_keys, fields)

        return cls(new_fields["coordinates"], new_fields["seqname"], **kwargs)

    @classmethod
    def default_for_key(cls, key):
        if key == "seqname":
            return ""
        if key == "coordinates":
            return Coordinates()

        raise ValueError("Unknown key: " + str(key))

    def get_3prime(self):
        if self._fields["coordinates"] is not None:
            return self._fields["coordinates"].get_3prime()
        return None

    def get_5prime(self):
        if self._fields["coordinates"] is not None:
            return self._fields["coordinates"].get_5prime()
        return None

    def strand(self):
        return self._fields["coordinates"].strand

    def coordinates(self):
        # type: () -> Coordinates
        return self._fields["coordinates"]

    def seqname(self):
        return self._fields["seqname"]

    def set_seqname(self, seqname):
        # type: (str) -> None
        self._fields["seqname"] = seqname

    def length(self):
        return self.coordinates().right - self.coordinates().left + 1


class Labels:

    def __init__(self, labels=None, name=None):
        # type: (List[Label], str) -> None

        if labels is None:
            labels = list()

        self._labels = copy.copy(labels)  # type: List[Label]
        self._labels_by_3p = {create_key_3prime_from_label(lab): lab for lab in labels}
        self.name = name

        if not isinstance(self._labels, list):
            self._labels = list(self._labels)

        # iterator
        self._iter_idx = 0
        self._iter_max = len(self._labels)

    def __iter__(self):
        # type: () -> Labels
        return self

    def __len__(self):
        return len(self._labels)

    def get_by_3prime_key(self, key):
        # type: (str) -> Union[Label, None]
        return self._labels_by_3p.get(key)

    def get_multiple_by_3prime_keys(self, keys):
        # type: (Iterable[str]) -> Labels

        labels = Labels(name=self.name)

        for l in self._labels:
            if create_key_3prime_from_label(l) in keys:
                labels.add(l)

        return labels

    def __next__(self):
        # type: () -> Label
        if self._iter_idx == self._iter_max:
            self._iter_idx = 0
            raise StopIteration

        element = self._labels[self._iter_idx]
        self._iter_idx += 1
        return element

    next = __next__  # python 2

    def add(self, label):
        self._labels.append(label)
        self._labels_by_3p[create_key_3prime_from_label(label)] = label

        self._iter_max = len(self._labels)

    def add_multiple(self, labels):
        # type: (Labels) -> None
        self._labels += labels._labels

        for l in labels._labels:
            self._labels_by_3p[create_key_3prime_from_label(l)] = l

        self._iter_max = len(self._labels)

    def update(self, new_labels):
        # type: (Labels) -> Labels

        dict_label_by_3prime_key = {
            create_key_3prime_from_label(l): l for l in self._labels
        }

        for l in new_labels:
            key = create_key_3prime_from_label(l)

            dict_label_by_3prime_key[key] = l

        return Labels([v for _, v in dict_label_by_3prime_key.items()], name=self.name)

    def to_string(self, shift_coordinates_by=0):

        out = ""
        for n in range(self._iter_max):
            out += self._labels[n].to_string(shift_coordinates_by) + "\n"

        return out

    def to_string_lst(self, shift_coordinates_by=0):

        out = "# GeneMark.hmm-2 LST format\n"
        out += "# GeneMark.hmm-2 prokaryotic version: 1.14\n"
        out += "# File with sequence: tmpseq.fna\n"
        out += "# File with native parameters: itr_1.mod\n"
        out += "# Native species name and build: gms2-training\n"
        out += "# File with MetaGeneMark parameters: /storage4/karl/sbsp/biogem/sbsp/bin_external/gms2/mgm_11.mod\n"
        out += "# translation table: 11\n"
        out += "# output date start: Mon Jun  8 09:26:44 2020\n\n"

        seqname_to_labels = dict()
        for l in self._labels:  # type: Label
            if l.seqname() not in seqname_to_labels:
                seqname_to_labels[l.seqname()] = list()

            seqname_to_labels[l.seqname()].append(l)

        for seqname, seqname_labels in seqname_to_labels.items():
            out += f"SequenceID: {seqname}\n"
            counter = 1

            for counter, l in enumerate(seqname_labels):
                out += str(counter)
                out += " " + str(l.strand())
                out += " " + str(l.left() + shift_coordinates_by)
                out += " " + str(l.right() + shift_coordinates_by)
                out += " " + str(l.right() - l.left() + 1)
                out += " " "nativebac" + " AGGAGG 6 1"
                out += " " + " ."

                out += "\n"

        return out

    def __getitem__(self, index):
        # type: (int) -> Label
        return self._labels[index]

    def __str__(self):
        return self.to_string()

    def get_labels_with_matching_3prime(self, labels_reference):
        # type: (Labels) -> Labels

        keys_3prime = {create_key_3prime_from_label(l) for l in labels_reference}
        return Labels(
            [l for l in self._labels if create_key_3prime_from_label(l) in keys_3prime]
        )

    def get_labels_without_matching_3prime(self, labels_reference):
        # type: (Labels) -> Labels

        keys_3prime = {create_key_3prime_from_label(l) for l in labels_reference}
        return Labels(
            [l for l in self._labels if create_key_3prime_from_label(l) not in keys_3prime]
        )

    def get_labels_on_strand(self, strand):
        # type: (str) -> Labels

        return Labels(
            [l for l in self if l.strand() == strand],
            name=self.name
        )

    def sort_by(self, coordinate, in_place=False):
        # type: (str, bool) -> Labels

        if in_place:
            raise NotImplemented()

        if coordinate == "left":
            return Labels(sorted(self._labels, key=lambda l: (l.seqname(), l.coordinates().left)), name=self.name)

        else:
            return Labels(sorted(self._labels, key=lambda l: (l.seqname(), l.coordinates().right)), name=self.name)

    def sort_by_attribute(self, attribute, in_place=False):
        # type: (str, bool) -> Labels
        if in_place:
            raise NotImplementedError()

        return Labels(sorted(self._labels, key=lambda l: l.get_attribute_value(attribute)), name=self.name)

    @staticmethod
    def _split_by_overlap_helper(labels, labels_reference, strand):
        # type: (Labels, Labels, str) -> [Labels, Labels]

        # assumption that labels are sorted by left coordinate

        labels_with_overlap = Labels(name=labels.name)
        labels_with_no_overlap = Labels(name=labels.name)

        key_to_index = {create_key_3prime_from_label(labels_reference[i]): i for i in range(len(labels_reference))}

        labels_combined = labels_reference
        labels_combined.add_multiple(Labels(
            [l for l in labels if create_key_3prime_from_label(l) not in key_to_index]
        ))

        if strand == "+":
            labels_combined = labels_combined.sort_by(
                "right")  # sorted(labels_combined, key=lambda l: l.coordinates().right)
        else:
            labels_combined = labels_combined.sort_by(
                "left")  # sorted(labels_combined, key=lambda l: l.coordinates().left)

        key_to_index = {create_key_3prime_from_label(labels_combined[i]): i for i in range(len(labels_combined))}

        for l in labels:
            key = create_key_3prime_from_label(l)

            index = key_to_index[key]

            if l.strand() == "+":
                if index == 0:
                    labels_with_no_overlap.add(l)
                else:
                    if l.coordinates().left - labels_combined[index - 1].coordinates().right <= 0:
                        labels_with_overlap.add(l)
                    else:
                        labels_with_no_overlap.add(l)

            # negative strand
            else:
                if index == len(labels_combined) - 1:
                    labels_with_no_overlap.add(l)
                else:
                    if labels_combined[index + 1].coordinates().left - l.coordinates().right <= 0:
                        labels_with_overlap.add(l)
                    else:
                        labels_with_no_overlap.add(l)

        return labels_with_overlap, labels_with_no_overlap

    def split_by_overlap_status(self, labels_reference):
        # type: (Labels) -> [Labels, Labels]

        labels_reference_sorted_by_left = labels_reference.sort_by(
            "left")  # sorted(labels_reference, key=lambda l: l.coordinates().left)
        labels_reference_sorted_by_right = labels_reference.sort_by(
            "right")  # sorted(labels_reference, key=lambda l: l.coordinates().right)

        labels_self_positive_sorted_by_left = self.get_labels_on_strand("+")
        labels_self_negative_sorted_by_right = self.get_labels_on_strand("-")

        labels_self_positive_overlap, labels_self_positive_no_overlap = Labels._split_by_overlap_helper(
            labels_self_positive_sorted_by_left, labels_reference_sorted_by_right, strand="+"
        )

        labels_self_negative_overlap, labels_self_negative_no_overlap = Labels._split_by_overlap_helper(
            labels_self_negative_sorted_by_right, labels_reference_sorted_by_left, strand="-"
        )

        labels_self_overlap = labels_self_positive_overlap
        labels_self_overlap.add_multiple(labels_self_negative_overlap)

        labels_self_no_overlap = labels_self_positive_no_overlap
        labels_self_no_overlap.add_multiple(labels_self_negative_no_overlap)

        return labels_self_overlap, labels_self_no_overlap

    @staticmethod
    def _compute_distance_to_upstream_genes_on_positive_strand(labels):
        # type: (Labels) -> Dict[str, Union[int, None]]

        dict_key_to_distance = dict()

        labels = labels.sort_by("left", in_place=False)

        for i in range(len(labels)):
            lab = labels[i]

            if lab.strand() == "-":
                continue

            key = create_gene_key_from_label(lab)

            # search previous 5 genes (for sanity) for the largest right coordinate
            largest_right = None
            for j in range(5):
                index_of_previous = i - j - 1
                if index_of_previous < 0:
                    break

                if largest_right is None:
                    largest_right = labels[index_of_previous].coordinates().right
                else:
                    if largest_right < labels[index_of_previous].coordinates().right:
                        largest_right = labels[index_of_previous].coordinates().right

            distance = lab.coordinates().left - largest_right if largest_right is not None else None

            dict_key_to_distance[key] = distance

        return dict_key_to_distance

    @staticmethod
    def _compute_distance_to_upstream_genes_on_negative_strand(labels):
        # type: (Labels) -> Dict[str, Union[int, None]]
        dict_key_to_distance = dict()

        labels = labels.sort_by("right", in_place=False)

        for i in range(len(labels)):
            i = len(labels) - i - 1  # to reverse index

            lab = labels[i]

            if lab.strand() == "+":
                continue

            key = create_gene_key_from_label(lab)

            # search previous 5 genes (for sanity) for the smallest left coordinate
            smallest_left = None
            for j in range(5):
                index_of_previous = i + j + 1
                if index_of_previous >= len(labels):
                    break

                if smallest_left is None:
                    smallest_left = labels[index_of_previous].coordinates().left
                else:
                    if smallest_left > labels[index_of_previous].coordinates().left:
                        smallest_left = labels[index_of_previous].coordinates().left

            distance = smallest_left - lab.coordinates().right if smallest_left is not None else None

            dict_key_to_distance[key] = distance

        return dict_key_to_distance

    def compute_distance_to_upstream_gene(self):
        # type: () -> Dict[str, int]

        labels_copy = Labels(copy.copy(self._labels))

        distances_pos_strand = Labels._compute_distance_to_upstream_genes_on_positive_strand(labels_copy)
        distances_neg_strand = Labels._compute_distance_to_upstream_genes_on_negative_strand(labels_copy)

        distances_neg_strand.update(distances_pos_strand)

        return distances_neg_strand


def create_gene_key_from_label(label, genome_name=None):
    # type: (sbsp_general.labels.Label, str) -> str

    return create_gene_key(genome_name, label.seqname(), label.coordinates().left, label.coordinates().right,
                           label.coordinates().strand)


def create_key_3prime_from_label(label, genome_name=None):
    # type: (Label, Union[str, None]) -> str
    if label.strand() == "+":
        return "{};{};{};{};{}".format(genome_name, label.seqname(), "",
                                       label.coordinates().right, label.strand())
    else:
        return "{};{};{};{};{}".format(genome_name, label.seqname(), label.coordinates().left,
                                       "", label.strand())


def shift_5prime(label, amount):
    # type: (Label, int) -> None

    """

    :param label: label
    :param amount: positive means downstream, negative means upstream
    :return:
    """

    if abs(amount) % 3 != 0:
        import logging
        logging.debug("Shifting 5prime by value ({}) not a multiple of 3".format(amount))

    if label.strand() == "+":
        label.coordinates().left += amount
    else:
        label.coordinates().right -= amount
