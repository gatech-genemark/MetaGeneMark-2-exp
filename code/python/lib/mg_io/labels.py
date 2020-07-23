import os
import re
import sys
from typing import *
from mg_general.general import get_value

sys.path.append(os.path.dirname(__file__) + "/..")  # add custom library directory to path

from mg_general.labels import Label, Coordinates, Labels


def write_string_to_file(astring, filename):
    with open(filename, "w") as f:
        f.write(astring)


def write_labels_to_file(labels, filename, shift_coordinates_by=1):
    out = labels.to_string(shift_coordinates_by=shift_coordinates_by)
    write_string_to_file(out, filename)


def create_attribute_dict(attribute_string, delimiter=";", key_value_delimiter="="):
    # type: (str, str) -> Dict[str, Any]

    attributes = dict()

    for current in attribute_string.strip().split(sep=delimiter):
        if len(current.strip().split(key_value_delimiter)) >= 2:
            k, v = current.strip().split(key_value_delimiter, maxsplit=1)
            attributes[k] = v

    return attributes


def read_labels_from_file(filename, shift=-1, name=None, **kwargs):
    # type: (str,  int, Union[str, None], Dict[str, Any]) -> Labels
    # FIXME: only supports gff

    ignore_frameshifted = get_value(kwargs, "ignore_frameshifted", True)
    ignore_partial = get_value(kwargs, "ignore_partial", True)
    tools = get_value(kwargs, "tools", None)
    key_value_delimiter = get_value(kwargs, "key_value_delimiter", "=", valid_type=str)

    labels = Labels(name=name)

    pattern = re.compile(r"([^\t]+)\t([^\t]+)\t(CDS)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t([^\t]+)\t([^\t]+)")

    with open(filename, "r") as f:

        for line in f:

            line = line.strip()

            m = pattern.match(line)
            if m:

                attributes = create_attribute_dict(m.group(9), key_value_delimiter=key_value_delimiter)

                label = Label.from_fields(
                    {
                        "left": int(m.group(4)) + shift,
                        "right": int(m.group(5)) + shift,
                        "strand": m.group(7),
                        "seqname": m.group(1),
                    },
                    attributes=attributes
                )

                if ignore_partial and label.is_partial():
                    continue

                if ignore_frameshifted and label.is_frameshifted():
                    continue

                if tools is not None and m.group(2) not in tools:
                    continue

                labels.add(label)
    return labels


def read_lst(pf_labels, shift=-1):
    # type: (str, int) -> Labels

    labels = list()

    # pattern = re.compile(r"([^\t]+)\t([^\t]+)\t(CDS)\t(\d+)\t(\d+)\t([^\t]+)\t([+-])\t([^\t]+)\t([^\t]+)")
    pattern = re.compile(r"([^\s]+)\s+([+-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(.+)$")

    # out = str(counter)
    # out += " " + str(l["strand"])
    # out += " " + str(l["left"] + shift)
    # out += " " + str(l["right"] + shift)
    # out += " " + str(l["right"] - l["left"] + 1)
    # out += " " "nativebac" + " AGGAGG 6 1"
    # out += " " + l["attributes"]["gene_type"]

    seqname = None

    with open(pf_labels, "r") as f:

        for line in f:

            line = line.strip()

            if line.startswith("SequenceID:"):
                seqname = line.split(":", maxsplit=1)[1].strip()
                continue
            elif len(line.strip()) == 0 or seqname is None:
                continue

            m = pattern.match(line)
            if m:
                attributes = m.group(6)

                label = {
                    "left": int(m.group(3)) + shift,
                    "right": int(m.group(4)) + shift,
                    "strand": m.group(2),
                    "seqname": seqname,
                    "attributes": attributes
                }

                labels.append(
                    Label(
                        Coordinates(
                            left=label["left"],
                            right=label["right"],
                            strand=label["strand"]
                        ),
                        seqname=seqname
                    )
                )

    return Labels(labels)
