import copy
from typing import *
import pandas as pd


class GenomeInfo:

    def __init__(self, name, genetic_code, attributes=None):
        # type: (str, int, Dict[str, Any]) -> None

        self.name = name
        self.genetic_code = genetic_code
        self.attributes = attributes


class GenomeInfoList:

    def __init__(self, list_genome_info):
        # type: (List[GenomeInfo]) -> None

        if list_genome_info is None:
            raise ValueError("List of genome names cannot be None")

        self._list_genome_info = copy.deepcopy(list_genome_info)

        # iterator
        self._iter_idx = 0
        self._iter_max = len(self._list_genome_info)

    @staticmethod
    def _make_list(data):
        if not isinstance(data, list):
            data = list(data)

        return data

    def __getitem__(self, item):
        return self._list_genome_info[item]

    def __iter__(self):
        # type: () -> GenomeInfoList
        return self

    def __len__(self):
        return len(self._list_genome_info)

    def __next__(self):
        # type: () -> GenomeInfo
        if self._iter_idx == self._iter_max:
            self._iter_idx = 0
            raise StopIteration

        element = self._list_genome_info[self._iter_idx]
        self._iter_idx += 1
        return element

    def next(self):
        return self.__next__()

    def to_file(self, pf_list):
        # type: (str) -> None

        def attributes_to_string(attributes):
            # type: (Dict[str, Any]) -> str
            if attributes is None or len(attributes) == 0:
                return "."
            else:
                return ";".join("{}={}".format(k, v) for k, v in attributes.items())

        df = pd.DataFrame({
            "gcfid": [x.name for x in self._list_genome_info],
            "genetic-code": [x.genetic_code for x in self._list_genome_info],
            "attributes": [attributes_to_string(x.attributes) for x in self._list_genome_info]
        })

        df.to_csv(pf_list, index=False)

    @classmethod
    def init_from_file(cls, pf_table):
        # type: (str) -> GenomeInfoList

        df = pd.read_csv(pf_table, header=0)

        list_genome_info = list()

        def parse_attributes(my_str):
            # type: (str) -> Dict[str, Any]

            loc_attributes = dict()
            for kv in my_str.strip().split(';'):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    loc_attributes[k] = v
            return loc_attributes

        for index, row in df.iterrows():

            gcfid = row["gcfid"]
            gcode = row["genetic-code"]
            attributes_str = str(row["attributes"])

            attributes = parse_attributes(attributes_str)

            list_genome_info.append(GenomeInfo(gcfid, gcode, attributes))

        return cls(list_genome_info)

