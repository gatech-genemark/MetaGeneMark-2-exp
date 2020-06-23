import logging
from typing import *

from mg_general.general import get_value, except_if_not_valid
from mg_general.labels import Labels, create_key_3prime_from_label, Label, create_gene_key_from_label

logger = logging.getLogger(__name__)


class LabelsComparisonDetailed:

    def __init__(self, labels_a, labels_b, **kwargs):
        # type: (Labels, Labels, Dict[str, Any]) -> None

        self.labels_a = labels_a
        self.labels_b = labels_b

        self.name_a = get_value(kwargs, "name_a", "a", default_if_none=True)
        self.name_b = get_value(kwargs, "name_b", "b", default_if_none=True)

        self.tag = get_value(kwargs, "tag", None)

        self.comparison = dict()
        self._compare_labels(**kwargs)

    def _compare_labels(self, **kwargs):
        # type: (Dict[str, Any]) -> None

        split_on_attributes = get_value(kwargs, "split_on_attributes", None)

        self.comparison = {
            "all": self._compare_labels_helper(self.labels_a, self.labels_b),
            "attribute": dict()
        }

        if split_on_attributes is not None:
            for attribute_name in split_on_attributes:
                attribute_to_labels_a, attribute_to_labels_b = LabelsComparisonDetailed._split_labels_by_attribute(
                    self.labels_a, self.labels_b, attribute_name
                )
                self.comparison["attribute"][attribute_name] = dict()

                for attribute_value in attribute_to_labels_b.keys():
                    self.comparison["attribute"][attribute_name][attribute_value] = self._compare_labels_helper(
                        self.labels_a, attribute_to_labels_b[attribute_value]
                    )

    @staticmethod
    def _split_labels_by_attribute(labels_a, labels_b, attribute_name):
        # type: (Labels, Labels, str) -> Tuple[Dict[str, Labels], Dict[str, Labels]]
        attribute_value_to_labels_a = LabelsComparisonDetailed._split_labels_by_attribute_helper(
            labels_a, attribute_name
        )

        attribute_value_to_labels_b = LabelsComparisonDetailed._split_labels_by_attribute_helper(
            labels_b, attribute_name
        )

        all_attribute_values = set(attribute_value_to_labels_a.keys()).union(attribute_value_to_labels_b.keys())

        for v in all_attribute_values:
            if v not in attribute_value_to_labels_a:
                attribute_value_to_labels_a[v] = Labels()
            if v not in attribute_value_to_labels_b:
                attribute_value_to_labels_b[v] = Labels()

        return attribute_value_to_labels_a, attribute_value_to_labels_b

    @staticmethod
    def _split_labels_by_attribute_helper(labels, attribute_name):
        # type: (Labels, str) -> Dict[str, Labels]

        attribute_value_to_labels = dict()

        for l in labels:
            attribute_value = l.get_attribute_value(attribute_name)

            if attribute_value is not None:
                if attribute_value not in attribute_value_to_labels:
                    attribute_value_to_labels[attribute_value] = Labels()

                attribute_value_to_labels[attribute_value].add(l)

        return attribute_value_to_labels

    @staticmethod
    def _compare_labels_helper(labels_a, labels_b, **kwargs):
        # type: (Labels, Labels, Dict[str, Any]) -> Dict[str, Any]

        comparison = dict()

        key_3prime_to_label_a = {create_key_3prime_from_label(l): l for l in labels_a}
        key_3prime_to_label_b = {create_key_3prime_from_label(l): l for l in labels_b}

        compare_3p = LabelsComparisonDetailed._split_by_match_3prime(
            key_3prime_to_label_a, key_3prime_to_label_b
        )

        compare_3p_5p = LabelsComparisonDetailed._split_by_match_5prime(
            compare_3p["match"]
        )

        stats = dict()
        stats["total-a"] = len(labels_a)
        stats["total-b"] = len(labels_b)
        stats["match-3p"] = len(compare_3p["match"])
        stats["unique-3p-a"] = len(compare_3p["unique-a"])
        stats["unique-3p-b"] = len(compare_3p["unique-b"])
        stats["match-3p-5p"] = len(compare_3p_5p["match"])

        comparison["stats"] = stats

        comparison["labels"] = dict()
        comparison["labels"]["match-3p-5p"] = {
            "a": Labels([x[0] for x in compare_3p_5p["match"].values()]),
            "b": Labels([x[1] for x in compare_3p_5p["match"].values()])
        }
        comparison["labels"]["match-3p"] = {
            "a": Labels([x[0] for x in compare_3p["match"].values()]),
            "b": Labels([x[1] for x in compare_3p["match"].values()])
        }
        comparison["labels"]["match-3p-not-5p"] = {
            "a": Labels([x for x in compare_3p_5p["unique-a"].values()]),
            "b": Labels([x for x in compare_3p_5p["unique-b"].values()])
        }
        return comparison

    @staticmethod
    def _split_by_match_3prime(key_3prime_to_label_a, key_3prime_to_label_b):
        # type: (Dict[str, Label], Dict[str, Label]) -> Dict[str, Dict[str, Any]]
        result = {
            "match": dict(),
            "unique-a": dict(),
            "unique-b": dict()
        }

        keys_match = set(key_3prime_to_label_a.keys()).intersection(set(key_3prime_to_label_b.keys()))
        keys_unique_a = set(key_3prime_to_label_a.keys()).difference(set(key_3prime_to_label_b.keys()))
        keys_unique_b = set(key_3prime_to_label_b.keys()).difference(set(key_3prime_to_label_a.keys()))

        result["match"] = {
            key: (key_3prime_to_label_a[key], key_3prime_to_label_b[key]) for key in keys_match
        }

        result["unique-a"] = {key: key_3prime_to_label_a[key] for key in keys_unique_a}
        result["unique-b"] = {key: key_3prime_to_label_b[key] for key in keys_unique_b}

        return result

    @staticmethod
    def _split_by_match_5prime(key_to_pair_3p):
        # type: (Dict[str, Tuple(Label, Label)]) -> Dict[str, Dict[str, Any]]
        result = {
            "match": dict(),
            "unique-a": dict(),
            "unique-b": dict()
        }

        for key, pair in key_to_pair_3p.items():
            l_a, l_b = pair

            key_a = create_gene_key_from_label(l_a)
            key_b = create_gene_key_from_label(l_b)

            if key_a == key_b:
                result["match"][key] = (l_a, l_b)
            else:
                result["unique-a"][key] = l_a
                result["unique-b"][key] = l_b

        return result


    def _parse_compp_output(self, output):
        # type: (str) -> None
        d = self._split_compp_output(output)

        # put into reasonably-named keys
        good_and_bad = [
            ("found", "found_in_A"),
            ("identical", "identical_in_A")
        ]

        for letter in [("a", "A"), ("b", "B")]:
            g = letter[0]       # good
            b = letter[1]       # bad

            good_and_bad += [
                ("in-{}".format(g), "in_{}".format(b)),
                ("long-in-{}".format(g), "long_in_{}".format(b)),
                ("short-in-{}".format(g), "short_in_{}".format(b)),
                ("unique-in-{}".format(g), "unique_in_{}".format(b)),
            ]

        self.stats = {goodname: d[badname] for goodname, badname in good_and_bad}

    def intersection(self, source):
        # type: (str) -> Labels

        except_if_not_valid(source, {"a", "b"})

        return self.comparison["all"]["labels"]["match-3p-5p"][source]

    def match_3p_5p(self, source):
        # type: (str) -> Labels
        except_if_not_valid(source, {"a", "b"})
        return self.comparison["all"]["labels"]["match-3p-5p"][source]

    def match_3p_not_5p(self, source):
        # type: (str) -> Labels
        except_if_not_valid(source, {"a", "b"})
        return self.comparison["all"]["labels"]["match-3p-not-5p"][source]

    def match_3p(self, source):
        # type: (str) -> Labels
        except_if_not_valid(source, {"a", "b"})
        return self.comparison["all"]["labels"]["match-3p"][source]




