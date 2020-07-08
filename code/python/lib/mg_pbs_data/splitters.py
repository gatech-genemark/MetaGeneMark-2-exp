import math
import os
import logging
import pandas as pd
from typing import *

from mg_container.genome_list import GenomeInfoList
from mg_general.general import get_value, all_elements_equal

log = logging.getLogger(__name__)

T = TypeVar('T')


def get_number_elements_per_split(total_elements, num_splits):
    # type: (int, int) -> int
    """Compute the number of elements in each split"""
    return math.ceil(total_elements / float(num_splits))


def generate_splits(elements, num_splits):
    # type: (Sized[T], int) -> Generator[T]
    """Evenly split elements into 'num_splits' groups, and yield one by one"""

    num_per_split = get_number_elements_per_split(len(elements), num_splits)

    num_elements = len(elements)
    start = 0

    while start < num_elements:
        end = min(num_elements, start + num_per_split)
        yield elements[start:end]
        start = end


def split_gil(data, num_splits, **kwargs):
    # type: (Dict[str, Any], int, Dict[str, Any]) -> List[Dict[str, Any]]
    pf_output_template = get_value(data, "pf_output_template", None, value_type=str)
    arg_name_pf_output = get_value(data, "arg_name_pf_output", "pf_output", value_type=str)
    arg_name_gil = get_value(kwargs, "arg_name_gil", "gil", value_type=str)
    arg_name_jobid = get_value(kwargs, "arg_name_jobid", None, value_type=str)

    gil = data


    list_of_list_of_gi = list()
    for i in range(num_splits):
        list_of_list_of_gi.append(list())

    for index, gi in enumerate(gil):
        index_of_list = index % num_splits

        list_of_list_of_gi[index_of_list].append(gi)

    list_output = list()
    for i in range(len(list_of_list_of_gi)):
        val = {
            arg_name_gil: GenomeInfoList(list_of_list_of_gi[i]),
        }

        if pf_output_template is not None:
            val[arg_name_pf_output]: pf_output_template.format(i)
        if arg_name_jobid is not None:
            val[arg_name_jobid] = i

        list_output.append(val)

    return list_output


def split_list(data, num_splits, **kwargs):
    # type: (Dict[str, Any], int, Dict[str, Any]) -> List[Dict[str, Any]]
    pf_output_template = get_value(data, "pf_output_template", None, value_type=str)
    arg_name_pf_output = get_value(data, "arg_name_pf_output", "pf_output", value_type=str)
    arg_name_data = get_value(kwargs, "arg_name_data", "data", value_type=str)
    arg_name_jobid = get_value(kwargs, "arg_name_jobid", None, value_type=str)


    list_of_lists = list()
    for i in range(num_splits):
        list_of_lists.append(list())

    for index, item in enumerate(data):
        index_of_list = index % num_splits

        list_of_lists[index_of_list].append(item)

    list_output = list()
    for i in range(len(list_of_lists)):
        val = {
            arg_name_data: list_of_lists[i],
        }

        if pf_output_template is not None:
            val[arg_name_pf_output]: pf_output_template.format(i)

        if arg_name_jobid is not None:
            val[arg_name_jobid] = i

        list_output.append(val)

    return list_output







# old ones

def split_list_DEPRECATED(data, num_splits, pd_work, **kwargs):
    # type: (Dict[str, Any], int, str, Dict[str, Any]) -> List[Dict[str, str]]

    list_pf_data = data["list_pf_data"]
    pf_output_template = data["pf_output_template"]

    list_splits = list()

    split_number = 1
    for v in list_pf_data:
        list_splits.append({"pf_data": v, "pf_output": pf_output_template.format(split_number)})
        split_number += 1

    return list_splits


def split_dict(data, num_splits, pd_work, **kwargs):
    # type: (Dict[str, Any], int, str, Dict[str, Any]) -> List[Dict[str, Any]]

    a_dict = data["dict"]  # type: Dict[str, Any]
    pf_output_template = data["pf_output_template"]

    list_splits = list()
    num_splits = min(num_splits, len(a_dict))
    for i in range(num_splits):
        list_splits.append(dict())
        list_splits[-1]["data"] = dict()

    index = 0

    for k, v in a_dict.items():
        list_splits[index % num_splits]["data"][k] = v
        index += 1

    for split_number in range(1, len(list_splits) + 1):
        list_splits[split_number - 1]["pf_output"] = pf_output_template.format(split_number)
        list_splits[split_number - 1]["msa_output_start"] = split_number
        list_splits[split_number - 1]["fn_tmp_prefix"] = split_number

    return list_splits


def split_genome_info_list(data, num_splits, pd_work, **kwargs):
    # type: (Dict[str, Any], int, str, Dict[str, Any]) -> List[Dict[str, Any]]

    genome_info_list = get_value(data, "gil", required=True)

    pf_output_template = get_value(data, "pf_output_template", "")

    if num_splits > len(genome_info_list):
        num_splits = len(genome_info_list)

    list_of_list_of_gi = list()
    for i in range(num_splits):
        list_of_list_of_gi.append(list())

    for index, gi in enumerate(genome_info_list):
        index_of_list = index % num_splits

        list_of_list_of_gi[index_of_list].append(gi)

    return [
        {
            "gil": GenomeInfoList(list_of_list_of_gi[i]),
            "pf_output": pf_output_template.format(i)
        } for i in range(len(list_of_list_of_gi))
    ]


def split_q3prime_files(data, num_splits, pd_work, **kwargs):
    # type: (Dict[str, Any], int, str, Dict[str, Any]) -> List[Dict[str, str]]

    file_number = 1

    list_splits = list()
    q3prime_to_list_pf = data["q3prime_to_list_pf"]
    pf_output_template = data["pf_output_template"]

    for q3prime_key in q3prime_to_list_pf.keys():
        list_pf = q3prime_to_list_pf[q3prime_key]
        list_splits.append({"list_pf_data": list_pf, "pf_output": pf_output_template.format(file_number),
                            "q3prime": q3prime_key,
                            "msa_output_start": file_number})

        file_number += 1

    return list_splits


def split_list_and_remerge_by_key(data, num_splits, pd_work, **kwargs):
    # type: (Dict[str, Any], int, str, Dict[str, Any]) -> List[Dict[str, str]]

    # def copy_files_to_new_dir_by_group(list_pf_data, pf_new_data_format, group_key):
    # type: # (List[str], str, str) -> List[str]

    # goal: move data from one directory to another, such that new files are one per group_key, and
    # this is done memory-efficiently

    # sketch: go through files one by one, create a new file every time we get to a new key.
    file_number = 1
    group_to_file_number = dict()

    list_pf_data = data["list_pf_data"]
    group_key = data["group_key"]
    pf_output_template = data["pf_output_template"]

    list_pf_new = list()

    for pf_old in list_pf_data:

        try:
            df_old = pd.read_csv(pf_old, header=0)

            # loop over groups
            if group_key in df_old:
                for name, df_group in df_old.groupby(group_key):

                    # if group hasn't been tapped yet, create a new file for it
                    if name not in group_to_file_number.keys():
                        pf_new = pf_output_template.format(file_number)
                        group_to_file_number[name] = file_number

                        list_pf_new.append({"pf_data": pf_new, "pf_output": pf_output_template.format(file_number),
                                            "msa_output_start": file_number})
                        file_number += 1

                        df_new = df_group
                    else:

                        curr_file_number = group_to_file_number[name]
                        pf_new = pf_output_template.format(curr_file_number)

                        df_new = pd.read_csv(pf_new, header=0)
                        df_new = df_new.append(df_group)

                    df_new.to_csv(pf_new, index=False)


        except IOError:
            pass

    return list_pf_new


def split_query_genomes_target_genomes_one_vs_group(data, num_splits, pd_work, **kwargs):
    # type: (Dict[str, str], int, str, Dict[str, Any]) -> List[Dict[str, str]]

    fn_q_split_formatted = get_value(kwargs, "fn_q_split_formatted", "q_split_{}.list")
    fn_t_split_formatted = get_value(kwargs, "fn_t_split_formatted", "t_split_{}.list")

    pf_q_list = data["pf_q_list"]
    pf_t_list = data["pf_t_list"]
    pf_output_template = data["pf_output_template"]

    list_pf_splits = list()

    q_list = GenomeInfoList.init_from_file(pf_q_list)
    t_list = GenomeInfoList.init_from_file(pf_t_list)

    split_number = 1
    for q_genome in q_list:

        # split
        for split in generate_splits(t_list, num_splits):
            pf_q_split = os.path.join(pd_work, fn_q_split_formatted.format(split_number))
            pf_t_split = os.path.join(pd_work, fn_t_split_formatted.format(split_number))

            q_split = GenomeInfoList([q_genome])
            t_split = GenomeInfoList(split)

            q_split.to_file(pf_q_split)
            t_split.to_file(pf_t_split)

            list_pf_splits.append({
                "pf_q_list": pf_q_split,
                "pf_t_list": pf_t_split,
                "pf_output": pf_output_template.format(split_number)
            })

            split_number += 1

    return list_pf_splits
