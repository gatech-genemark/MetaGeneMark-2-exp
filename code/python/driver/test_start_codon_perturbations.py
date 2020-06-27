# Author: Karl Gemayel
# Created: 6/24/20, 11:29 AM
import copy
import logging
import argparse
from tempfile import mktemp, mkstemp

import numpy as np
import pandas as pd
from typing import *
from tqdm import tqdm
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_container.gms2_mod import GMS2Mod
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import get_value, os_join, fix_names
from mg_general.labels_comparison_detailed import get_gene_start_difference
from mg_io.general import remove_p
from mg_models.shelf import run_gms2_prediction_with_model
from mg_parallelization.generic_threading import run_n_per_thread
from mg_viz import sns

parser = argparse.ArgumentParser("Test effect of perturbing start codon probabilities"
                                 " on the set of verified genes.")

parser.add_argument('--pf-gil', required=True, help="Genome list")

parser.add_argument('--alpha-mean', default=0, type=float, help="Mean of perturbation")
parser.add_argument('--alpha-std', default=0.1, type=float, help="Standard deviation of perturbation")
parser.add_argument('--num-perturbations', default=10, type=int, help="Number of experiments.")

parser.add_argument('--from-existing', default=False, action='store_true')
add_env_args_to_parser(parser)
parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment.init_from_argparse(parsed_args)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger


def perturb_start_codons(gms2_mod, p):
    # type: (GMS2Mod, float) -> None
    normalizer = 0

    start_codons = {"ATG", "GTG", "TTG"}
    for codon in start_codons:
        gms2_mod.items[codon] = max(float(gms2_mod.items[codon]) + p, 0.00000001)
        normalizer += gms2_mod.items[codon]

    # normalize to probabilities
    for codon in start_codons:
        gms2_mod.items[codon] /= normalizer


def test_start_codon_perturbation_for_gi(env, gi, p, gms2_mod):
    # type: (Environment, GenomeInfo, float) -> Dict[str, Any]

    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_verified = os_join(env["pd-data"], gi.name, "verified.gff")

    [f, pf_new_mod] = mkstemp(".mod")


    [f, pf_new_pred] = mkstemp(".gff")


    # perturb and write new model to file
    new_mod = copy.deepcopy(gms2_mod)
    perturb_start_codons(new_mod, p)
    new_mod.to_file(pf_new_mod)

    # run gms2 prediction step with new model
    run_gms2_prediction_with_model(env, pf_sequence, pf_new_mod, pf_new_pred)

    # compare with verified
    error = get_gene_start_difference(pf_verified, pf_new_pred)

    remove_p(pf_new_mod)
    remove_p(pf_new_pred)

    return {
        "Genome": gi.name,
        "Perturbation": p,
        **{f"Original {s}": gms2_mod.items[s] for s in {"ATG", "GTG", "TTG"}},
        **{f"Perturbed {s}": new_mod.items[s] for s in {"ATG", "GTG", "TTG"}},
        **error
    }

def run_perturbations(env, gi, perturbations, gms2_mod):
    list_entries = list()

    for p in tqdm(perturbations, f"{gi.name}", leave=False, total=len(perturbations)):
        entry = test_start_codon_perturbation_for_gi(env, gi, p, gms2_mod)
        entry["Perturbation"] = p
        entry["Genome"] = gi.name
        list_entries.append(entry)

    return list_entries

def test_start_codon_perturbations_for_gi(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> pd.DataFrame

    alpha_mean = get_value(kwargs, "alpha_mean", default=0, valid_type=float)
    alpha_std = get_value(kwargs, "alpha_std", default=0.1, valid_type=float)
    num_perturbations = get_value(kwargs, "num_perturbations", 20, valid_type=int)

    # read model file
    pf_mod = os_join(env["pd-runs"], gi.name, "gms2", "GMS2.mod")
    gms2_mod = GMS2Mod.init_from_file(pf_mod)

    # sample perturbations
    perturbations = np.random.normal(alpha_mean, alpha_std, num_perturbations)

    list_entries = list()

    # for p in tqdm(perturbations, f"{gi.name}", leave=False, total=len(perturbations)):
    #     entry = test_start_codon_perturbation_for_gi(env, gi, p, gms2_mod)
    #     list_entries.append(entry)

    list_entries = run_n_per_thread(
        [x for x in perturbations],
        test_start_codon_perturbation_for_gi,
        "p",
        {"gi": gi, "env": env, "gms2_mod": gms2_mod}
    )

    return pd.DataFrame(list_entries)


def test_start_codon_perturbations(env, gil, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, Any]) -> None

    from_existing = get_value(kwargs, "from_existing", False)

    if not from_existing:
        list_df = list()
        counter = 0
        for gi in tqdm(gil, "Genomes", total=len(gil)):
            list_df.append(test_start_codon_perturbations_for_gi(env, gi, **kwargs))
            counter += 1
            # if counter == 2:
            #     break

        df = pd.concat(list_df, ignore_index=True, sort=False)

        df["Genome"] = df.apply(fix_names, axis=1)
        df.to_csv(os_join(env["pd-work"], "summary.csv"), index=False)
    else:
        df = pd.read_csv(os_join(env["pd-work"], "summary.csv"))

    sns.catplot(df, "Genome", "Error")
    sns.lmplot(df, "Perturbation", "Error", hue="Genome", sns_kwargs={"lowess": True})

    import seaborn
    fig, ax = plt.subplots()
    for g in set(df["Genome"]):
        seaborn.regplot("Perturbation", "Original ATG", data=df[df["Genome"] == g], ax=ax, lowess=True, label=g)
        seaborn.regplot("Perturbation", "Perturbed ATG", data=df[df["Genome"] == g], ax=ax, lowess=True, label=g)
    plt.show()

    for s in {"ATG", "GTG", "TTG"}:
        df[s] = df[f"Perturbed {s}"] - df[f"Original {s}"]
        sns.lmplot(df, s, "Error", lowess=True, hue="Genome")

    df["Total"] = abs(df["ATG"]) + abs(df["GTG"]) + abs(df["TTG"])
    sns.lmplot(df, "Total", "Error", hue="Genome", lowess=True)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_gil)
    test_start_codon_perturbations(
        env, gil,
        alpha_mean=args.alpha_mean,
        alpha_std=args.alpha_std,
        num_perturbations=args.num_perturbations,
        from_existing=args.from_existing
    )


if __name__ == "__main__":
    main(my_env, parsed_args)
