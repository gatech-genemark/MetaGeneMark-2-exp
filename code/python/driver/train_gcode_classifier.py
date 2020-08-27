# Author: Karl Gemayel
# Created: 8/27/20, 11:53 AM

import logging
import argparse
import os
from collections import defaultdict
import random

import pandas as pd
from typing import *
import numpy as np

# noinspection All
from sklearn.preprocessing import StandardScaler

import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_io.general import save_obj, load_obj
from mg_viz import sns

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-data', required=True)
parser.add_argument('--pf-checkpoint')

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


def stratified_group_k_fold(X, y, groups, k, seed=None):
    labels_num = np.max(y) + 1
    y_counts_per_group = defaultdict(lambda: np.zeros(labels_num))
    y_distr = Counter()
    for label, g in zip(y, groups):
        y_counts_per_group[g][label] += 1
        y_distr[label] += 1

    y_counts_per_fold = defaultdict(lambda: np.zeros(labels_num))
    groups_per_fold = defaultdict(set)

    def eval_y_counts_per_fold(y_counts, fold):
        y_counts_per_fold[fold] += y_counts
        std_per_label = []
        for label in range(labels_num):
            label_std = np.std([y_counts_per_fold[i][label] / y_distr[label] for i in range(k)])
            std_per_label.append(label_std)
        y_counts_per_fold[fold] -= y_counts
        return np.mean(std_per_label)

    groups_and_y_counts = list(y_counts_per_group.items())
    random.Random(seed).shuffle(groups_and_y_counts)

    for g, y_counts in sorted(groups_and_y_counts, key=lambda x: -np.std(x[1])):
        best_fold = None
        min_eval = None
        for i in range(k):
            fold_eval = eval_y_counts_per_fold(y_counts, i)
            if min_eval is None or fold_eval < min_eval:
                min_eval = fold_eval
                best_fold = i
        y_counts_per_fold[best_fold] += y_counts
        groups_per_fold[best_fold].add(g)

    all_groups = set(groups)
    for i in range(k):
        train_groups = all_groups - groups_per_fold[i]
        test_groups = groups_per_fold[i]

        train_indices = [i for i, g in enumerate(groups) if g in train_groups]
        test_indices = [i for i, g in enumerate(groups) if g in test_groups]

        yield train_indices, test_indices

def get_distribution(y_vals):
    y_distr = Counter(y_vals)
    y_vals_sum = sum(y_distr.values())
    return [f'{y_distr[i] / y_vals_sum:.2%}' for i in range(np.max(y_vals) + 1)]




def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    from sklearn import datasets
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import GridSearchCV

    df = pd.read_csv(args.pf_data)
    # df = df.sample(100000)
    df["True Gcode"].replace({11: 0, 4: 1}, inplace=True)
    y = df["True Gcode"].values.astype(int)
    groups = df["Genome"].values
    df.drop(["Genome", "True Gcode"], axis=1, inplace=True)
    df.fillna(0, inplace=True)
    X = df.values

    pf_checkpoint = args.pf_checkpoint
    if not pf_checkpoint or not os.path.isfile(pf_checkpoint):
        # Define a pipeline to search for the best combination of PCA truncation
        # and classifier regularization.
        pca = PCA()
        # set the tolerance to a large value to make the example faster
        logistic = LogisticRegression(max_iter=10000, tol=0.1, solver='sag', verbose=5)
        pipe = Pipeline(steps=[('scaler', StandardScaler()), ('logistic', logistic)])

        # df = df.sample(100000)

        # Parameters of pipelines can be set using ‘__’ separated parameter names:
        param_grid = {
            # 'pca__n_components': [5,],
            'logistic__C': np.logspace(-4, 4, 4)
        }

        search = GridSearchCV(
            pipe, param_grid, n_jobs=-1, cv=stratified_group_k_fold(X, y, groups, k=5),
            scoring='f1', verbose=5)

        logger.critical("Fitting")
        search.fit(X, y)
        logger.critical("Prediction")
        # grid = search.cv_results_
        model = LogisticRegression(max_iter=10000, tol=0.0001, verbose=1, solver='sag', C=search.best_params_["logistic__C"])
        model.fit(X, y)
        y_pred = model.predict(X)
        logger.critical("Done")
        # print("Best parameter (CV score=%0.3f):" % search.best_score_)
        # print(search.best_params_)

        if pf_checkpoint:
            save_obj([model, y_pred, y, groups], pf_checkpoint)
    else:
        model, y_pred, y, groups = load_obj(pf_checkpoint)

    distrs = [get_distribution(y)]
    index = ['training set']
    y_pred = model.predict(X)

    # Plot the PCA spectrum
    for fold_ind, (dev_ind, val_ind) in enumerate(stratified_group_k_fold(X, y, groups, k=5)):
        dev_y, val_y = y[dev_ind], y[val_ind]
        dev_groups, val_groups = groups[dev_ind], groups[val_ind]

        assert len(set(dev_groups) & set(val_groups)) == 0

        distrs.append(get_distribution(dev_y))
        index.append(f'development set - fold {fold_ind}')
        distrs.append(get_distribution(val_y))
        index.append(f'validation set - fold {fold_ind}')

    print(pd.DataFrame(distrs, index=index, columns=[f'Label {l}' for l in range(np.max(y) + 1)]).to_csv())


    df["Target"] = y
    df["Prediction"] = y_pred

    df["Match"] = df["Target"] == df["Prediction"]

    df_stats = df.groupby(["Target", "Sequence Length"], as_index=False).mean()
    sns.scatterplot(df_stats, "Sequence Length", "Match", hue="Target")
    sns.lineplot(df_stats, "Sequence Length", "Match", hue="Target")
    sns.lmplot(df_stats, "Sequence Length", "Match", hue="Target", sns_kwargs={"lowess": True,
                                                                               "scatter_kws": {"s": 0}})





if __name__ == "__main__":
    main(my_env, parsed_args)
