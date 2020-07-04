# Author: Karl Gemayel
# Created: 6/22/20, 11:26 AM

import logging
import math
from typing import *


from typing import *
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.axes
from pandas.plotting import scatter_matrix
from mg_general.general import get_value
import seaborn as sns

log = logging.getLogger(__name__)


def add_identity(axes, *line_args, **line_kwargs):
    # type: (matplotlib.axes.Axes, List[str], Dict[str, Any]) -> matplotlib.axes.Axes
    identity, = axes.plot([], [], *line_args, **line_kwargs)

    def callback(l_axes):
        low_x, high_x = l_axes.get_xlim()
        low_y, high_y = l_axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])

    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

def save_figure(figure_options, fig=None):
    # type: (FigureOptions, plt.Figure) -> None
    if figure_options is not None and figure_options.save_fig is not None:
        if fig:
            # fig.tight_layout()
            plt.savefig(figure_options.save_fig) #, bbox_index="tight")
        else:
            plt.savefig(figure_options.save_fig , bbox_index="tight")


class FigureOptions:

    def __init__(self, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, save_fig=None, **kwargs):
        # type: (str, str, str, Tuple(int, int), Tuple(int, int), str, Dict[str, Any]) -> None
        self.balanced = get_value(kwargs, "balanced", False)
        self.legend_title = get_value(kwargs, "legend_title", None)

        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xlim = xlim
        self.ylim = ylim
        self.save_fig = save_fig

    @staticmethod
    def set_properties_for_axis(axis, figure_options):
        # type: (matplotlib.axes.Axes, FigureOptions) -> None

        if figure_options is None:
            return

        if figure_options.title:
            axis.set_title(figure_options.title)
        if figure_options.xlabel:
            axis.set_xlabel(figure_options.xlabel)
        if figure_options.ylabel:
            axis.set_ylabel(figure_options.ylabel)

        if figure_options.xlim:
            axis.set_xlim(*figure_options.xlim)
        if figure_options.ylim:
            axis.set_ylim(*figure_options.ylim)

        if figure_options.balanced:
            min_xy = min(axis.get_xlim()[0], axis.get_ylim()[0])
            max_xy = max(axis.get_xlim()[1], axis.get_ylim()[1])

            axis.set_xlim(min_xy, max_xy)
            axis.set_ylim(min_xy, max_xy)

        if figure_options.legend_title:
            axis.legend().texts[0].set_text(figure_options.legend_title)



def plot_scatter_for_dataframe_columns(df, column_names, figure_options=None, **kwargs):
    # type: (pd.DataFrame, list[str], FigureOptions, Dict[basestring, Any]) -> None

    if len(column_names) != 2:
        raise ValueError("Scatter plot can only be done for 2 columns at a time")

    fig, ax = plt.subplots()

    plot_identity = get_value(kwargs, "plot_identity", False)
    color_by = get_value(kwargs, "color_by_value", None)
    hold = get_value(kwargs, "hold", False)

    column_x = column_names[0]
    column_y = column_names[1]

    plot_kws = {"s": 10, "linewidth": 0, "alpha": 0.3}
    if color_by is not None:
        ax = sns.scatterplot(df[column_x], df[column_y], hue=df[color_by], **plot_kws)
    else:
        ax = sns.scatterplot(df[column_x], df[column_y], **plot_kws)

    # if color_by_value:
    #     fig, ax = plt.subplots()
    #     group = df.groupby(color_by_value)
    #     # df.groupby(color_by_value).plot(column_names[0], column_names[1], s=2, ax=ax, kind="scatter", legend=True)  # type: matplotlib.axes.Axes
    #     for name, df_group in group:
    #         ax.plot(df_group[column_names[0]], df_group[column_names[1]], label=name, linestyle="", marker="o", ms=3)  # type: matplotlib.axes.Axes
    #     ax.legend()
    # else:
    #     ax = df.plot.scatter(column_names[0], column_names[1], s=2)  # type: matplotlib.axes.Axes

    if plot_identity:
        add_identity(ax, color="r", ls="--")

    FigureOptions.set_properties_for_axis(ax, figure_options)

    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig)

    plt.show()

def plot_scatter_matrix_for_dataframe_columns(df_data, column_names, figure_options=None):
    # type: (pd.DataFrame, Iterable[str], FigureOptions) -> None

    df_features = df_data[column_names]

    sm = scatter_matrix(df_features, diagonal="kde", figsize=(10, 10))
    # Change label rotation

    [s.xaxis.label.set_rotation(45) for s in sm.reshape(-1)]
    [s.yaxis.label.set_rotation(0) for s in sm.reshape(-1)]

    # May need to offset label when rotating to prevent overlap of figure
    [s.get_yaxis().set_label_coords(-0.3, 0.5) for s in sm.reshape(-1)]

    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_index="tight")

    plt.show()


def plot_scatter_with_error_for_dataframe_columns(df, column_names, column_yerr, figure_options=None):
    # type: (pd.DataFrame, List[str], str, FigureOptions) -> None

    fig, ax = plt.subplots()

    ax.errorbar(df[column_names[0]], df[column_names[1]], df[column_yerr], linestyle='None', marker='^')

    FigureOptions.set_properties_for_axis(ax, figure_options)

    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig)

    plt.show()



def plot_hist_from_raw_data(mylist, bins=None, figure_options=None):
    # type: (List, FigureOptions) -> None

    fig, ax = plt.subplots()
    ax.hist(mylist, bins=bins)

    FigureOptions.set_properties_for_axis(ax, figure_options)

    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig)

    plt.show()


def plot_hist_from_list(mylist, figure_options=None):
    # type: (List, FigureOptions) -> None

    fig, ax = plt.subplots()
    ax.hist(mylist, bins=range(min(mylist), max(mylist)))

    FigureOptions.set_properties_for_axis(ax, figure_options)

    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig)

    plt.show()


def jitter(df, column_names, jitter="auto"):
    # type: (pd.DataFrame, List[str], Union[str, float]) -> None

    n = len(df)
    for c in column_names:
        if jitter == "auto":
            max_value = max(df[c])
            min_value = min(df[c])

            jitter_mag = (max_value - min_value) * 0.001
        else:
            jitter_mag = float(jitter)

        df[c] = np.random.uniform(-jitter_mag, +jitter_mag, n) + df[c]

def plot_scatter_matrix(df_data, column_names, color_by, figure_options=None, **kwargs):
    if color_by is not None:
        df_features = df_data[column_names + [color_by]]
    else:
        df_features = df_data[column_names]

    should_jitter = get_value(kwargs, "jitter", False)
    if should_jitter:
        jitter(df_features, column_names)

    fig, ax = plt.subplots()

    sns.set(style="ticks")

    if color_by is not None:
        ax = sns.pairplot(df_features, hue=color_by, plot_kws={"s": 10, "linewidth": 0, "alpha": 0.3})
    else:
        ax = sns.pairplot(df_features, plot_kws={"s": 10})

    # for lh in ax._legend.legendHandles:
    #     lh.set_alpha(1)
    #     lh._sizes = [50]

        # sm = scatter_matrix(df_features, diagonal="kde", figsize=(10, 10))
    # # Change label rotation
    #
    # [s.xaxis.label.set_rotation(45) for s in sm.reshape(-1)]
    # [s.yaxis.label.set_rotation(0) for s in sm.reshape(-1)]
    #
    # # May need to offset label when rotating to prevent overlap of figure
    # [s.get_yaxis().set_label_coords(-0.3, 0.5) for s in sm.reshape(-1)]
    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_index="tight")

    plt.show()

def plot_scatter_multiple_columns(df_data, column_x, list_column_y, list_labels=None, figure_options=None):
    fig, ax = plt.subplots()
    sns.set(style="ticks")

    for column_y, name in zip(list_column_y, list_labels):
        sns.scatterplot(df_data[column_x], df_data[column_y], ax=ax, label=name)

    # ax[0].legend(list_labels)
    FigureOptions.set_properties_for_axis(ax, figure_options)


    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_index="tight")

    plt.show()

def plot_scatter(df_data, column_x, column_y, color_by=None, figure_options=None, hold=False, **kwargs):

    if not hold:
        fig, ax = plt.subplots()

    sns.set(style="ticks")

    if color_by is not None:
        ax = sns.scatterplot(df_data[column_x], df_data[column_y], hue=color_by)
    else:
        ax = sns.scatterplot(df_data[column_x], df_data[column_y])

    FigureOptions.set_properties_for_axis(ax, figure_options)



    plot_identity = get_value(kwargs, "plot_identity", False)
    hv_dp = get_value(kwargs, "hv_dp", None)

    if plot_identity:
        add_identity(ax, color="r", ls="--")

    if hv_dp:
        xmin = ax.get_xlim()[0]
        ymin = ax.get_ylim()[0]

        x_val = hv_dp[0]
        y_val = hv_dp[1]

        ax.hlines(y=y_val, xmin=xmin, xmax=x_val, color='b')
        ax.vlines(x=x_val, ymin=ymin, ymax=y_val, color='b')

    if figure_options is not None and figure_options.save_fig is not None:
        plt.savefig(figure_options.save_fig, bbox_index="tight")

    plt.show()


def square_subplots(num_items):
    num_rows = int(math.sqrt(num_items))
    num_cols = math.ceil(num_items / float(num_rows))
    return num_rows, num_cols