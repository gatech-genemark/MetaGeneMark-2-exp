import logging
from typing import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


from sbsp_general.general import get_value
from sbsp_viz.general import FigureOptions, add_identity, save_figure

logger = logging.getLogger(__name__)


def kdeplot(df, x, y, hue=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None
    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())

    _, ax = plt.subplots()
    y_df = None if y is None else df[y]

    g = sns.kdeplot(df[x], y_df, legend=False, **sns_kwargs)

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    FigureOptions.set_properties_for_axis(ax, figure_options)
    save_figure(figure_options)
    plt.show()

def scatterplot(df, x, y, hue=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None

    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    ax = get_value(kwargs, "ax", None)

    identity = get_value(kwargs, "identity", False)

    if not ax:
        _, ax = plt.subplots()

    g = sns.scatterplot(x=x, y=y, hue=hue, data=df, linewidth=0, **sns_kwargs)

    if identity:
        add_identity(ax, color="r", ls="--")

    FigureOptions.set_properties_for_axis(ax, figure_options)
    legend = get_value(kwargs, "legend", "full")
    legend_loc = get_value(kwargs, "legend_loc", None)
    if hue is not None and legend:
        title = get_value(kwargs, "legend_title", None)
        if not legend_loc:
            plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), title=title)
        else:
            plt.legend(loc=legend_loc)


    save_figure(figure_options)
    plt.show()


def lineplot(df, x, y, hue=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None

    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    ax = get_value(kwargs, "ax", None)
    show = get_value(kwargs, "show", ax is None)
    legend = get_value(kwargs, "legend", "full")
    legend_loc = get_value(kwargs, "legend_loc", None)
    legend_ncol = get_value(kwargs, "legend_ncol", 1)

    identity = get_value(kwargs, "identity", False)

    if not ax:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    g = sns.lineplot(x=x, y=y, hue=hue, data=df, ax=ax, legend=legend, **sns_kwargs)

    if identity:
        add_identity(ax, color="r", ls="--")

    FigureOptions.set_properties_for_axis(ax, figure_options)
    if hue is not None and legend:
        title = get_value(kwargs, "legend_title", None)
        if not legend_loc:
            plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), title=title, ncol=legend_ncol)
        else:
            plt.legend(loc=legend_loc, ncol=legend_ncol, title=title)
        if title is not None and len(title)  == 0:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles=handles[1:], labels=labels[1:], ncol=legend_ncol)

    if show:
        save_figure(figure_options, fig)
        plt.show()



def catplot(df, x, y, hue=None, kind="box", figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None
    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())

    g = sns.catplot(x=x, y=y, data=df, kind=kind, hue=hue, legend=False, aspect=1.5, **sns_kwargs)

    if kind == "point":
        plt.setp(g.ax.lines, linewidth=1)  # set lw for all lines of g axes
        # plt.setp(g.ax.lines, markersize=0)  # set lw for all lines of g axes

    FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    legend = get_value(kwargs, "legend", "full")
    legend_loc = get_value(kwargs, "legend_loc", None)
    if hue is not None and legend:
        title = get_value(kwargs, "legend_title", None)
        if not legend_loc:
            plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), title=title)
        else:
            plt.legend(loc=legend_loc)

    save_figure(figure_options)
    plt.show()


def lmplot(df, x, y, hue=None, figure_options=None, **kwargs):
    # type: (pd.DataFrame, str, str, Union[str, None], FigureOptions, Dict[str, Any]) -> None

    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    if "aspect" not in sns_kwargs:
        sns_kwargs["aspect"] = 2

    g = sns.lmplot(x=x, y=y, hue=hue, data=df, legend=False, **sns_kwargs)

    FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    legend = get_value(kwargs, "legend", "full")
    legend_loc = get_value(kwargs, "legend_loc", None)
    if hue is not None and legend:
        title = get_value(kwargs, "legend_title", None)
        if not legend_loc:
            g.axes[0][0].legend(loc='center left', bbox_to_anchor=(1.05, 0.5), title=title)
        else:
            g.axes[0][0].legend(loc=legend_loc)

    save_figure(figure_options, fig=g.fig)
    plt.subplots_adjust(right=1)
    plt.show()
    return g


def distplot(df, x, figure_options=None, **kwargs):
    _, ax = plt.subplots()

    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    if "kde" not in sns_kwargs:
        sns_kwargs["kde"] = True

    g = sns.distplot(df[x], bins=50, **sns_kwargs)

    FigureOptions.set_properties_for_axis(g.axes, figure_options)
    save_figure(figure_options)
    plt.show()


def jointplot(df, x, y, hue=None, figure_options=None, **kwargs):
    _, ax = plt.subplots()
    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    # g = sns.lmplot(x=x, y=y, hue=hue, data=df, aspect=2, legend=False, ci=None)
    g = sns.jointplot(x, y, data=df, **sns_kwargs)

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    # FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    save_figure(figure_options)
    plt.show()


def tsplot(df, x, y, hue=None, figure_options=None, **kwargs):
    _, ax = plt.subplots()
    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    # g = sns.lmplot(x=x, y=y, hue=hue, data=df, aspect=2, legend=False, ci=None)
    sns.tsplot(df[y].values, df[x].values, **sns_kwargs)

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    # FigureOptions.set_properties_for_axis(g.axes[0][0], figure_options)
    save_figure(figure_options)
    plt.show()

def barplot(df, x, y, hue, figure_options=None, **kwargs):
    sns_kwargs = get_value(kwargs, "sns_kwargs", dict())
    ax = get_value(kwargs, "ax", None)

    g = sns.barplot(x=x, y=y, data=df, hue=hue,  ax=ax, **sns_kwargs)

    if hue is not None:
        plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))

    FigureOptions.set_properties_for_axis(g, figure_options)
    plt.tight_layout()
    save_figure(figure_options)
    # plt.tight_layout(rect=[-0.3,0,1,1.2])
    plt.show()




