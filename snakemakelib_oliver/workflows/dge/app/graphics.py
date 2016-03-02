# Copyright (C) 2015 by Per Unneberg
import pandas as pd
from math import log10
from snakemakelib.application import SampleApplication, PlatformUnitApplication
from snakemakelib.io import IOTarget, IOAggregateTarget
from snakemakelib.graphics import scatter, points, tooltips, facet_grid, colorbrewer, mlines, lines
from blaze import Data, append, odo, DataFrame
from snakemakelib_oliver.odo import fastqc
from snakemakelib.odo import cutadapt
from bokeh.charts import Scatter
from bokeh.plotting import figure, gridplot

__all__ = ['qc_cutadapt_plot_metrics',]


DEFAULT_TOOLS = "pan,wheel_zoom,box_zoom,box_select,reset,save,hover,resize"

# Plotting functions
def qc_cutadapt_plot_metrics(df, **kwargs):
    df.set_index(['SM', 'PU', 'PlatformUnit', 'statistic'], inplace=True)
    df.sortlevel(inplace=True)
    df = df.loc[pd.IndexSlice[:, :, :, "Reads percent"], :].reset_index()
    from bokeh.charts import Scatter
    p = Scatter(df, x="PlatformUnit", y="value",
                color="statistic", legend="top_right",
                title="Cutadapt metrics", ylabel="% reads with adapter")
    return p
