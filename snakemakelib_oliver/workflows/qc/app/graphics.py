# Copyright (C) 2015 by Per Unneberg
import pandas as pd
from math import log10
from snakemakelib.application import SampleApplication, PlatformUnitApplication
from snakemakelib.io import IOTarget, IOAggregateTarget
from snakemakelib.graphics import scatter, points, tooltips, facet_grid, colorbrewer, mlines, lines
from blaze import Data, append, odo, DataFrame
from snakemakelib_oliver.odo import fastqc
from bokeh.charts import Scatter
from bokeh.plotting import figure, gridplot
