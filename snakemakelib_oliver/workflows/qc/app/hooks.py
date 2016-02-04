# basic post processing hooks
import re
import pandas
from snakemakelib_oliver.odo import fastqc

__all__ = ['qc_fastqc_summary_hook']

def get_base(x):
    return re.sub('\_fastqc.*', '', x)

def qc_fastqc_summary_hook(df, **kwargs):
    df['PU'] = df['PU'].map(get_base)
    df.drop('PlatformUnit', axis=1)
    df_wide = df.reset_index().pivot_table(values=["flag"], index=["SM", "PU"], columns="statistic", aggfunc=lambda x: x)
    df_wide.columns = df_wide.columns.droplevel()
    return df_wide


