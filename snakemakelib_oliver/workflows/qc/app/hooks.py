# basic post processing hooks
import re
import pandas

from snakemakelib.odo import cutadapt
from snakemakelib_oliver.odo import fastqc, bowtie2, htseq

__all__ = ['qc_fastqc_summary_hook',
           'qc_cutadapt_post_processing_hook',
           'qc_bowtie2_post_processing_hook',
           'qc_htseq_post_processing_hook']


def get_base(x):
    return re.sub('\_fastqc.*', '', x)


def qc_fastqc_summary_hook(df, **kwargs):
    df['PU'] = df['PU'].map(get_base)
    df.drop('PlatformUnit', axis=1)
    df_wide = df.reset_index().pivot_table(values=["flag"], index=["SM", "PU"], columns="statistic", aggfunc=lambda x: x)
    df_wide.columns = df_wide.columns.droplevel()
    return df_wide


def qc_cutadapt_post_processing_hook(df, **kwargs):
    df_wide = df.reset_index().pivot_table(values=["value"], index=["SM", "PU", "PlatformUnit"], columns="statistic")
    df_wide.columns = df_wide.columns.droplevel()
    df_wide["Reads percent"] = 100.0 * df_wide["Reads with adapters"] /\
        df_wide["Total reads processed"]
    df = df_wide.stack()
    df.name = "value"
    return df


def qc_bowtie2_post_processing_hook(df, **kwargs):
    df_wide = df.reset_index().pivot_table(values=['counts'], index=["SM", "PU"], columns="statistic")
    df_wide.columns = df_wide.columns.droplevel()
    df_wide['Percent uniquely aligned'] = 100 * df_wide['Number Uniquely Aligned'] / df_wide['Number of Reads']
    return df_wide


def qc_htseq_post_processing_hook(df, **kwargs):
    df_wide = df.reset_index().pivot_table(values=['count'], index=["FBgn"], columns="SM")
    #df_wide.columns = df_wide.columns.droplevel()
    return df
