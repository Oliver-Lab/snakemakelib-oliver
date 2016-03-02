""" basic post processing hooks """

from snakemakelib_oliver.odo import bowtie2, htseq

__all__ = ['dge_htseq_post_processing_hook', 'dge_bowtie2_post_processing_hook']


def dge_bowtie2_post_processing_hook(df, **kwargs):
    df_wide = df.reset_index().pivot_table(values=['counts'], index=["SM", "PU"], columns="statistic")
    df_wide.columns = df_wide.columns.droplevel()
    df_wide['Percent uniquely aligned'] = 100 * df_wide['Number Uniquely Aligned'] / df_wide['Number of Reads']
    return df_wide


def dge_htseq_post_processing_hook(df, **kwargs):
    df_wide = df.reset_index().pivot_table(values='count', index="FBgn", columns="SM")
    #df_wide.columns = df_wide.columns.droplevel()
    return df_wide
