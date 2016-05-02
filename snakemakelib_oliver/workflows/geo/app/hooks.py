""" basic post processing hooks """

from snakemakelib_oliver.odo import geo

__all__ = ['geo_summary_table_post_processing_hook']


def geo_summary_table_post_processing_hook(df, **kwargs):
    return df.drop(['SM', 'PU', 'PlatformUnit'], axis=1)
