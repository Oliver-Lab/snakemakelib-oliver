import re
from blaze import resource, DataFrame
import pandas as pd
from snakemakelib.odo.pandas import annotate_by_uri


@resource.register('.+htseq.counts')
@annotate_by_uri
def resource_fastqc_summary(uri, **kwargs):
    with open(uri):
        data = pd.read_csv(uri, sep="\t", header=None, names=["FBgn", "count"], index_col=["FBgn"])
    return DataFrame(data)
