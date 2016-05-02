import re
from blaze import resource, DataFrame
import pandas as pd
from snakemakelib.odo.pandas import annotate_by_uri


@resource.register('.+fastq.summary')
@annotate_by_uri
def resource_fastqc_summary(uri, **kwargs):
    with open(uri):
        data = pd.read_csv(uri, sep=",", index_col=["fileName"])
    return DataFrame(data)
