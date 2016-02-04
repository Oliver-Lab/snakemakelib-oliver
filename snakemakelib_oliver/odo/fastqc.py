from blaze import resource
import pandas as pd
from snakemakelib.odo.pandas import annotate_by_uri

@resource.register('.+\_fastqc/summary.txt')
@annotate_by_uri
def resource_fastqc_summary(uri, **kwargs):
    with open(uri):
        data = pd.read_csv(uri, sep="\t", header=None, comment="#",
                           names = ["flag", "statistic", "sample"],
                           index_col=["sample"])
    return data
