import re
from blaze import resource, DataFrame
import pandas as pd
from snakemakelib.odo.pandas import annotate_by_uri

@resource.register('.+\_fastqc/summary.txt')
@annotate_by_uri
def resource_fastqc_summary(uri, **kwargs):
    with open(uri):
        data = pd.read_csv(uri, sep="\t", header=None, comment="#",
                           names=["flag", "statistic", "sample"],
                           index_col=["sample"])
    return data


def _data_reader(uri):
    with open(uri) as fh:
        data = "".join(fh)

    # Split sections
    sections = re.split('\n>>', data)
    clean = [x for x in sections if x.strip() != 'END_MODULE' and not x.startswith('##')]

    # Build dataframe for each section and add to dictionary
    dd = dict()
    for section in clean:
        rows = section.split('\n')
        name = rows[0].split('\t')[0]

        # Header rows start with "#"
        if rows[1].startswith('#'):

            # Drop Total deduplicated
            if rows[1].startswith('#Total Deduplicated Percentage'):
                header = rows[2].lstrip('#').split('\t')
                data = [x.split('\t') for x in rows[3:]]
            else:
                header = rows[1].lstrip('#').split('\t')
                data = [x.split('\t') for x in rows[2:]]
            index = header[0]
        else:
            header = None
            data = [x.split('\t') for x in rows[1:]]
            index = None

        # Make dataframe
        df = DataFrame(data)
        df.columns = header
        if index is not None:
            df.set_index(index, inplace=True)
        dd[name] = df

    return dd


@resource.register('.+\_fastqc/fastqc_data.txt')
@annotate_by_uri
def resource_fastqc_data(uri, key='Basic Statistics', **kwargs):
    data = _data_reader(uri)
    return data[key]

