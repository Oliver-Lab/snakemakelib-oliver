from blaze import DataFrame, odo
from snakemakelib_oliver.odo import fastqc
import pytest

@pytest.fixture(scope="module")
def fastqc_fq_summary_data(tmpdir_factory):
    fn = tmpdir_factory.mktemp('data').join('sample_1_fastqc/summary.txt')
    fn.write("""
PASS	Basic Statistics	DamIDseq_1_1.fastq.gz
PASS	Per base sequence quality	DamIDseq_1_1.fastq.gz
PASS	Per tile sequence quality	DamIDseq_1_1.fastq.gz
PASS	Per sequence quality scores	DamIDseq_1_1.fastq.gz
FAIL	Per base sequence content	DamIDseq_1_1.fastq.gz
FAIL	Per sequence GC content	DamIDseq_1_1.fastq.gz
PASS	Per base N content	DamIDseq_1_1.fastq.gz
PASS	Sequence Length Distribution	DamIDseq_1_1.fastq.gz
WARN	Sequence Duplication Levels	DamIDseq_1_1.fastq.gz
FAIL	Overrepresented sequences	DamIDseq_1_1.fastq.gz
PASS	Adapter Content	DamIDseq_1_1.fastq.gz
FAIL	Kmer Content	DamIDseq_1_1.fastq.gz
""")
    return fn


def test_fastqc_fq(fastqc_fq_summary_data):
    df = odo(str(fastqc_fq_summary_data), DataFrame)
    assert 1 == 1
    print(df)
