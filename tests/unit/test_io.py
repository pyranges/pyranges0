import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import pyranges as pr

ensembl_gtf = "tests/unit/test_data/ensembl.gtf"


def test_read_gtf():
    gr = pr.read_gtf(ensembl_gtf, full=True)

    assert len(gr.columns) == 26

    df = gr.df
    transcript = df.iloc[1]
    assert transcript["tag"] == "basic"

    exon = df[df["exon_id"] == "ENSE00003812156"].iloc[0]
    assert exon["tag"] == "basic"

    gr = pr.read_gtf(ensembl_gtf, full=True, duplicate_attr=True)
    print(gr.columns)
    assert len(gr.columns) == 26

    df = gr.df
    transcript = df.iloc[1]
    assert transcript["tag"] == "basic"

    exon = df[df["exon_id"] == "ENSE00003812156"].iloc[0]
    assert exon["tag"] == "CCDS,basic"
    # assert list(gr.df.columns[:4]) == "Chromosome Start End Strand".split()


def test_read_gtf_multiple_chunks():
    skiprows = pr.readers.skiprows(ensembl_gtf)
    gr = pr.readers.read_gtf_full(ensembl_gtf, chunksize=2, skiprows=skiprows)
    gr2 = pr.readers.read_gtf_full(ensembl_gtf, skiprows=skiprows)

    df = gr.df
    df2 = gr2.df
    # reading in chunks might mix types because reading a chunk of . 0 will make each value a str
    # however, when reading a chunk of 0 0 the type will be int
    df.loc[:, "Frame"] = df["Frame"].astype(str)

    df.Feature = np.array(df["Feature"].astype(str))
    df2.Feature = np.array(df2["Feature"].astype(str))

    assert_frame_equal(df, df2, check_dtype=False)


def test_read_gff3():
    gr = pr.read_gff3("tests/unit/test_data/gencode.gff3")

    assert len(gr.columns) == 26
    # assert list(gr.df.columns[:4]) == "Chromosome Start End Strand".split()


def test_read_bed():
    pr.read_bed("pyranges/example_data/chipseq.bed")


@pytest.mark.parametrize("as_df", [True, False])
@pytest.mark.parametrize("nrows", [1, 2])
def test_read_restricted_gtf_forwards_output_type_and_row_limit(as_df, nrows):
    ranges = pr.read_gtf(ensembl_gtf, full=False, as_df=as_df, nrows=nrows)

    assert isinstance(ranges, pd.DataFrame if as_df else pr.PyRanges)
    frame = ranges if as_df else ranges.df
    assert len(frame) == nrows
    assert list(frame.Start) == [11868] * nrows
    assert list(frame.End) == [14409] * nrows
    assert list(frame.Feature) == ["gene", "transcript"][:nrows]


def test_read_restricted_gtf_with_zero_rows():
    frame = pr.read_gtf(ensembl_gtf, full=False, as_df=True, nrows=0)
    assert isinstance(frame, pd.DataFrame)
    assert frame.empty
    assert {"Chromosome", "Start", "End", "Feature"}.issubset(frame.columns)

    ranges = pr.read_gtf(ensembl_gtf, full=False, nrows=0)
    assert isinstance(ranges, pr.PyRanges)
    assert ranges.empty
