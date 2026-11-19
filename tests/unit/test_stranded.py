import numpy as np

import pyranges as pr

np.random.seed(0)


def test_stranded():
    cpg = pr.data.cpg()
    exons = pr.data.exons()

    j = cpg.join(exons)

    assert j.stranded

    j.Strand = "."

    assert not j.stranded

    j.Strand = np.random.choice("+ -".split(), size=len(j))

    assert j.stranded

    for _, df in j:
        assert len(df.Strand.drop_duplicates()) == 1


def test_unstrand():
    exons = pr.data.exons()

    cpg = pr.data.cpg()
    print(exons)
    print(exons.columns)
    x = exons.unstrand()
    print(x.stranded)
    print(x)
    print(x.columns)
    for _, df in x:
        assert not df.index.duplicated().sum()

    # x = pr.PyRanges(x.df.reset_index(drop=True))
    cpg.join(x)


def test_dataframe_round_trip_preserves_stranded_merge():
    ranges = pr.from_string("Chromosome Start End Strand\nchr1 0 10 +\nchr1 5 15 -")
    frame = ranges.df
    assert "." in frame.Strand.cat.categories

    reconstructed = pr.PyRanges(frame)
    merged = reconstructed.merge(strand=True, count=True).df
    assert merged[["Start", "End", "Strand", "Count"]].values.tolist() == [[0, 10, "+", 1], [5, 15, "-", 1]]


def test_dataframe_round_trip_with_unknown_strand_is_unstranded():
    ranges = pr.from_string("Chromosome Start End Strand\nchr1 0 10 +\nchr1 5 15 .")
    reconstructed = pr.PyRanges(ranges.df)
    assert list(reconstructed.Strand) == ["+", "."]
    assert not reconstructed.stranded
