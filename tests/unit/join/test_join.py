from pandas.testing import assert_series_equal

import pyranges as pr


# with slack
def test_join_with_slack():
    gr1 = pr.PyRanges(chromosomes="chr1", starts=[0], ends=[10], strands="+")
    gr2 = pr.PyRanges(chromosomes="chr1", starts=[15], ends=[20], strands="+")

    result = gr1.join(gr2, slack=10)
    df = result.df
    print(df)
    assert not df.empty


def test_join_without_reordering():
    f1 = pr.from_dict(
        {
            "Chromosome": ["chr1", "chr1", "chr1"],
            "Start": [3, 8, 5],
            "End": [6, 9, 7],
            "Name": ["interval1", "interval3", "interval2"],
        }
    )
    f2 = pr.from_dict({"Chromosome": ["chr1", "chr1"], "Start": [1, 6], "End": [2, 7], "Name": ["a", "b"]})

    lj = f1.join(f2, how="left", preserve_order=True)
    assert_series_equal(lj.Name, f1.Name)

    rj = f1.join(f2, how="right", preserve_order=True)
    assert_series_equal(rj.Name_b, f2.Name, check_names=False)

    oj = f1.join(f2, how="outer", preserve_order=True)
    assert list(oj.Name) == ["interval1", "interval3", "interval2", "-1"]
    assert list(oj.Name_b) == ["-1", "-1", "b", "a"]


def test_join_right_and_outer_keep_chromosomes_missing_in_self():
    f1 = pr.from_dict({"Chromosome": ["chr1", "chr2"], "Start": [1, 5], "End": [4, 9], "Name": ["a1", "a2"]})
    f2 = pr.from_dict({"Chromosome": ["chr1", "chr3"], "Start": [2, 7], "End": [3, 8], "Name": ["b1", "b3"]})

    right = f1.join(f2, how="right").df
    assert list(right.Chromosome) == ["chr1", "chr3"]
    assert list(right.Name) == ["a1", "-1"]
    assert list(right.Name_b) == ["b1", "b3"]
    assert list(right.Start) == [1, -1]
    assert list(right.Start_b) == [2, 7]

    outer = f1.join(f2, how="outer").df
    assert list(outer.Chromosome) == ["chr1", "chr2", "chr3"]
    assert list(outer.Name) == ["a1", "a2", "-1"]
    assert list(outer.Name_b) == ["b1", "-1", "b3"]


def test_join_outer_without_common_chromosomes():
    f1 = pr.from_dict({"Chromosome": ["chr2"], "Start": [1], "End": [4], "Name": ["a1"]})
    f2 = pr.from_dict({"Chromosome": ["chr1"], "Start": [2], "End": [3], "Name": ["b1"]})

    outer = f1.join(f2, how="outer").df
    assert list(zip(outer.Chromosome, outer.Name, outer.Name_b)) == [("chr1", "-1", "b1"), ("chr2", "a1", "-1")]


def test_join_right_and_outer_stranded():
    f1 = pr.from_dict(
        {"Chromosome": ["chr1"] * 2, "Start": [1, 1], "End": [4, 4], "Strand": ["+", "-"], "Name": ["a1", "a2"]}
    )
    f2 = pr.from_dict(
        {"Chromosome": ["chr1"] * 2, "Start": [2, 70], "End": [3, 80], "Strand": ["+", "+"], "Name": ["b1", "b2"]}
    )

    def names(**kwargs):
        df = f1.join(f2, **kwargs).df
        return sorted(zip(df.Name, df.Name_b))

    unstranded = [("-1", "b2"), ("a1", "b1"), ("a2", "b1")]
    assert names(how="right") == unstranded
    assert names(how="outer") == unstranded

    assert names(how="right", strandedness="same") == [("-1", "b2"), ("a1", "b1")]
    assert names(how="outer", strandedness="same") == [("-1", "b2"), ("a1", "b1"), ("a2", "-1")]

    assert names(how="right", strandedness="opposite") == [("-1", "b2"), ("a2", "b1")]
    assert names(how="outer", strandedness="opposite") == [("-1", "b2"), ("a1", "-1"), ("a2", "b1")]
