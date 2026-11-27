import numpy as np
import pandas as pd

import pyranges as pr


def _example():
    return pd.DataFrame(
        {
            "Chromosome": ["chr1"] * 4 + ["chr2"] * 2,
            "Start": [100, 200, 300, 400, 100, 200],
            "End": [150, 250, 350, 450, 150, 250],
            "Strand": ["+"] * 4 + ["-"] * 2,
            "Gene": ["A"] * 4 + ["B"] * 2,
            "PValue": [1e-5, 0.2, 0.5, 0.9, 0.01, 0.04],
        }
    )


def test_simes_matches_min_n_p_over_rank():
    result = pr.stats.simes(_example(), "Gene", "PValue")

    # Simes (1986): min_i n * p_(i) / i, i.e. R min(p.adjust(p, "BH"))
    np.testing.assert_allclose(result.Simes, [4e-5, 0.02])


def test_simes_keep_position_returns_same_p_values():
    df = _example()

    result = pr.stats.simes(df, "Gene", "PValue", keep_position=True)

    expected = pr.stats.simes(df, "Gene", "PValue")
    assert list(result.Gene) == list(expected.Gene)
    np.testing.assert_allclose(result.Simes, expected.Simes)
    assert list(result.Start) == [100, 100]
    assert list(result.End) == [450, 250]
    assert list(result.Strand) == ["+", "-"]
