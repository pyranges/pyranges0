import numpy as np
import pandas as pd
import pytest

import pyranges as pr


def test_fdr_enforces_benjamini_hochberg_monotonicity():
    p_vals = pd.Series([0.0039591368855297175, 0.0037600512992788937, 0.0075061166500909205])

    result = pr.stats.fdr(p_vals)

    expected = np.array([0.005938705328294576, 0.005938705328294576, 0.0075061166500909205])
    assert isinstance(result, pd.Series)
    np.testing.assert_allclose(result.to_numpy(), expected)
    assert result.iloc[1] <= result.iloc[0] <= result.iloc[2]


def test_fdr_matches_scipy_benjamini_hochberg():
    scipy_stats = pytest.importorskip("scipy.stats")
    p_vals = np.array([0.01, 0.01, 0.05, 0.9, 0.8, 0.002, 0.003, 0.001])

    result = pr.stats.fdr(p_vals)
    expected = scipy_stats.false_discovery_control(p_vals, method="bh")

    np.testing.assert_allclose(result, expected)


def test_fdr_handles_ties_and_clipping_matching_bh():
    p_vals = np.array([0.01, 0.01, 0.05, 0.9, 0.8])

    result = pr.stats.fdr(p_vals)

    expected = np.array([0.025, 0.025, 0.08333333333333333, 0.9, 0.9])
    np.testing.assert_allclose(result, expected)


def test_fdr_accepts_list_and_preserves_nans_and_empty():
    list_result = pr.stats.fdr([0.01, 0.02, 0.05])
    np.testing.assert_allclose(list_result, [0.03, 0.03, 0.05])

    nan_result = pr.stats.fdr(pd.Series([0.01, np.nan, 0.02], index=["a", "b", "c"]))
    assert isinstance(nan_result, pd.Series)
    assert list(nan_result.index) == ["a", "b", "c"]
    assert np.isnan(nan_result.loc["b"])
    np.testing.assert_allclose(nan_result.loc[["a", "c"]].to_numpy(), [0.02, 0.02])

    empty_result = pr.stats.fdr(np.array([], dtype=float))
    assert empty_result.size == 0
