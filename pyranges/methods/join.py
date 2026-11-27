import numpy as np
import pandas as pd
from ncls import NCLS  # type: ignore


def _both_indexes(scdf, ocdf, how=False, **kwargs):
    assert (how in "containment first last left".split() + [False, None]) or isinstance(how, int)
    starts = scdf.Start.values
    ends = scdf.End.values
    indexes = scdf.index.values

    it = NCLS(ocdf.Start.values, ocdf.End.values, ocdf.index.values)

    if not how:
        _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)
    elif how == "containment":
        _self_indexes, _other_indexes = it.all_containments_both(starts, ends, indexes)
    elif how == "first":
        _self_indexes, _other_indexes = it.first_overlap_both(starts, ends, indexes)
    elif how == "last":
        _self_indexes, _other_indexes = it.last_overlap_both(starts, ends, indexes)
    elif how == "left":
        _self_indexes, _other_indexes = it.all_overlaps_both(starts, ends, indexes)

        missing_in_s = scdf.index.difference(_self_indexes)
        filler_o = np.ones(len(missing_in_s), dtype=int) * -1

        _self_indexes = np.concatenate([_self_indexes, missing_in_s])
        _other_indexes = np.concatenate([_other_indexes, filler_o])

        # sort indexes if the order should be preserved
        if kwargs["preserve_order"]:
            _sort_indexes = np.argsort(_self_indexes, kind="stable")
            _self_indexes = _self_indexes[_sort_indexes]
            _other_indexes = _other_indexes[_sort_indexes]

    return _self_indexes, _other_indexes


def null_types(h):
    h2 = h.copy()
    for n, d in zip(h, h.dtypes):
        if n in ["Chromosome", "Strand"]:
            continue

        d = str(d)

        # for bools, turn into a category so we can add an
        # "unknown" null value ("-1")
        if "bool" in d:
            h2[n] = h2[n].astype("category")
            d = "category"

        # use numerical -1 for numbers
        # or string "-1" for strings, objects, and cats
        if "int" in d or "float" in d:
            null = -1
        elif "str" in d or d == "object":
            null = "-1"
        elif d == "category":
            tmp_cat = h2[n].copy()
            tmp_cat = tmp_cat.cat.add_categories("-1")
            h2[n] = tmp_cat
            null = "-1"
        else:
            raise Exception("Unknown dtype {} in a column {}".format(d, n))

        h2[n] = null

    return h2


def _both_dfs(scdf, ocdf, how=False, **kwargs):
    _self_indexes, _other_indexes = _both_indexes(scdf, ocdf, how, **kwargs)

    if how == "left":
        oh = null_types(ocdf.head(1))
        oh.index = [-1]

        ocdf = pd.concat([ocdf, oh])

        scdf = scdf.reindex(_self_indexes)
        ocdf = ocdf.reindex(_other_indexes)

        if "Strand" in scdf and "Strand" in ocdf:
            x = ocdf.index.values == -1
            ocdf.loc[x, "Strand"] = scdf[x].Strand.values

    else:
        scdf = scdf.reindex(_self_indexes)
        ocdf = ocdf.reindex(_other_indexes)

    return scdf, ocdf


def _join_columns(scdf, ocdf, **kwargs):
    if not kwargs.get("new_pos"):
        suffix = kwargs.get("suffix", "_b")
    else:
        suffix = kwargs.get("suffixes", "_a _b".split())[1]

    nix = pd.Index(range(len(scdf)))
    scdf.index = nix
    ocdf.index = nix

    ocdf = ocdf.drop("Chromosome", axis=1)

    df = scdf.join(ocdf, rsuffix=suffix)

    if kwargs.get("report_overlap"):
        df["Overlap"] = df[["End", "End" + suffix]].min(axis=1) - df[["Start", "Start" + suffix]].max(axis=1)

    return df


def _write_both(scdf, ocdf, **kwargs):
    if scdf.empty or ocdf.empty:
        if kwargs["how"] == "left" and ocdf.empty:
            ocdf = null_types(kwargs["example_header_other"])
        else:
            return None

    scdf, ocdf = _both_dfs(scdf, ocdf, **kwargs)

    return _join_columns(scdf, ocdf, **kwargs)


def _write_missing_in_self(ocdf, **kwargs):
    """Join the rows of other without overlap with null values for the columns of self."""

    header = null_types(kwargs["example_header_self"])
    scdf = header.loc[header.index.repeat(len(ocdf))].copy()
    scdf["Chromosome"] = ocdf.Chromosome.values

    if "Strand" in scdf:
        scdf["Strand"] = ocdf.Strand.values if "Strand" in ocdf else "."

    return _join_columns(scdf, ocdf.copy(), **kwargs)
