#!/usr/bin/env python3
"""Step 5 (interim): readable type-level summaries that do not depend on any boundary.

results/tables/vmh_supertype_calls.csv   one row per curated gene: for each of the nine VMH
    supertypes the percentage of cells expressing beyond background (10x v3, both sexes) and
    the call; plus the same percentage in females and males pooled over VMH supertypes
results/tables/vmh_type_counts.csv       neurons of each VMH supertype in each spatial dataset,
    and how many of them carry the CCF "VMH" label
"""
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import DATASETS, load_cells  # noqa: E402
from paths import CONFIG, RESULTS  # noqa: E402

T = RESULTS / "tables"
CODE = {"broad": "B", "subset": "S", "sparse": "s", "glia_level": "g", "background": ".", "too_few_cells": "?"}


def main():
    e = pd.read_parquet(T / "expr_by_type.parquet")
    for c in ("level", "type", "chemistry", "sex", "gene", "call"):
        e[c] = e[c].astype(str)
    e = e[(e.level == "supertype") & e.type.str.contains(" VMH ") & (e.chemistry == "10Xv3")]
    kind = {}
    for name, k in (("neuropeptide_ligands", "ligand"), ("neuropeptide_receptors", "receptor"), ("context_genes", "context")):
        for g in pd.read_csv(CONFIG / f"{name}.csv").gene:
            kind.setdefault(g, k)
    a = e[e.sex == "all"]
    pct = (a.pivot(index="gene", columns="type", values="excess_frac") * 100).round().astype(int)
    call = a.pivot(index="gene", columns="type", values="call")
    out = pd.DataFrame({"kind": pd.Series(kind)})
    out["n_supertypes_expressing"] = call.isin(["broad", "subset", "sparse"]).sum(axis=1)
    for t in pct.columns:
        out[t[5:].replace(" Glut", "")] = pct[t].astype(str) + call[t].map(CODE)
    for sex in ("F", "M"):                       # cell-weighted over the nine supertypes
        s = e[e.sex == sex]
        w = s.assign(x=s.excess_frac * s.n_cells).groupby("gene")[["x", "n_cells"]].sum()
        out[f"pct_{sex}_all_VMH_types"] = (100 * w.x / w.n_cells).round(1)
    out = out.sort_values(["kind", "n_supertypes_expressing"], ascending=[True, False])
    out.index.name = "gene"
    out.to_csv(T / "vmh_supertype_calls.csv")
    rows = []
    for ds in DATASETS:
        c = load_cells(ds)
        n = c[c.neuron & c.vmh_type]
        g = n.groupby("supertype").agg(n_neurons=("cell_label", "size"), n_with_ccf_vmh_label=("in_ccf_vmh", "sum")).reset_index()
        g.insert(0, "dataset", ds)
        rows.append(g)
    pd.concat(rows).to_csv(T / "vmh_type_counts.csv", index=False)
    print(out.groupby("kind").n_supertypes_expressing.apply(lambda x: f"{(x > 0).sum()} of {len(x)} genes expressed in at least one VMH supertype").to_string())


if __name__ == "__main__":
    main()
