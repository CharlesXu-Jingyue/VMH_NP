#!/usr/bin/env python3
"""Step 6: peptide and receptor expression by VMH compartment (parts 1 and 2).

Compartments are the aligned unified-atlas labels of step 3: VMHDM, VMHC, VMHVL, the
undivided ends of the nucleus (VMH) and the shell (VMHSh, kept separate). "Core" pools the
first four. Two estimates are given for every curated gene and compartment:

  projected   sum over cell types of (share of the compartment's neurons in that type, from
              the spatial data) x (fraction of that type expressing beyond background, from the
              10x v3 data of step 2; zero where step 2 did not call the type expressing). Types are clusters where the cluster has >= 20 cells in
              the 10x group, otherwise the supertype, otherwise the subclass. Available for
              every gene, and separately for female and male 10x cells. Genes that most of
              the compartment's types detect only at the non-neuronal level are called
              "ubiquitous" (glia_level_share), since the test cannot place them.
  measured    for genes on the brain's own spatial panel: fraction of the compartment's neurons
              above a background threshold, the 95th percentile of that gene in non-neuronal
              cells of the same brain. Ties at the cutoff can leave fewer than 5% of reference cells above it.

Outputs (results/tables)
  compartment_composition.csv   neurons per brain, compartment and supertype
  gene_by_compartment.csv       both estimates for every gene, brain and compartment
  gene_by_ap.csv                projected fraction in the anterior, middle and posterior thirds of the core
  sex_differences.csv           female vs male 10x donors, per VMH supertype and gene
and <data_root>/derived/cluster_projection.parquet (per-cluster values used for projection).
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import COMPARTMENTS, CORE, DETECTION_MIN, D, DATASETS, NON_NEURONAL, SCREENED, load_cells  # noqa: E402
from paths import CONFIG, RESULTS  # noqa: E402

T = RESULTS / "tables"


def rank_sum_exact(a, b):
    """Two-sided exact permutation p-value of the Wilcoxon rank-sum statistic, with midranks for ties.

    Enumerates every split of the pooled values (at most a few thousand here), so the result does not
    depend on which approximation a library version picks for small or tied samples."""
    r = rankdata(np.r_[a, b])
    na, n = len(a), len(a) + len(b)
    centre = na * (n + 1) / 2
    sums = np.array([r[list(c)].sum() for c in combinations(range(n), na)])
    return float((np.abs(sums - centre) >= abs(r[:na].sum() - centre) - 1e-9).mean())
MIN_10X = 20


def gene_kinds():
    kind = {}
    for name, k in (("neuropeptide_ligands", "ligand"), ("neuropeptide_receptors", "receptor"), ("context_genes", "context")):
        for g in pd.read_csv(CONFIG / f"{name}.csv").gene:
            kind.setdefault(g, k)
    return kind


def type_matrices():
    """{(sex, level): DataFrame type x gene} for excess fraction and mean log2(CPM+1), 10x v3."""
    e = pd.read_parquet(T / "expr_by_type.parquet")
    for c in ("level", "type", "chemistry", "sex", "gene", "call"):
        e[c] = e[c].astype(str)
    e = e[(e.chemistry == "10Xv3") & (e.n_cells >= MIN_10X)]
    out = {}
    # a type contributes only where step 2 called it expressing; sub-threshold residue is set to zero
    e["gated"] = e.excess_frac.where(e.call.isin(["broad", "subset", "sparse"]), 0.0)
    e["glia"] = (e.call == "glia_level").astype(float)
    for (sex, level), d in e.groupby(["sex", "level"]):
        out[(sex, level)] = (d.pivot(index="type", columns="gene", values="gated"),
                             d.pivot(index="type", columns="gene", values="mean_log2cpm"),
                             d.pivot(index="type", columns="gene", values="glia"))
    return out


def resolve(mats, sex, cluster, supertype, subclass):
    """Expression vectors for one spatial cell type, with the fallback chain; returns (frac, mean, source)."""
    chain = [(sex, "cluster", cluster), (sex, "supertype", supertype)]
    if sex != "all":
        chain += [("all", "cluster", cluster), ("all", "supertype", supertype)]
    chain += [("all", "subclass", subclass)]
    for s, level, key in chain:
        f, m, gl = mats.get((s, level), (None, None, None))
        if f is not None and key in f.index:
            return f.loc[key], m.loc[key], gl.loc[key], f"{level}:{s}"
    return None, None, None, "none"


def labelled_neurons(ds, shifts):
    c = load_cells(ds)
    lab = pd.read_parquet(D / "atlas_labels.parquet")
    c = c.merge(lab[lab.dataset == ds][["cell_label", "kim_aligned"]], on="cell_label", how="left")
    c["compartment"] = c.kim_aligned.where(c.kim_aligned.isin(COMPARTMENTS), "outside")
    c["ap_aligned"] = c.x_ccf + shifts.loc[ds, "shift_ap_mm"]
    return c


def projection(neurons, mats, genes):
    """Projected fraction and mean per group of neurons, for sex in all/F/M."""
    g = neurons.groupby(["cluster", "supertype", "subclass"]).size().rename("n").reset_index()
    res, cover = {}, 0.0
    for sex in ("all", "F", "M"):
        num_f, num_m, num_g, den = np.zeros(len(genes)), np.zeros(len(genes)), np.zeros(len(genes)), 0.0
        for r in g.itertuples():
            f, m, gl, src = resolve(mats, sex, r.cluster, r.supertype, r.subclass)
            if f is None:
                continue
            num_f += r.n * f.reindex(genes).fillna(0).values
            num_m += r.n * m.reindex(genes).fillna(0).values
            num_g += r.n * gl.reindex(genes).fillna(0).values
            den += r.n
            if sex == "all" and src.startswith("cluster"):
                cover += r.n
        res[sex] = (num_f / max(den, 1), num_m / max(den, 1), num_g / max(den, 1))
    return res, cover / max(g.n.sum(), 1)


def main():
    T.mkdir(parents=True, exist_ok=True)
    kind = gene_kinds()
    genes = list(kind)
    mats = type_matrices()
    shifts = pd.read_csv(T / "atlas_alignment.csv").set_index("dataset")
    comp_rows, gene_rows, ap_rows, floor_rows, proj_types = [], [], [], [], {}
    for ds, meta in DATASETS.items():
        c = labelled_neurons(ds, shifts)
        z = np.load(D / f"{meta['tag']}_expr.npz", allow_pickle=True)
        X, panel = z["X"], list(z["genes"])
        assert (z["cell_label"] == c.cell_label.values).all()
        nn = c["class"].str.contains(NON_NEURONAL).values
        thr = {g: float(np.percentile(X[nn, panel.index(g)].astype(float), 95)) for g in genes if g in panel}
        n = c[c.neuron]
        Xn = X[c.neuron.values]
        groups = {k: (n.compartment == k).values for k in COMPARTMENTS}
        groups["core"] = n.compartment.isin(CORE).values
        # Blank probes target no transcript. Put through the same rule as a gene (own 95th percentile of
        # non-neuronal cells, then the fraction of a compartment's neurons above it), they show how large a
        # measured fraction arises with nothing to detect. Only the Allen panel carries them.
        blanks = [j for j, g in enumerate(panel) if g.lower().startswith("blank")]
        floor = {}
        if blanks:
            bthr = np.percentile(X[nn][:, blanks].astype(float), 95, axis=0)
            for comp, m in groups.items():
                if m.any():
                    f = (Xn[m][:, blanks].astype(float) > bthr).mean(axis=0)
                    floor[comp] = float(np.percentile(f, 95))
                    floor_rows.append(dict(dataset=ds, compartment=comp, n_neurons=int(m.sum()), n_blank_probes=len(blanks),
                                           n_zero_cutoff=int((bthr == 0).sum()), blank_frac_median=float(np.median(f)),
                                           blank_frac_p95=floor[comp], blank_frac_max=float(f.max()),
                                           n_blanks_reaching_detection_min=int((f >= DETECTION_MIN).sum())))
        for comp, m in groups.items():
            sub = n[m]
            if not len(sub):
                continue
            vc = sub.supertype.astype(str).value_counts()
            vc = vc.rename_axis("st").reset_index(name="k").sort_values(["k", "st"], ascending=[False, True], kind="stable")   # ties by name
            for st, k in zip(vc.st, vc.k):
                comp_rows.append(dict(dataset=ds, compartment=comp, supertype=st, n_neurons=int(k), share=k / len(sub)))
            res, cover = projection(sub, mats, genes)
            for j, g in enumerate(genes):
                row = dict(dataset=ds, sex_of_brain=n.donor_sex.iloc[0], compartment=comp, gene=g, kind=kind[g], n_neurons=len(sub),
                           proj_frac=res["all"][0][j], proj_frac_F=res["F"][0][j], proj_frac_M=res["M"][0][j],
                           proj_mean_log2cpm=res["all"][1][j], glia_level_share=res["all"][2][j],
                           cluster_level_coverage=cover, measured=g in panel)
                if g in panel:
                    v = Xn[m, panel.index(g)]
                    row.update(meas_frac=float((v.astype(float) > thr[g]).mean()), meas_mean=float(v.mean()), meas_threshold=thr[g],
                               meas_blank_floor=floor.get(comp, np.nan))
                gene_rows.append(row)
        # anterior / middle / posterior thirds of the core, equal neuron counts
        core = n[groups["core"] & n.ap_aligned.notna().values]
        q = core.ap_aligned.quantile([1 / 3, 2 / 3]).values
        third = np.select([core.ap_aligned < q[0], core.ap_aligned < q[1]], ["anterior", "middle"], "posterior")
        for name in ("anterior", "middle", "posterior"):
            sub = core[third == name]
            res, _ = projection(sub, mats, genes)
            for j, g in enumerate(genes):
                ap_rows.append(dict(dataset=ds, third=name, ap_min_mm=sub.ap_aligned.min(), ap_max_mm=sub.ap_aligned.max(),
                                    gene=g, kind=kind[g], n_neurons=len(sub), proj_frac=res["all"][0][j]))
        # per-type vectors actually used, for the explorer
        for r in n.groupby(["cluster", "supertype", "subclass"]).size().reset_index().itertuples():
            if r.cluster not in proj_types:
                f, mlog, _, src = resolve(mats, "all", r.cluster, r.supertype, r.subclass)
                if f is not None:
                    proj_types[r.cluster] = (f.reindex(genes).fillna(0).values, mlog.reindex(genes).fillna(0).values, src)
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(T / "compartment_composition.csv", index=False, float_format="%.4g")
    gb = pd.DataFrame(gene_rows)
    # "ubiquitous": most neurons are of types that detect the gene but not above the non-neuronal level
    gb["call"] = np.select([gb.proj_frac >= 0.5, gb.proj_frac >= 0.2, gb.proj_frac >= 0.1, gb.glia_level_share >= 0.5],
                           ["broad", "subset", "sparse", "ubiquitous"], "low")
    gb["projected_hit"] = gb.proj_frac >= DETECTION_MIN
    gb["measured_hit"] = (gb.meas_frac >= DETECTION_MIN).astype("boolean").where(gb.measured)
    measured_hit = gb.measured_hit.fillna(False).to_numpy(dtype=bool)
    # A measured fraction stands on its own only where blank probes exist and it exceeds what they reach in
    # that compartment. Without that it still counts as agreement when the projection reaches 10% as well.
    supported = measured_hit & (gb.meas_frac > gb.meas_blank_floor).to_numpy(dtype=bool)
    gb["measured_supported"] = pd.Series(supported, index=gb.index).astype("boolean").where(gb.measured)
    gb["evidence"] = np.select(
        [~gb.measured, gb.projected_hit & measured_hit, supported,
         measured_hit & gb.meas_blank_floor.notna(), measured_hit, gb.projected_hit],
        ["off panel", "both above", "measured only", "measured only, within blank noise",
         "measured only, no blank control", "projected only"], default="both below")
    gb["counts"] = gb.compartment.isin(SCREENED) & (gb.projected_hit | supported)
    pd.DataFrame(floor_rows).to_csv(T / "blank_probe_floor.csv", index=False, float_format="%.4g")
    gb.to_csv(T / "gene_by_compartment.csv", index=False, float_format="%.4g")
    pd.DataFrame(ap_rows).to_csv(T / "gene_by_ap.csv", index=False, float_format="%.4g")
    cols = {"cluster": list(proj_types), "source": [v[2] for v in proj_types.values()]}
    for j, g in enumerate(genes):
        cols[f"frac:{g}"] = [v[0][j] for v in proj_types.values()]
        cols[f"mean:{g}"] = [v[1][j] for v in proj_types.values()]
    pd.DataFrame(cols).to_parquet(D / "cluster_projection.parquet", index=False)

    # agreement between the two estimates where both exist
    v = gb[gb.measured & gb.compartment.isin(COMPARTMENTS) & gb.kind.isin(["ligand", "receptor"])]
    for ds, d in v.groupby("dataset"):
        rho = d[["proj_frac", "meas_frac"]].corr(method="spearman").iloc[0, 1]
        print(f"{ds}: {d.gene.nunique()} measured genes; Spearman(projected, measured) across gene x compartment = {rho:.2f}")

    sex_differences(kind)


def sex_differences(kind):
    """Female vs male 10x v3 donors for each VMH supertype: exact rank-sum on per-donor detection fractions."""
    cells = pd.read_parquet(D / "tenx_cells.parquet").reset_index(drop=True)
    X = sp.load_npz(D / "tenx_expr.npz").tocsc()
    genes = pd.read_csv(D / "tenx_genes.csv").gene.tolist()
    e = pd.read_parquet(T / "expr_by_type.parquet")
    for c in ("level", "type", "chemistry", "sex", "gene", "call"):
        e[c] = e[c].astype(str)
    e = e[(e.level == "supertype") & (e.chemistry == "10Xv3") & e.type.str.contains(" VMH ")]
    piv = {s: e[e.sex == s].set_index(["type", "gene"]) for s in ("F", "M")}
    v3 = cells[(cells.library_method == "10Xv3") & cells.supertype.str.contains(" VMH ", na=False)]
    rows = []
    for st, d in v3.groupby("supertype"):
        # Keep the library coverage criterion, then pool eligible cells by biological donor.
        eligible = d[d.groupby("library_label").cell_label.transform("size") >= 15]
        donors = eligible.groupby(["donor_sex", "donor_label"]).groups
        if any(sum(sex == key[0] for key in donors) < 3 for sex in ("F", "M")):
            continue
        libraries = eligible.groupby("donor_sex").library_label.nunique()
        tested_cells = eligible.groupby("donor_sex").size()
        for j, g in enumerate(genes):
            if kind.get(g) not in ("ligand", "receptor") or (st, g) not in piv["F"].index or (st, g) not in piv["M"].index:
                continue
            f, m = piv["F"].loc[(st, g)], piv["M"].loc[(st, g)]
            if max(f.excess_frac, m.excess_frac) < 0.1:
                continue
            col = X[:, j]
            # integer count over integer size: equal fractions are exactly equal, so ties do not depend on the platform
            fractions = {key: int((col[idx] > 0).sum()) / len(idx) for key, idx in donors.items()}
            a = [v for (sex, _), v in fractions.items() if sex == "F"]
            b = [v for (sex, _), v in fractions.items() if sex == "M"]
            p = rank_sum_exact(a, b)
            rows.append(dict(supertype=st, gene=g, kind=kind[g], n_donors_F=len(a), n_donors_M=len(b),
                             n_libraries_F=libraries["F"], n_libraries_M=libraries["M"],
                             tested_cells_F=tested_cells["F"], tested_cells_M=tested_cells["M"],
                             n_cells_F=int(f.n_cells), n_cells_M=int(m.n_cells), excess_frac_F=float(f.excess_frac),
                             excess_frac_M=float(m.excess_frac), difference_F_minus_M=float(f.excess_frac - m.excess_frac),
                             log2cpm_F=float(f.mean_log2cpm), log2cpm_M=float(m.mean_log2cpm), p_value=p))
    s = pd.DataFrame(rows)
    # Benjamini-Hochberg across all tests
    order = np.argsort(s.p_value.values, kind="stable")
    ranked = s.p_value.values[order] * len(s) / (np.arange(len(s)) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    s.loc[s.index[order], "q_value"] = np.minimum(q, 1.0)
    s = s.sort_values(["q_value", "p_value", "supertype", "gene"], kind="stable")
    s.to_csv(T / "sex_differences.csv", index=False, float_format="%.4g")
    hit = s[(s.q_value < 0.05) & (s.difference_F_minus_M.abs() >= 0.15)]
    print(f"sex: {len(s)} tests over {s.supertype.nunique()} supertypes; {len(hit)} with q < 0.05 and |difference| >= 0.15")


if __name__ == "__main__":
    main()
