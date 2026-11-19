#!/usr/bin/env python3
"""Step 8: report figures (results/figures).

  compartment_ligands.png, compartment_receptors.png   percent of neurons expressing, by compartment
  supertype_ligands.png, supertype_receptors.png       percent of cells expressing beyond background, by supertype
  supertype_by_compartment.png                          where each VMH supertype sits
  projected_vs_measured.png                             agreement of the two estimates
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import COMPARTMENTS as COMPS, COMPARTMENT_COLOURS, COMPARTMENT_NAMES, DATASETS, VMH_SUPERTYPES  # noqa: E402
from paths import RESULTS  # noqa: E402

T, F = RESULTS / "tables", RESULTS / "figures"
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e0"
SEQ = LinearSegmentedColormap.from_list("blue", ["#f3f6fb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
SHORT = {"VMHDM": "DM", "VMHC": "C", "VMHVL": "VL", "VMHSh": "shell", "VMH": "ends"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": MUTED,
                     "xtick.color": MUTED, "ytick.color": INK, "text.color": INK, "figure.facecolor": SURFACE,
                     "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE})


def heat(ax, M, title, annotate=True):
    """M: DataFrame of fractions (0-1); NaN cells are left empty."""
    ax.imshow(np.ma.masked_invalid(M.values), cmap=SEQ, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(M.shape[1]), M.columns, fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(range(M.shape[0]), M.index, fontsize=8, style="italic")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks(np.arange(-.5, M.shape[1]), minor=True)
    ax.set_yticks(np.arange(-.5, M.shape[0]), minor=True)
    ax.grid(which="minor", color=SURFACE, lw=2)
    ax.tick_params(which="minor", length=0)
    if annotate:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M.values[i, j]
                if np.isfinite(v) and v >= 0.095:
                    ax.text(j, i, f"{100 * v:.0f}", ha="center", va="center", fontsize=7, color="#ffffff" if v >= 0.5 else INK)
    ax.set_title(title, fontsize=8.5, color=INK, pad=24, loc="left")


def compartment_figure(gb, kind, fname):
    g = gb[gb.kind == kind]
    a = g[g.dataset == "allen_m"].pivot(index="gene", columns="compartment", values="proj_frac")
    keep = g.loc[g.compartment.isin(COMPS) & (g.projected_hit | g.measured_hit.fillna(False)), "gene"].unique()
    order = a.loc[keep, "core"].sort_values(ascending=False).index
    panels = []
    for ds, label in (("allen_m", "Allen M"), ("zhuang1_f", "Zhuang 1 F"), ("zhuang2_m", "Zhuang 2 M")):
        for column, method in (("proj_frac", "Projected"), ("meas_frac", "Measured")):
            values = g[g.dataset == ds].pivot(index="gene", columns="compartment", values=column)
            panels.append((f"{method}\n{label}", values.loc[order, COMPS]))
    fig, axs = plt.subplots(1, len(panels), figsize=(15, 1.7 + 0.24 * len(order)), sharey=True)
    for ax, (title, M) in zip(axs, panels):
        heat(ax, M.rename(columns=SHORT), title)
    what = "Neuropeptide precursors" if kind == "ligand" else "Neuropeptide receptors"
    fig.suptitle(f"{what}: transcript detection (%) by VMH compartment", x=0.01, ha="left", fontsize=11)
    fig.text(0.01, 0.005, "Genes reaching 10% by either method in any brain and compartment. Projected: cell-type mix of the compartment\n"
             "x 10x expression of each type. Measured: neurons above the 95th percentile of non-neuronal cells; empty where the\n"
             "gene is not on that panel. Values under 10 are not printed.", fontsize=7, color=MUTED, va="bottom")
    fig.tight_layout(rect=(0, 0.05 + 0.35 / max(len(order), 6), 1, 0.97))
    fig.savefig(F / fname, dpi=150)
    plt.close(fig)
    return list(order)


def supertype_figure(e, kind_of, kind, fname):
    d = e[(e.level == "supertype") & e.type.isin(VMH_SUPERTYPES) & (e.chemistry == "10Xv3") & (e.sex == "all")]
    d = d[d.gene.map(kind_of) == kind]
    val = d.pivot(index="gene", columns="type", values="excess_frac")[VMH_SUPERTYPES]
    call = d.pivot(index="gene", columns="type", values="call")[VMH_SUPERTYPES]
    val = val.where(call.isin(["broad", "subset", "sparse"]))
    keep = val.index[val.notna().any(axis=1)]
    order = val.loc[keep].fillna(0).max(axis=1).sort_values(ascending=False).index
    M = val.loc[order]
    M.columns = [c[9:].replace(" Glut", "") for c in M.columns]
    fig, ax = plt.subplots(figsize=(6.4, 1.5 + 0.24 * len(order)))
    heat(ax, M, "")
    plt.setp(ax.get_xticklabels(), rotation=35, ha="left")
    what = "Neuropeptide precursors" if kind == "ligand" else "Neuropeptide receptors"
    fig.suptitle(f"{what}: percent of cells expressing beyond background, by VMH supertype", x=0.01, ha="left", fontsize=11)
    fig.text(0.01, 0.005, "10x v3 cells, both sexes. Empty: not called expressing in that supertype. Values under 10 are not printed.",
             fontsize=7, color=MUTED, va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(F / fname, dpi=150)
    plt.close(fig)
    return list(order)


def composition_figure():
    t = pd.read_csv(T / "atlas_supertype_by_compartment.csv")
    t = t[t.labels == "aligned"]
    cols = ["VMHDM", "VMHC", "VMHVL", "VMH", "VMHSh", "outside"]
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6), sharey=True)
    for ax, (ds, meta) in zip(axs, DATASETS.items()):
        d = t[t.dataset == ds].set_index("supertype").reindex(VMH_SUPERTYPES)[cols].fillna(0)
        n = d.sum(axis=1)
        share = d.div(n.where(n > 0), axis=0)
        left = np.zeros(len(d))
        for c in cols:
            ax.barh(range(len(d)), share[c].fillna(0), left=left, height=0.62, color=COMPARTMENT_COLOURS[c],
                    edgecolor=SURFACE, linewidth=1.5, label=COMPARTMENT_NAMES[c])
            left += share[c].fillna(0).values
        ax.set_yticks(range(len(d)), [f"{s[9:].replace(' Glut', '')}" for s in d.index])
        for i, k in enumerate(n):
            ax.text(1.02, i, f"{int(k)}", va="center", fontsize=7, color=MUTED)
        ax.set_xlim(0, 1.0)
        ax.set_xticks([0, .5, 1], ["0", "50", "100%"])
        ax.invert_yaxis()
        ax.set_title(meta["label"], fontsize=9, loc="left")
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=6, fontsize=8, frameon=False)
    fig.suptitle("Where the neurons of each VMH supertype sit (aligned atlas labels); neuron counts at right", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    fig.savefig(F / "supertype_by_compartment.png", dpi=150)
    plt.close(fig)


def validation_figure(gb):
    # the five atlas compartments (the three divisions, the undivided ends and the shell), as in step 6
    v = gb[gb.measured & gb.compartment.isin(COMPS) & gb.kind.isin(["ligand", "receptor"])]
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.7), sharex=True, sharey=True)
    for ax, (ds, meta) in zip(axs, DATASETS.items()):
        d = v[v.dataset == ds]
        rho = d[["proj_frac", "meas_frac"]].corr(method="spearman").iloc[0, 1]
        ax.plot([0, 1], [0, 1], color=GRID, lw=1)
        ax.scatter(d.proj_frac, d.meas_frac, s=22, color="#2a78d6", edgecolor=SURFACE, linewidth=0.8)
        ax.set_title(f"{meta['label']}\n{d.gene.nunique()} genes x 5 compartments, Spearman {rho:.2f}", fontsize=9, loc="left")
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
        ax.set_xticks([0, .5, 1], ["0", "50", "100%"]); ax.set_yticks([0, .5, 1], ["0", "50", "100%"])
        ax.grid(color=GRID, lw=1); ax.set_axisbelow(True)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xlabel("projected from cell types (10x)")
    axs[0].set_ylabel("measured in situ")
    fig.suptitle("Projected against measured fraction of neurons expressing, for genes on each spatial panel", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(F / "projected_vs_measured.png", dpi=150)
    plt.close(fig)


def main():
    F.mkdir(parents=True, exist_ok=True)
    gb = pd.read_csv(T / "gene_by_compartment.csv")
    e = pd.read_parquet(T / "expr_by_type.parquet")
    for c in ("level", "type", "chemistry", "sex", "gene", "call"):
        e[c] = e[c].astype(str)
    kind_of = gb.drop_duplicates("gene").set_index("gene").kind
    for kind in ("ligand", "receptor"):
        a = compartment_figure(gb, kind, f"compartment_{kind}s.png")
        b = supertype_figure(e, kind_of, kind, f"supertype_{kind}s.png")
        print(f"{kind}: {len(a)} genes in the compartment figure, {len(b)} in the supertype figure")
    composition_figure()
    validation_figure(gb)


if __name__ == "__main__":
    main()
