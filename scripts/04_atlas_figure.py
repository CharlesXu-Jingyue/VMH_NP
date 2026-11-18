#!/usr/bin/env python3
"""Step 4: figures of the reference atlas.

results/figures/atlas_alignment_<dataset>.png   aligned, interpolated compartments over the cell
    types, one figure per brain. Grey rings: neurons the same unified atlas
    places in VMH before the shift, at the released coordinates (not the Allen CCF boundary).
results/figures/atlas_interpolation.png         the interpolated atlas itself: at each plate (over the
    published drawing), halfway between plates, and cut along the anterior-posterior axis next to
    the nearest-plate version."""
import importlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import COMPARTMENT_COLOURS as _CC, D, DATASETS, load_cells  # noqa: E402

COMPARTMENT_COLOURS = {k: v for k, v in _CC.items() if k != "outside"}
from paths import RESULTS  # noqa: E402

SUPERTYPE_COLOURS = {
    "0568 VMH Nr5a1 Glut_3": "#3b6fb6", "0569 VMH Nr5a1 Glut_4": "#58c4e0", "0570 VMH Nr5a1 Glut_5": "#1e8a3c",
    "0567 VMH Nr5a1 Glut_2": "#c9b22b", "0571 VMH Nr5a1 Glut_6": "#a8328a", "0566 VMH Nr5a1 Glut_1": "#7a7a7a",
    "0564 VMH Fezf1 Glut_2": "#f08030", "0565 VMH Fezf1 Glut_3": "#d1301c", "0563 VMH Fezf1 Glut_1": "#111111",
}

def main():
    out = RESULTS / "figures"
    out.mkdir(parents=True, exist_ok=True)
    outlines = json.load(open(D / "atlas_outlines.json"))["datasets"]
    lab = pd.read_parquet(D / "atlas_labels.parquet")
    for ds, meta in DATASETS.items():
        c = load_cells(ds).merge(lab[lab.dataset == ds][["cell_label", "kim_raw"]], on="cell_label")
        zs = c.groupby("brain_section_label")[meta["z"]].median().sort_values()
        nv = c[c.neuron & c.vmh_type].groupby("brain_section_label").size()
        secs = [s for s in zs.index if nv.get(s, 0) >= 30]
        ncol = min(len(secs), 4)
        nrow = int(np.ceil(len(secs) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 4.3 * nrow), squeeze=False)
        v = c[c.neuron & c.vmh_type]
        x0, x1, y0, y1 = v.px.quantile(.005) - .25, v.px.quantile(.995) + .25, v.py.quantile(.005) - .25, v.py.quantile(.995) + .25
        for ax, s in zip(axs.ravel(), secs):
            d = c[c.brain_section_label == s]
            n = d[d.neuron]
            ax.scatter(n.px, n.py, s=2, c="#dcdcdc", linewidths=0)
            raw = n[n.kim_raw.isin(COMPARTMENT_COLOURS)]
            ax.scatter(raw.px, raw.py, s=16, facecolors="none", edgecolors="#9a9a9a", linewidths=0.35)
            for st, col in SUPERTYPE_COLOURS.items():
                e = n[n.supertype == st]
                ax.scatter(e.px, e.py, s=7, c=col, linewidths=0)
            for name, lines in outlines[ds].get(s, {}).items():
                for line in lines:
                    p = np.array(line)
                    ax.plot(p[:, 0], p[:, 1], c=COMPARTMENT_COLOURS[name], lw=1.8 if name != "VMHSh" else 1.2,
                            ls="-" if name != "VMHSh" else (0, (3, 2)))
            ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.set_aspect("equal")
            ax.set_title(f"{s}   {int(nv.get(s, 0))} VMH-type neurons", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
        for ax in axs.ravel()[len(secs):]:
            ax.axis("off")
        handles = [plt.Line2D([], [], c=col, lw=2, label=name) for name, col in COMPARTMENT_COLOURS.items()]
        handles += [plt.Line2D([], [], marker="o", ls="", c=col, ms=5, label=st[5:]) for st, col in SUPERTYPE_COLOURS.items()]
        handles += [plt.Line2D([], [], marker="o", ls="", mfc="none", mec="#9a9a9a", ms=6, label="same atlas before the shift (not Allen CCF)")]
        fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8, frameon=False)
        fig.suptitle(f"{meta['label']}: unified atlas compartments (interpolated between plates) after rigid alignment, over VMH supertypes", fontsize=11)
        fig.tight_layout(rect=(0, 0.07, 1, 0.97))
        fig.savefig(out / f"atlas_alignment_{ds}.png", dpi=110)
        plt.close(fig)
        print(f"{ds}: {len(secs)} sections drawn")
    interpolation_figure(out)


def interpolation_figure(out):
    """The interpolated VMH on its own, against the plates it was built from."""
    from matplotlib.colors import ListedColormap, to_rgba
    from skimage import measure
    ra = importlib.import_module("03_reference_atlas")
    V = ra.load_volume()
    A = ra.SmoothAtlas(V)
    order = ["VMHDM", "VMHC", "VMHVL", "VMH", "VMHSh"]
    code = {n: i + 1 for i, n in enumerate(order)}
    fill = ListedColormap([(1, 1, 1, 0)] + [to_rgba(COMPARTMENT_COLOURS[n], 0.30) for n in order])
    solid = ListedColormap([(1, 1, 1, 1)] + [to_rgba(COMPARTMENT_COLOURS[n], 0.85) for n in order])
    r0, r1, q0, q1 = 600, 740, 455, 575                       # left VMH, 10 um pixels (DV rows, ML columns)
    rr, qq = np.meshgrid(np.arange(r0, r1, 0.5), np.arange(q0, q1, 0.5), indexing="ij")

    def smooth_codes(kf, r, q):
        """Compartment code of the interpolated atlas at plate position kf and pixel position (r, q)."""
        xyz = np.c_[np.ravel(((kf - 1) * 10 + ra.AP0 + 0.5) * 0.01), np.ravel((r + 0.5) * 0.01), np.ravel((q + 0.5) * 0.01)]
        M = A.margins(xyz)
        lab = np.where(M.max(0) > 0, np.array([code[n] for n in A.names])[M.argmax(0)], 0)
        return lab.reshape(np.shape(r)), M.reshape((len(A.names),) + np.shape(r))

    def plate_codes(k, rows, cols):
        P = V[k - ra.K0][rows, cols]
        lab = np.zeros(P.shape, dtype=int)
        for i, n in ra.VMH_IDS.items():
            lab[P == i] = code[n]
        return lab

    plates = list(range(55, 67))
    fig = plt.figure(figsize=(21, 11.6))
    gs = fig.add_gridspec(3, 12, height_ratios=[1, 1, 1.45], hspace=0.30, wspace=0.06)
    for j, k in enumerate(plates):
        for row, kf in ((0, float(k)), (1, k + 0.499)):      # just short of halfway, where the end plates change labelling
            ax = fig.add_subplot(gs[row, j])
            ax.set_xticks([]); ax.set_yticks([])
            if row == 1 and j == len(plates) - 1:
                ax.axis("off")
                h = [plt.Line2D([], [], c=COMPARTMENT_COLOURS[n], lw=2.5, label=lab_) for n, lab_ in
                     (("VMHDM", "dorsomedial"), ("VMHC", "central"), ("VMHVL", "ventrolateral"), ("VMH", "undivided end"), ("VMHSh", "shell"))]
                ax.legend(handles=h, loc="center", fontsize=9, frameon=False, title="compartment", title_fontsize=9)
                continue
            if row == 0:
                ax.imshow(plate_codes(k, slice(r0, r1), slice(q0, q1)), cmap=fill, vmin=0, vmax=5, extent=(q0, q1, r1, r0),
                          interpolation="nearest")
            _, M = smooth_codes(np.full(rr.shape, kf), rr, qq)
            for n, m in zip(A.names, M):
                for c in measure.find_contours(m, 0.0):
                    ax.plot(q0 + c[:, 1] * 0.5 + 0.5, r0 + c[:, 0] * 0.5 + 0.5, c=COMPARTMENT_COLOURS[n], lw=1.5)
            ax.set_xlim(q0, q1); ax.set_ylim(r1, r0); ax.set_aspect("equal")
            ax.set_title(f"plate {k}" if row == 0 else f"halfway, {k} to {k + 1}", fontsize=9)
            if j == 0:
                ax.set_ylabel("at the plates: published drawing (fill)\nand interpolated atlas (lines)" if row == 0
                              else "between plates:\ninterpolated atlas", fontsize=9)
    # cuts along the anterior-posterior axis: nearest plate against interpolated
    core = np.argwhere(np.isin(V[60 - ra.K0][r0:r1, q0:q1], [i for i, n in ra.VMH_IDS.items() if n in ra.CORE]))
    rc, qc = r0 + int(np.median(core[:, 0])), q0 + int(np.median(core[:, 1]))
    kfs = np.arange(53.5, 68.01, 0.02)
    near_k = np.clip(np.round(kfs).astype(int), ra.K0, ra.K1)
    cuts = [("horizontal cut", "medial to lateral", np.arange(q1 - 0.5, q0, -0.5), lambda k: plate_codes(k, rc, slice(q0, q1))[::-1],
             lambda kk, u: smooth_codes(kk, np.full(kk.shape, float(rc)), u)[0]),
            ("sagittal cut", "dorsal to ventral", np.arange(r0, r1, 0.5), lambda k: plate_codes(k, slice(r0, r1), qc),
             lambda kk, u: smooth_codes(kk, u, np.full(kk.shape, float(qc)))[0])]
    for ci, (name, ulabel, us, nearest, smooth) in enumerate(cuts):
        kk, uu = np.meshgrid(kfs, us, indexing="ij")
        imgs = (np.stack([np.repeat(nearest(k), 2)[:len(us)] for k in near_k]), smooth(kk, uu))
        for cj, (img, title) in enumerate(zip(imgs, ("nearest plate", "interpolated"))):
            ax = fig.add_subplot(gs[2, (ci * 2 + cj) * 3:(ci * 2 + cj) * 3 + 3])
            ax.imshow(img.T, cmap=solid, vmin=0, vmax=5, aspect="auto", origin="upper", interpolation="nearest",
                      extent=(kfs[0], kfs[-1], len(us) * 0.005, 0))
            for k in plates:
                ax.axvline(k, c="#00000030", lw=0.6)
            ax.set_xticks(plates[::2]); ax.set_xlabel("plate (100 um apart, anterior to posterior)", fontsize=9)
            ax.set_ylabel(f"{ulabel}, mm", fontsize=9)
            ax.set_title(f"{name}: {title}", fontsize=10)
    fig.suptitle("Unified atlas VMH, left side: published plates and the interpolation between them", fontsize=12, y=0.955)
    fig.savefig(out / "atlas_interpolation.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("atlas_interpolation.png drawn")


if __name__ == "__main__":
    main()
