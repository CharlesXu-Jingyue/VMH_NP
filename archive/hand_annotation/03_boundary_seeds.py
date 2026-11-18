#!/usr/bin/env python3
"""Step 3: data-driven seed boundaries for the VMH and its subdivisions on every spatial
section, as a starting point for manual correction in the boundary editor.

Seeds are computed per section from kernel densities of neurons (sigma = SIGMA_MM):
  VMH     region where neurons of the two VMH subclasses are the local majority
  VMHvl   part of that region where Esr1-rich VMH clusters are the local majority
          (Esr1 status of a cluster comes from the 10x data, step 2)
  VMHdm / VMHc   the rest, split at the midpoint of its extent along the axis pointing
          away from the VMHvl seed. This split is a geometric placeholder only.

Output: results/boundaries/seed_polygons.json and a preview figure.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from skimage import measure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import RESULTS, data_root  # noqa: E402

D = data_root() / "derived"
OUT = RESULTS / "boundaries"
PIX, SIGMA_MM = 0.02, 0.06
NON_NEURONAL = "Astro|Oligo|OPC|Vascular|Immune|OEC"
VMH_SUBCLASSES = ("128 VMH Fezf1 Glut", "129 VMH Nr5a1 Glut")
ESR1_RICH = 0.5        # cluster-level excess fraction of Esr1 in 10x
DATASETS = {
    "allen_m": dict(tag="merfish", label="Allen MERFISH, male", xy=("x_reconstructed", "y_reconstructed")),
    "zhuang1_f": dict(tag="zhuang_abca_1", label="Zhuang MERFISH 1, female", xy=("x", "y")),
    "zhuang2_m": dict(tag="zhuang_abca_2", label="Zhuang MERFISH 2, male", xy=("x", "y")),
}


def load_cells(ds):
    c = pd.read_parquet(D / f"{DATASETS[ds]['tag']}_cells.parquet")
    xk, yk = DATASETS[ds]["xy"]
    c = c.rename(columns={xk: "px", yk: "py"})
    c["neuron"] = ~c["class"].str.contains(NON_NEURONAL)
    c["vmh_type"] = c.subclass.isin(VMH_SUBCLASSES)
    c["in_ccf_vmh"] = c.parcellation_structure == "VMH"
    return c


def esr1_rich_clusters():
    e = pd.read_parquet(RESULTS / "tables" / "expr_by_type.parquet")
    e = e[(e.gene == "Esr1") & (e.chemistry == "10Xv3") & (e.sex == "all")]
    cl = e[e.level == "cluster"].set_index("type")
    st = e[e.level == "supertype"].set_index("type").excess_frac
    return cl, st


def polygon(mask, x0, y0):
    """Longest outline of a boolean grid mask as a simplified polygon in mm."""
    cs = measure.find_contours(np.pad(mask, 1).astype(float), 0.5)
    if not cs:
        return None
    c = max(cs, key=len) - 1
    c = measure.approximate_polygon(c, tolerance=0.6)
    if len(c) < 4:
        return None
    return [[round(float(x0 + p[0] * PIX), 4), round(float(y0 + p[1] * PIX), 4)] for p in c[:-1]]


def seeds_for_section(d, esr1_cluster, esr1_super):
    n = d[d.neuron]
    x0, y0 = d.px.min() - 0.1, d.py.min() - 0.1
    shape = (int((d.px.max() + 0.1 - x0) / PIX) + 1, int((d.py.max() + 0.1 - y0) / PIX) + 1)

    def dens(sub, w=None):
        g = np.zeros(shape)
        np.add.at(g, (((sub.px - x0) / PIX).astype(int), ((sub.py - y0) / PIX).astype(int)), 1.0 if w is None else w)
        return ndimage.gaussian_filter(g, SIGMA_MM / PIX)

    v = n[n.vmh_type]
    # Esr1 status: cluster level when the cluster has enough 10x cells, else supertype
    ex = v.cluster.map(esr1_cluster.excess_frac.where(esr1_cluster.n_cells >= 20))
    ex = ex.fillna(v.supertype.map(esr1_super)).fillna(0.0)
    rich = (ex >= ESR1_RICH).values.astype(float)
    r_all, r_vmh, r_vl = dens(n), dens(v), dens(v, rich)
    floor = 3.0 / (2 * np.pi * (SIGMA_MM / PIX) ** 2)               # about three VMH neurons within one sigma
    vmh = (r_vmh >= 0.5 * r_all) & (r_vmh >= floor)
    vmh = ndimage.binary_opening(vmh, iterations=1)
    # separate the hemispheres at the third ventricle, located by its tanycytes and ependymal cells
    wall = d[d.subclass.str.contains("Tanycyte|Ependymal", na=False)]
    mid_mm = wall.px.median() if len(wall) >= 10 else None
    if mid_mm is not None:
        mi = int((mid_mm - x0) / PIX)
        vmh[max(mi, 0):mi + 1, :] = False
    lab, k = ndimage.label(vmh)
    out = []
    comps = []
    for i in range(1, k + 1):
        m = lab == i
        if m.sum() * PIX ** 2 < 0.03:
            continue
        comps.append(ndimage.binary_fill_holes(m))
    if mid_mm is not None:
        mid = (mid_mm - x0) / PIX
    else:
        mid = np.mean([ndimage.center_of_mass(m)[0] for m in comps]) if len(comps) >= 2 else shape[0] / 2
    for m in comps:
        vl = m & (r_vl >= 0.5 * r_vmh)
        vl = ndimage.binary_opening(vl, iterations=1)
        if vl.any():
            l2, _ = ndimage.label(vl)
            sizes = ndimage.sum(vl, l2, range(1, l2.max() + 1))
            vl = l2 == (1 + int(np.argmax(sizes)))
            if vl.sum() * PIX ** 2 < 0.01:
                vl[:] = False
        rest = m & ~vl
        ii, jj = np.nonzero(rest)
        if len(ii) == 0:
            continue
        c_rest = np.array([ii.mean(), jj.mean()])
        if vl.any():
            axis = np.array(ndimage.center_of_mass(vl)) - c_rest
        else:                                      # no VMHvl seed here: assume it lies ventrolateral
            axis = np.array([np.sign(c_rest[0] - mid) or 1.0, 1.0])
        axis = axis / np.linalg.norm(axis)
        proj = (ii - c_rest[0]) * axis[0] + (jj - c_rest[1]) * axis[1]
        cut = 0.5 * (proj.min() + proj.max())
        dm = np.zeros_like(m); dm[ii[proj < cut], jj[proj < cut]] = True
        ce = np.zeros_like(m); ce[ii[proj >= cut], jj[proj >= cut]] = True
        for label, mask in (("VMHdm", dm), ("VMHc", ce), ("VMHvl", vl)):
            if mask.any():
                p = polygon(mask, x0, y0)
                if p:
                    out.append(dict(label=label, points=p))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cl, st = esr1_rich_clusters()
    seeds = {}
    for ds in DATASETS:
        c = load_cells(ds)
        seeds[ds] = {}
        for sec, d in c.groupby("brain_section_label"):
            polys = seeds_for_section(d, cl, st)
            seeds[ds][sec] = polys
        n = sum(len(v) for v in seeds[ds].values())
        print(f"{ds}: {len(seeds[ds])} sections, {n} seed polygons")
    with open(OUT / "seed_polygons.json", "w") as fh:
        json.dump(dict(version=1, units="mm, section coordinates of each dataset", datasets=seeds), fh)


if __name__ == "__main__":
    main()
