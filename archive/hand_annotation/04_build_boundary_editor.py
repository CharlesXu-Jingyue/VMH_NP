#!/usr/bin/env python3
"""Step 4: build the self-contained boundary editor (results/boundaries/boundary_editor.html).

The page shows every spatial section around the VMH, coloured by cell type or by a marker
gene, with the seed regions from step 3 as editable polygons. The exported
vmh_boundaries.json goes to config/ and defines "VMH" and its subdivisions downstream.
"""
import base64
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from skimage import measure

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import RESULTS, data_root  # noqa: E402

spec = importlib.util.spec_from_file_location("seeds", HERE / "03_boundary_seeds.py")
seeds_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(seeds_mod)

D = data_root() / "derived"
OUT = RESULTS / "boundaries"
VMH_SUPERTYPES = {   # display order and colours
    "0568 VMH Nr5a1 Glut_3": "#3b6fb6", "0569 VMH Nr5a1 Glut_4": "#58c4e0", "0570 VMH Nr5a1 Glut_5": "#1e8a3c",
    "0567 VMH Nr5a1 Glut_2": "#c9b22b", "0571 VMH Nr5a1 Glut_6": "#a8328a", "0566 VMH Nr5a1 Glut_1": "#7a7a7a",
    "0564 VMH Fezf1 Glut_2": "#f08030", "0565 VMH Fezf1 Glut_3": "#d1301c", "0563 VMH Fezf1 Glut_1": "#111111",
}
NEIGHBOUR_COLOURS = ["#b9d3ee", "#c8e6c9", "#e6ccf0", "#f6d6b8", "#d9d2a6", "#bfe3e0", "#f4c2cf", "#cfcfe8"]
GENES = {   # marker layers offered in the editor; "measured" where the panel has the gene, else cluster-projected
    "allen_m": ["Esr1", "Nr5a1", "Fezf1", "Lepr", "Satb2", "Tac1", "Adcyap1", "Oxtr", "Dlk1", "Sox14",
                "Cartpt", "Gal", "Pdyn", "Penk", "Nts", "Slc17a6", "Slc32a1"],
    "zhuang1_f": ["Esr1", "Nr5a1", "Fezf1", "Lepr", "Satb2", "Ar", "Tac1", "Adcyap1", "Oxtr", "Cckbr", "Crhr2", "Npy1r",
                  "Dlk1", "Cartpt", "Slc17a6", "Slc32a1"],
    "zhuang2_m": ["Esr1", "Nr5a1", "Fezf1", "Lepr", "Satb2", "Ar", "Tac1", "Adcyap1", "Oxtr", "Cckbr", "Crhr2", "Npy1r",
                  "Dlk1", "Cartpt", "Slc17a6", "Slc32a1"],
}
ZKEY = {"allen_m": "z_reconstructed", "zhuang1_f": "z", "zhuang2_m": "z"}


def enc(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def ccf_outlines():
    """VMH outline of the CCF annotation on each Allen MERFISH section, in mm."""
    A = np.load(D / "anatomy.npz", allow_pickle=True)
    crop, org, vmh = A["merfish_crop"], A["merfish_crop_origin"], int(A["vmh_index"])
    out = {}
    for k in range(crop.shape[2]):
        lines = []
        for c in measure.find_contours(np.pad(crop[:, :, k] == vmh, 1).astype(float), 0.5):
            c = measure.approximate_polygon(c - 1, tolerance=0.5)
            # voxel centres sit at integer multiples of 10 um
            lines.append([[round(float((p[0] + org[0]) * 0.01), 4), round(float((p[1] + org[1]) * 0.01), 4)] for p in c])
        out[round(float((k + org[2]) * 0.2), 1)] = lines
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = json.load(open(OUT / "seed_polygons.json"))["datasets"]
    expr = pd.read_parquet(RESULTS / "tables" / "expr_by_type.parquet")
    expr = expr[(expr.chemistry == "10Xv3") & (expr.sex == "all") & (expr.n_cells >= 20)]
    proj_cluster = expr[expr.level == "cluster"].pivot(index="type", columns="gene", values="mean_log2cpm")
    proj_super = expr[expr.level == "supertype"].pivot(index="type", columns="gene", values="mean_log2cpm")
    cells = {ds: seeds_mod.load_cells(ds) for ds in seeds_mod.DATASETS}
    # neighbouring subclasses worth their own colour: the most frequent non-VMH neuronal subclasses overall
    allc = pd.concat([c[c.neuron & ~c.vmh_type] for c in cells.values()])
    neigh = list(allc.subclass.value_counts().index[:len(NEIGHBOUR_COLOURS)])
    types = ([dict(name="non-neuronal", color="#ececec", layer=0), dict(name="other neurons", color="#cfcfcf", layer=1)]
             + [dict(name=s[4:], color=col, layer=2) for s, col in zip(neigh, NEIGHBOUR_COLOURS)]
             + [dict(name=s[5:], color=col, layer=3) for s, col in VMH_SUPERTYPES.items()])
    code = {**{s: 2 + i for i, s in enumerate(neigh)}, **{s: 2 + len(neigh) + i for i, s in enumerate(VMH_SUPERTYPES)}}
    ccf = ccf_outlines()
    data = dict(source="Allen Brain Cell Atlas release 20250531; seeds from scripts/03_boundary_seeds.py",
                types=types, datasets={})
    for ds, meta in seeds_mod.DATASETS.items():
        c = cells[ds]
        z = np.load(D / f"{meta['tag']}_expr.npz", allow_pickle=True)
        X, panel = z["X"], list(z["genes"])
        assert (z["cell_label"] == c.cell_label.values).all()
        tcode = np.where(c.neuron, 1, 0).astype(np.uint8)
        for key, col in (("subclass", neigh), ("supertype", VMH_SUPERTYPES)):
            m = c[key].map(code)
            tcode = np.where(m.notna() & c.neuron, m.fillna(0).astype(int), tcode).astype(np.uint8)
        layers, ginfo = {}, []
        for g in GENES[ds]:
            if g in panel:
                v, kind, unit = X[:, panel.index(g)].astype(float), "measured", "log2 (in situ)"
            else:
                v = c.cluster.map(proj_cluster[g]) if g in proj_cluster else pd.Series(np.nan, index=c.index)
                v = v.fillna(c.supertype.map(proj_super[g])).fillna(0.0).values.astype(float)
                kind, unit = "projected", "log2(CPM+1), cluster mean"
            top = float(np.percentile(v[v > 0], 99.5)) if (v > 0).any() else 1.0
            layers[g] = np.round(255 * np.clip(v / top, 0, 1)).astype(np.uint8)
            ginfo.append(dict(name=g, kind=kind, max=top, unit=unit))
        zsec = c.groupby("brain_section_label")[ZKEY[ds]].median().sort_values()
        sections = []
        for sec, zval in zsec.items():
            m = (c.brain_section_label == sec).values
            d = c[m]
            sections.append(dict(
                id=sec, pos=f"z {zval:.2f} mm", n=int(m.sum()), nv=int((d.neuron & d.vmh_type).sum()),
                x=enc(np.round(d.px.values * 1000).astype(np.int16)), y=enc(np.round(d.py.values * 1000).astype(np.int16)),
                t=enc(tcode[m]), f=enc(d.in_ccf_vmh.values.astype(np.uint8)),
                g={g: enc(layers[g][m]) for g in layers},
                ccf=ccf.get(round(float(zval), 1), []) if ds == "allen_m" else [],
                seeds=seeds[ds].get(sec, [])))
        data["datasets"][ds] = dict(label=meta["label"], genes=ginfo, sections=sections)
        print(f"{ds}: {len(sections)} sections, {len(c)} cells, {len(ginfo)} gene layers "
              f"({sum(g['kind'] == 'measured' for g in ginfo)} measured)")
    html = (HERE / "boundary_editor_template.html").read_text().replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")))
    out = OUT / "boundary_editor.html"
    out.write_text(html)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
