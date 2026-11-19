#!/usr/bin/env python3
"""Step 9: build the explorer.

docs/index.html   the page
docs/data.js      its data, loaded as a script
(docs/doc.js, the text of the Documentation tab, is built by step 11)

The files are static: they work opened from the folder, and from any web host that serves the
docs/ folder (GitHub Pages: branch main, folder /docs).

The page shows, for any curated gene: a map of neurons around the VMH in each spatial brain
(measured in situ where the gene is on that brain's panel, otherwise projected from the cell's
cluster), the aligned atlas compartments, expression by compartment and by supertype, the
gene's ligand or receptor partners, and the part 3 source tables.
"""
import base64
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import (COMPARTMENTS, COMPARTMENT_COLOURS, COMPARTMENT_NAMES, D, DATASETS, NON_NEURONAL, SCREENED, VMH_SUPERTYPES,  # noqa: E402
                    load_cells)
from paths import CONFIG, REPO, RESULTS  # noqa: E402

T, OUT = RESULTS / "tables", REPO / "docs"
COMPS = COMPARTMENTS + ["core"]
COMP_CODE = {"outside": 0, "VMHDM": 1, "VMHC": 2, "VMHVL": 3, "VMH": 4, "VMHSh": 5, "no_ccf": 6}
CROP_MM = 0.35
MIN_VMH_TYPE = 30


def enc(a):
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def r1(x):
    return None if x is None or not np.isfinite(x) else round(float(x) * 100, 1)


def gene_tables():
    lig = pd.read_csv(CONFIG / "neuropeptide_ligands.csv").fillna("")
    rec = pd.read_csv(CONFIG / "neuropeptide_receptors.csv").fillna("")
    ctx = pd.read_csv(CONFIG / "context_genes.csv").fillna("")
    pairs = pd.read_csv(CONFIG / "ligand_receptor_pairs.csv")
    genes = {}
    for r in lig.itertuples():
        genes[r.gene] = dict(kind="ligand", family=r.family, what=r.peptides, tag=r.source_class, note=r.note,
                             partners=list(dict.fromkeys(pairs[pairs.ligand_gene == r.gene].receptor_gene)))
        conditional = pairs[(pairs.ligand_gene == r.gene) & pairs.qualification.notna()]
        genes[r.gene]["pairNotes"] = "; ".join(f"{p.peptide}, {p.receptor_gene}: {p.qualification}"
                                               for p in conditional.itertuples())
    for r in rec.itertuples():
        genes.setdefault(r.gene, dict(kind="receptor", family=r.family, what=r.principal_coupling, tag="", note=r.note,
                                      partners=list(dict.fromkeys(pairs[pairs.receptor_gene == r.gene].ligand_gene))))
    for r in ctx.itertuples():
        genes.setdefault(r.gene, dict(kind="context", family="marker", what=r.role, tag="", note="", partners=[]))
    return genes


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = gene_tables()
    names = list(genes)

    # ---- per-gene numbers
    gb = pd.read_csv(T / "gene_by_compartment.csv")
    for ds in DATASETS:
        d = gb[gb.dataset == ds]
        P = d.pivot(index="gene", columns="compartment", values="proj_frac")
        M = d.pivot(index="gene", columns="compartment", values="meas_frac")
        C = d.pivot(index="gene", columns="compartment", values="call")
        for g in names:
            genes[g].setdefault("proj", {})[ds] = [r1(P.loc[g, c]) for c in COMPS]
            genes[g].setdefault("meas", {})[ds] = ([r1(M.loc[g, c]) for c in COMPS] if bool(d[d.gene == g].measured.iloc[0]) else None)
            gene_rows = d[(d.gene == g) & d.compartment.isin(COMPARTMENTS)]
            disagreements = gene_rows[gene_rows.evidence.str.startswith(("measured only", "projected only"))]
            genes[g].setdefault("discordant", {})[ds] = dict(zip(disagreements.compartment, disagreements.evidence))
            if (C.loc[g, COMPARTMENTS] == "ubiquitous").any():
                genes[g]["ubiquitous"] = True
    e = pd.read_parquet(T / "expr_by_type.parquet")
    for c in ("level", "type", "chemistry", "sex", "gene", "call"):
        e[c] = e[c].astype(str)
    s = e[(e.level == "supertype") & e.type.isin(VMH_SUPERTYPES) & (e.chemistry == "10Xv3") & (e.sex == "all")]
    X = s.pivot(index="gene", columns="type", values="excess_frac").reindex(columns=VMH_SUPERTYPES)
    K = s.pivot(index="gene", columns="type", values="call").reindex(columns=VMH_SUPERTYPES)
    code = {"broad": "B", "subset": "S", "sparse": "s", "glia_level": "g", "background": "", "too_few_cells": "?"}
    for g in names:
        genes[g]["st"] = [r1(v) for v in X.loc[g]]
        genes[g]["stc"] = [code.get(v, "?") for v in K.loc[g].fillna("too_few_cells")]
    src = pd.read_csv(T / "ligand_sources.csv").fillna("").set_index("gene")
    con = pd.read_csv(T / "receptor_ligand_concordance.csv").fillna("").set_index("receptor")
    for g in names:
        if g in src.index and genes[g]["kind"] == "ligand":
            r = src.loc[g]
            genes[g]["src"] = dict(local=bool(r.local), afferent=r.afferent_regions, endocrine=bool(r.endocrine),
                                   top=r.top_regions_brainwide, provisional=r.provisional_regions)
        if g in con.index and genes[g]["kind"] == "receptor":
            r = con.loc[g]
            genes[g]["src"] = dict(inVmh=bool(r.in_vmh), verdict=r.ligand_sources, local=r.local_ligands,
                                   afferent=r.afferent_ligands, endocrine=r.endocrine_ligands, conditional=r.conditional_partners)
        if "src" in genes[g]:
            genes[g]["src"].update(evidence=r.spatial_evidence, hit=bool(r.spatial_hit), peak=r1(r.peak_fraction),
                                   brain=r.peak_dataset, compartment=r.peak_compartment, method=r.peak_method)

    # ---- spatial cells
    lab = pd.read_parquet(D / "atlas_labels.parquet")
    shifts = pd.read_csv(T / "atlas_alignment.csv").set_index("dataset")
    outlines = json.load(open(D / "atlas_outlines.json"))["datasets"]
    proj = pd.read_parquet(D / "cluster_projection.parquet").set_index("cluster")
    st_code = {s: i + 1 for i, s in enumerate(VMH_SUPERTYPES)}
    clusters, brains = {}, {}
    for ds, meta in DATASETS.items():
        c = load_cells(ds).merge(lab[lab.dataset == ds][["cell_label", "kim_aligned"]], on="cell_label", how="left")
        z = np.load(D / f"{meta['tag']}_expr.npz", allow_pickle=True)
        Xp, panel = z["X"], list(z["genes"])
        assert (z["cell_label"] == c.cell_label.values).all()
        nn = c["class"].str.contains(NON_NEURONAL).values
        c["comp"] = c.kim_aligned.where(c.kim_aligned.isin(COMP_CODE), "outside")
        c["ap"] = c.x_ccf + shifts.loc[ds, "shift_ap_mm"]
        nv = c[c.neuron & c.vmh_type].groupby("brain_section_label").size()
        secs = [s for s in nv.index if nv[s] >= MIN_VMH_TYPE]
        inv = c[c.comp.isin(COMPARTMENTS)]
        x0, x1, y0, y1 = inv.px.min() - CROP_MM, inv.px.max() + CROP_MM, inv.py.min() - CROP_MM, inv.py.max() + CROP_MM
        keep = (c.neuron & c.brain_section_label.isin(secs) & c.px.between(x0, x1) & c.py.between(y0, y1)).values
        k = c[keep].copy()
        k["row"] = np.flatnonzero(keep)
        order = k.groupby("brain_section_label").ap.median().sort_values().index      # anterior to posterior
        k["sec_rank"] = k.brain_section_label.map({s: i for i, s in enumerate(order)})
        k = k.sort_values(["sec_rank", "row"])
        for cl in k.cluster.unique():
            clusters.setdefault(cl, len(clusters))
        measured = {}
        for g in names:
            if g in panel:
                col = Xp[:, panel.index(g)].astype(float)
                thr = float(np.percentile(col[nn], 95))
                v = col[k.row.values]
                top = float(np.percentile(v[v > 0], 99.5)) if (v > 0).any() else 1.0
                top = max(top, thr * 1.05, 1e-6)
                measured[g] = dict(v=enc(np.round(255 * np.clip(v / top, 0, 1)).astype(np.uint8)), max=top, thr=thr,
                                   above=enc(np.packbits(v > thr, bitorder="little")))
        sections, off = [], 0
        for sname in order:
            d = k[k.brain_section_label == sname]
            sections.append(dict(id=sname, short=sname.rsplit(".", 1)[-1], ap=round(float(d.ap.median()), 2), n=len(d), off=off,
                                 nv=int(d.vmh_type.sum()), outlines=outlines[ds].get(sname, {})))
            off += len(d)
        brains[ds] = dict(
            label=meta["label"], sex=str(c.donor_sex.iloc[0]), n=len(k), sections=sections,
            x=enc(np.round(k.px.values * 1000).astype(np.int16)), y=enc(np.round(k.py.values * 1000).astype(np.int16)),
            cl=enc(k.cluster.map(clusters).values.astype(np.uint16)),
            st=enc(k.supertype.map(st_code).fillna(0).values.astype(np.uint8)),
            cp=enc(k.comp.map(COMP_CODE).values.astype(np.uint8)), measured=measured,
            shift=[float(shifts.loc[ds, f"shift_{a}_mm"]) for a in ("ap", "dv", "ml")])
        print(f"{ds}: {len(k)} neurons in {len(sections)} sections; {len(measured)} measured genes")
    cl_names = sorted(clusters, key=clusters.get)
    P = proj.reindex(cl_names)
    projected = {g: enc(np.round(100 * P[f"frac:{g}"].fillna(0).values).astype(np.uint8)) for g in names}
    # where each cluster's projected value comes from: its own sequenced cells, its supertype's or its
    # subclass's (too few cells of its own), or nowhere (not in the hypothalamus sequencing data)
    proj_src = enc(P.source.str.split(":").str[0].map({"cluster": 0, "supertype": 1, "subclass": 2}).fillna(3).values.astype(np.uint8))

    # what blank probes reach in each compartment (percent of neurons), for the brains whose panel has them
    fl = pd.read_csv(T / "blank_probe_floor.csv")
    blank_floor = {ds: {r.compartment: round(100 * r.blank_frac_p95, 1) for r in d.itertuples()} for ds, d in fl.groupby("dataset")}
    data = dict(
        screened=SCREENED, blankFloor=blank_floor,
        built="Allen Brain Cell Atlas release 20250531; unified atlas of Chon et al. 2019",
        genes=genes, order=names, comps=COMPS, compNames=[COMPARTMENT_NAMES.get(c, "core (three compartments and ends)") for c in COMPS],
        compColours=COMPARTMENT_COLOURS, compCode=COMP_CODE,
        supertypes=[s[9:].replace(" Glut", "") for s in VMH_SUPERTYPES], clusters=cl_names, projected=projected, projSrc=proj_src, brains=brains)
    blob = json.dumps(data, separators=(",", ":"))
    tpl = (HERE / "explorer_template.html").read_text()
    (OUT / "data.js").write_text("window.__DATA__ = " + blob + ";\n")
    (OUT / "index.html").write_text(tpl)
    (OUT / ".nojekyll").write_text("")          # GitHub Pages: serve the files as they are
    print(f"wrote {OUT}: index.html and data.js {len(blob) / 1e6:.1f} MB, {len(names)} genes, {len(cl_names)} clusters")


if __name__ == "__main__":
    main()
