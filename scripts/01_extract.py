#!/usr/bin/env python3
"""Step 1: reduce the Allen Brain Cell Atlas downloads to compact VMH-centred extracts.

Reads <data_root>/abc_atlas, writes <data_root>/derived. Each stage is skipped when its
outputs already exist (delete the files to recompute), and stages are small enough to
run one at a time in a memory-limited (about 3 GB) or time-limited environment.

Stage groups
  merfish_*   Allen MERFISH (male, 500 target genes and 50 blank controls): cells around VMH with all panel genes,
              plus whole-brain counts of cells per (cluster, CCF structure)
  zhuang*_*   Zhuang MERFISH brains 1 (female) and 2 (male), 1,122 genes: same
  tenx_*      10x v2/v3 hypothalamus matrices: curated genes for every QC-passed cell
  clusters    whole-brain per-cluster mean expression of curated genes; taxonomy table
  anatomy     CCF annotation crops and depth-inside-VMH for every extracted cell
"""
import argparse
import glob
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import CONFIG, abc_file, abc_root, data_root  # noqa: E402

OUT = data_root() / "derived"
VMH = "VMH"
PAD_MM = 0.6          # in-plane margin kept around the VMH bounding box
VMH_SUBCLASSES = ("128 VMH Fezf1 Glut", "129 VMH Nr5a1 Glut")
MIN_VMH_TYPE = 20     # VMH-subclass neurons needed for a section to count as containing VMH
FORCE = False
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def curated_genes():
    genes = []
    for name in ("neuropeptide_ligands", "neuropeptide_receptors", "context_genes"):
        genes += list(pd.read_csv(CONFIG / f"{name}.csv").gene)
    return list(dict.fromkeys(genes))


def vmh_parcellation_index():
    f = abc_file("Allen-CCF-2020", "parcellation_to_parcellation_term_membership")
    m = pd.read_csv(f)
    idx = m.loc[(m.parcellation_term_acronym == VMH) & (m.parcellation_term_set_name == "structure"),
                "parcellation_index"].unique()
    assert len(idx) == 1, idx
    return int(idx[0]), m


def vmh_window(df, xcol, ycol, zcol):
    """Sections and in-plane window that contain the molecularly defined VMH.

    The CCF-registered VMH label turned out to be offset from the territory of the VMH cell
    types by one to three sections, so sections are chosen by where neurons of the two VMH
    subclasses actually are: every section with at least MIN_VMH_TYPE such cells, plus one
    flanking section on each side. The window is their bounding box (1st-99th percentile)
    padded by PAD_MM, also covering CCF-labelled VMH cells in those sections.
    """
    is_type = df.subclass.isin(VMH_SUBCLASSES)
    per_sec = df[is_type].groupby("brain_section_label").size()
    core = per_sec[per_sec >= MIN_VMH_TYPE].index
    zsec = df.groupby("brain_section_label")[zcol].median().sort_values()
    order = list(zsec.index)
    pos = sorted(order.index(s) for s in core)
    keep = order[max(pos[0] - 1, 0):pos[-1] + 2]
    d = df[df.brain_section_label.isin(keep)]
    ref = d[d.subclass.isin(VMH_SUBCLASSES) | (d.parcellation_structure == VMH)]
    x0, x1 = ref[xcol].quantile(0.01) - PAD_MM, ref[xcol].quantile(0.99) + PAD_MM
    y0, y1 = ref[ycol].quantile(0.01) - PAD_MM, ref[ycol].quantile(0.99) + PAD_MM
    out = d[d[xcol].between(x0, x1) & d[ycol].between(y0, y1)].reset_index(drop=True)
    return out, keep, list(core), (x0, x1, y0, y1)


def read_csr_rows(h5, rows):
    """Rows (sorted, unique) of an AnnData CSR matrix stored in an open h5py file."""
    X = h5["X"]
    n_var = int(X.attrs["shape"][1])
    indptr = X["indptr"][:].astype(np.int64)
    blocks = []
    # read in contiguous spans of at most ~40M stored values
    start = 0
    while start < len(rows):
        end = start + 1
        while end < len(rows) and indptr[rows[end] + 1] - indptr[rows[start]] < 40_000_000:
            end += 1
        r0, r1 = rows[start], rows[end - 1] + 1
        lo, hi = indptr[r0], indptr[r1]
        data = X["data"][lo:hi]
        indices = X["indices"][lo:hi]
        block = sp.csr_matrix((data, indices, indptr[r0:r1 + 1] - lo), shape=(r1 - r0, n_var))
        blocks.append(block[rows[start:end] - r0])
        start = end
    return sp.vstack(blocks).tocsr()


# ----------------------------------------------------------------------------- MERFISH
def stage_merfish_meta():
    near_out, xtab_out = OUT / "_merfish_window.parquet", OUT / "merfish_cluster_by_structure.parquet"
    if near_out.exists() and xtab_out.exists() and not FORCE:
        return log("merfish_meta: outputs exist, skipping")
    meta = abc_file("MERFISH-C57BL6J-638850-CCF", "cell_metadata_with_parcellation_annotation")
    use = ["cell_label", "brain_section_label", "cluster_alias", "average_correlation_score", "donor_sex",
           "neurotransmitter", "class", "subclass", "supertype", "cluster",
           "x_reconstructed", "y_reconstructed", "z_reconstructed", "x_ccf", "y_ccf", "z_ccf",
           "parcellation_index", "parcellation_division", "parcellation_structure"]
    keep, xtab = [], []
    for ch in pd.read_csv(meta, usecols=use, dtype={"cell_label": str}, chunksize=400_000):
        xtab.append(ch.groupby(["cluster_alias", "parcellation_division", "parcellation_structure"], dropna=False).size())
        keep.append(ch[ch.z_reconstructed.between(5.0, 8.6)])          # generous slab around the hypothalamus
    xt = pd.concat(xtab).groupby(level=[0, 1, 2], dropna=False).sum().rename("n_cells").reset_index()
    xt.to_parquet(xtab_out)
    near = pd.concat(keep, ignore_index=True)
    cells, secs, core, (x0, x1, y0, y1) = vmh_window(near, "x_reconstructed", "y_reconstructed", "z_reconstructed")
    cells.to_parquet(near_out)
    log(f"merfish_meta: sections kept {secs}; with VMH cell types {core}; window x[{x0:.2f},{x1:.2f}] "
        f"y[{y0:.2f},{y1:.2f}]; {len(cells)} cells kept, {int((cells.parcellation_structure == VMH).sum())} "
        f"CCF-labelled VMH, {int(cells.subclass.isin(VMH_SUBCLASSES).sum())} of VMH subclasses")


def stage_merfish_expr():
    cells_out, expr_out = OUT / "merfish_cells.parquet", OUT / "merfish_expr.npz"
    if cells_out.exists() and expr_out.exists() and not FORCE:
        return log("merfish_expr: outputs exist, skipping")
    cells = pd.read_parquet(OUT / "_merfish_window.parquet")
    with h5py.File(abc_file("MERFISH-C57BL6J-638850", "C57BL6J-638850/log2"), "r") as h5:
        labels = pd.Index(h5["obs/cell_label"][:].astype(str))
        genes = h5["var/gene_symbol"][:].astype(str)
        rows = labels.get_indexer(cells.cell_label)
        assert (rows >= 0).all(), "MERFISH cells missing from the expression matrix"
        order = np.argsort(rows)
        X = read_csr_rows(h5, rows[order])
    cells = cells.iloc[order].reset_index(drop=True)
    cells.to_parquet(cells_out)
    np.savez_compressed(expr_out, X=X.toarray().astype(np.float32), genes=genes, cell_label=cells.cell_label.values)
    log(f"merfish_expr: wrote {X.shape[0]} cells x {X.shape[1]} genes")


# ------------------------------------------------------------------------------ Zhuang
def zhuang_meta(z):
    tag = z.lower().replace("-", "_")
    near_out, xtab_out = OUT / f"_{tag}_window.parquet", OUT / f"{tag}_cluster_by_structure.parquet"
    if near_out.exists() and xtab_out.exists() and not FORCE:
        return log(f"{z} meta: outputs exist, skipping")
    vmh_idx, memb = vmh_parcellation_index()
    struct = (memb[memb.parcellation_term_set_name == "structure"]
              .drop_duplicates("parcellation_index").set_index("parcellation_index").parcellation_term_acronym)
    division = (memb[memb.parcellation_term_set_name == "division"]
                .drop_duplicates("parcellation_index").set_index("parcellation_index").parcellation_term_acronym)
    meta = pd.read_csv(abc_file(z, "cell_metadata_with_cluster_annotation"), dtype={"cell_label": str},
                       usecols=["cell_label", "brain_section_label", "donor_sex", "cluster_alias", "x", "y", "z",
                                "subclass_confidence_score", "cluster_confidence_score", "high_quality_transfer",
                                "neurotransmitter", "class", "subclass", "supertype", "cluster"])
    ccf = pd.read_csv(abc_file(z + "-CCF", "ccf_coordinates"), dtype={"cell_label": str})
    ccf = ccf.rename(columns={"x": "x_ccf", "y": "y_ccf", "z": "z_ccf"})
    meta = meta.merge(ccf, on="cell_label", how="left")
    meta["parcellation_structure"] = meta.parcellation_index.map(struct)
    meta["parcellation_division"] = meta.parcellation_index.map(division)
    (meta.groupby(["cluster_alias", "parcellation_division", "parcellation_structure"], dropna=False).size()
         .rename("n_cells").reset_index().to_parquet(xtab_out))
    cells, secs, core, (x0, x1, y0, y1) = vmh_window(meta, "x", "y", "z")
    cells.to_parquet(near_out)
    log(f"{z} meta: sex {sorted(meta.donor_sex.unique())}; {len(meta)} cells, {int(meta.parcellation_index.notna().sum())} "
        f"with CCF coordinates; sections kept {secs}; with VMH cell types {core}; {len(cells)} cells kept, "
        f"{int((cells.parcellation_index == vmh_idx).sum())} CCF-labelled VMH, "
        f"{int(cells.subclass.isin(VMH_SUBCLASSES).sum())} of VMH subclasses")


def zhuang_expr(z):
    tag = z.lower().replace("-", "_")
    cells_out, expr_out = OUT / f"{tag}_cells.parquet", OUT / f"{tag}_expr.npz"
    if cells_out.exists() and expr_out.exists() and not FORCE:
        return log(f"{z} expr: outputs exist, skipping")
    cells = pd.read_parquet(OUT / f"_{tag}_window.parquet")
    with h5py.File(abc_file(z, f"{z}/log2"), "r") as h5:
        labels = pd.Index(h5["obs/cell_label"][:].astype(str))
        genes = h5["var/gene_symbol"][:].astype(str)
        rows = labels.get_indexer(cells.cell_label)
        assert (rows >= 0).all()
        order = np.argsort(rows)
        rows = rows[order]
        X = np.empty((len(rows), len(genes)), dtype=np.float32)
        start = 0                       # dense, chunk-compressed: read contiguous row spans, then pick rows
        while start < len(rows):
            end = start + 1
            while end < len(rows) and rows[end] - rows[start] < 60_000:
                end += 1
            r0, r1 = rows[start], rows[end - 1] + 1
            X[start:end] = h5["X"][r0:r1][rows[start:end] - r0]
            start = end
    cells = cells.iloc[order].reset_index(drop=True)
    cells.to_parquet(cells_out)
    np.savez_compressed(expr_out, X=X, genes=genes, cell_label=cells.cell_label.values)
    log(f"{z} expr: wrote {X.shape[0]} cells x {X.shape[1]} genes")


# --------------------------------------------------------------------------------- 10x
TENX = {"WMB-10Xv3-HY": ("WMB-10Xv3", "WMB-10Xv3-HY/log2"), "WMB-10Xv2-HY": ("WMB-10Xv2", "WMB-10Xv2-HY/log2")}


def stage_tenx_meta():
    meta_out, xtab_out = OUT / "_tenx_meta.parquet", OUT / "tenx_cluster_by_matrix.parquet"
    if meta_out.exists() and xtab_out.exists() and not FORCE:
        return log("tenx_meta: outputs exist, skipping")
    use = ["cell_label", "library_label", "feature_matrix_label", "library_method", "region_of_interest_acronym",
           "donor_label", "donor_sex", "cluster_alias", "neurotransmitter", "class", "subclass", "supertype", "cluster"]
    keep, xtab = [], []
    for ch in pd.read_csv(abc_file("WMB-10X", "cell_metadata_with_cluster_annotation"), usecols=use,
                          dtype={"cell_label": str}, chunksize=400_000):
        xtab.append(ch.groupby(["cluster_alias", "feature_matrix_label"]).size())
        keep.append(ch[ch.feature_matrix_label.isin(TENX)])
    pd.concat(xtab).groupby(level=[0, 1]).sum().rename("n_cells").reset_index().to_parquet(xtab_out)
    meta = pd.concat(keep, ignore_index=True)
    meta.to_parquet(meta_out)
    log(f"tenx_meta: {len(meta)} QC-passed cells in the hypothalamus matrices")


def tenx_matrix(name):
    part = OUT / f"_tenx_{name}.npz"
    if part.exists() and not FORCE:
        return log(f"tenx {name}: output exists, skipping")
    directory, fname = TENX[name]
    want = curated_genes()
    blocks = []
    with h5py.File(abc_file(directory, fname), "r") as h5:
        labels = h5["obs/cell_label"][:].astype(str)
        symbols = h5["var/gene_symbol"][:].astype(str)
        cols = np.array([np.flatnonzero(symbols == g)[0] for g in want])
        X = h5["X"]
        indptr = X["indptr"][:].astype(np.int64)
        n = len(labels)
        lut = np.full(len(symbols), -1, dtype=np.int32)
        lut[cols] = np.arange(len(cols))
        step = 4000
        for r0 in range(0, n, step):
            r1 = min(r0 + step, n)
            lo, hi = indptr[r0], indptr[r1]
            newcol = lut[X["indices"][lo:hi]]
            sel = newcol >= 0
            data = X["data"][lo:hi][sel]
            rowid = np.repeat(np.arange(r1 - r0), np.diff(indptr[r0:r1 + 1]))[sel]   # row of every kept value
            blocks.append(sp.csr_matrix((data, (rowid, newcol[sel])), shape=(r1 - r0, len(cols)), dtype=np.float32))
    M = sp.vstack(blocks).tocsr()
    np.savez_compressed(part, data=M.data, indices=M.indices, indptr=M.indptr, shape=M.shape,
                        cell_label=labels, genes=np.array(want))
    log(f"tenx {name}: read {n} cells x {len(want)} genes")


def stage_tenx_merge():
    cells_out, expr_out = OUT / "tenx_cells.parquet", OUT / "tenx_expr.npz"
    if cells_out.exists() and expr_out.exists() and not FORCE:
        return log("tenx_merge: outputs exist, skipping")
    meta = pd.read_parquet(OUT / "_tenx_meta.parquet").set_index("cell_label")
    obs, mats, genes = [], [], None
    for name in TENX:
        z = np.load(OUT / f"_tenx_{name}.npz", allow_pickle=False)
        mats.append(sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])))
        obs.append(pd.DataFrame({"cell_label": z["cell_label"], "matrix": name}))
        genes = z["genes"]
    obs = pd.concat(obs, ignore_index=True)
    X = sp.vstack(mats).tocsr()
    ok = obs.cell_label.isin(meta.index).values
    n_all = len(obs)
    obs, X = obs[ok].reset_index(drop=True), X[np.flatnonzero(ok)]
    cells = obs.join(meta, on="cell_label")
    assert (cells.matrix == cells.feature_matrix_label).all()
    cells.drop(columns="matrix").to_parquet(cells_out)
    sp.save_npz(expr_out, X)
    pd.Series(genes).to_csv(OUT / "tenx_genes.csv", index=False, header=["gene"])
    log(f"tenx_merge: {X.shape[0]} of {n_all} cells in the matrices pass QC; {X.shape[1]} genes; {X.nnz} stored values")


# ---------------------------------------------------------------------------- clusters
def stage_clusters():
    means_out, tax_out = OUT / "cluster_means.parquet", OUT / "cluster_taxonomy.parquet"
    if means_out.exists() and tax_out.exists() and not FORCE:
        return log("clusters: outputs exist, skipping")
    piv = pd.read_csv(abc_file("WMB-taxonomy", "cluster_to_cluster_annotation_membership_pivoted"))
    clu = pd.read_csv(abc_file("WMB-taxonomy", "cluster"))
    tax = clu.merge(piv, on="cluster_alias")
    tax.to_parquet(tax_out)
    gene = pd.read_csv(abc_file("WMB-10X", "gene"))
    want = curated_genes()
    sub = gene[gene.gene_symbol.isin(want)]
    assert sub.gene_symbol.is_unique and len(sub) == len(want), "curated symbols must map to exactly one gene"
    ens = sub.set_index("gene_symbol").gene_identifier.reindex(want)
    with h5py.File(abc_file("WMB-10X", "precomputed_stats_ABC_revision_230821"), "r") as h5:
        c2r = json.loads(h5["cluster_to_row"][()])
        cols = json.loads(h5["col_names"][()])
        pos = {g: i for i, g in enumerate(cols)}
        ci = np.array([pos[e] for e in ens])
        n = h5["sum"].shape[0]
        # n_cells is 1 per row: the "sum" rows are cluster means of log2(CPM+1)
        M = np.vstack([h5["sum"][r:r + 400][:, ci] for r in range(0, n, 400)])
        assert (h5["n_cells"][:] == 1).all()
    df = pd.DataFrame(M.astype(np.float32), columns=want)
    df.insert(0, "cluster_label", pd.Series(c2r).sort_values().index)
    df.to_parquet(means_out)
    log(f"clusters: wrote means for {df.shape[0]} clusters x {len(want)} genes; taxonomy {tax.shape}")


# ----------------------------------------------------------------------------- anatomy
def stage_anatomy():
    import nibabel as nib
    from scipy import ndimage
    out = OUT / "anatomy.npz"
    depth_out = OUT / "cell_depth.parquet"
    if out.exists() and depth_out.exists() and not FORCE:
        return log("anatomy: outputs exist, skipping")
    vmh_idx, memb = vmh_parcellation_index()
    base = abc_root()
    # (a) annotation resampled into the Allen MERFISH section space (10 um in plane, one slice per section)
    f = glob.glob(str(base / "image_volumes/MERFISH-C57BL6J-638850-CCF/*/resampled_annotation.nii.gz"))[0]
    img = nib.load(f)
    a = np.asanyarray(img.dataobj).astype(np.int32)
    zz = np.unique(np.where(a == vmh_idx)[2])
    z0, z1 = zz.min() - 1, zz.max() + 2
    m = np.where((a[:, :, zz.min():zz.max() + 1] == vmh_idx).any(axis=2))
    pad = int(PAD_MM / 0.01)
    i0, i1 = max(m[0].min() - pad, 0), min(m[0].max() + pad, a.shape[0])
    j0, j1 = max(m[1].min() - pad, 0), min(m[1].max() + pad, a.shape[1])
    crop = a[i0:i1, j0:j1, z0:z1].astype(np.int16)
    # (b) depth inside VMH in CCF space (10 um isotropic): positive inside, negative outside.
    # The 10 um volume is too large to hold in memory, so locate VMH at 25 um first.
    f25 = glob.glob(str(base / "image_volumes/Allen-CCF-2020/*/annotation_25.nii.gz"))[0]
    A25 = np.asanyarray(nib.load(f25).dataobj)
    w = np.where(A25 == vmh_idx)
    p = 80                                                    # 0.8 mm margin, in 10 um voxels
    f10 = glob.glob(str(base / "image_volumes/Allen-CCF-2020/*/annotation_10.nii.gz"))[0]
    img10 = nib.load(f10)
    shape10 = img10.shape
    lo = [max(int(w[k].min() * 2.5) - p, 0) for k in range(3)]
    hi = [min(int((w[k].max() + 1) * 2.5) + p, shape10[k]) for k in range(3)]
    sub = np.asarray(img10.dataobj[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]).astype(np.int32)
    del A25
    mask = sub == vmh_idx
    assert mask.any()
    depth = (ndimage.distance_transform_edt(mask) - ndimage.distance_transform_edt(~mask)) * 0.01   # mm
    names = memb[memb.parcellation_term_set_name == "structure"].drop_duplicates("parcellation_index")
    np.savez_compressed(out, merfish_crop=crop, merfish_crop_origin=np.array([i0, j0, z0]),
                        merfish_affine=img.affine, ccf_sub=sub.astype(np.int16), ccf_sub_origin=np.array(lo),
                        ccf_affine=img10.affine, ccf_depth=depth.astype(np.float32), vmh_index=vmh_idx,
                        parcellation_index=names.parcellation_index.values,
                        parcellation_acronym=names.parcellation_term_acronym.values.astype(str))
    log(f"anatomy: merfish crop {crop.shape}, ccf sub-volume {sub.shape}, VMH volume "
        f"{mask.sum() * 1e-6:.3f} mm^3; 10 um affine:\n{img10.affine}")
    # sample depth and the annotation label at every extracted cell
    inv = np.linalg.inv(img10.affine)
    rows = []
    for tag in ("merfish", "zhuang_abca_1", "zhuang_abca_2"):
        fcell = OUT / f"{tag}_cells.parquet"
        if not fcell.exists():
            continue
        c = pd.read_parquet(fcell, columns=["cell_label", "x_ccf", "y_ccf", "z_ccf", "parcellation_index"])
        xyz = c[["x_ccf", "y_ccf", "z_ccf"]].values.astype(float)
        best = None
        # voxel index from mm: plain axis order, or through the NIfTI affine; keep whichever
        # reproduces the released per-cell parcellation labels
        for mode in ("plain", "affine"):
            if mode == "plain":
                ijk = np.floor(xyz / 0.01)
            else:
                ijk = np.round((np.c_[xyz, np.ones(len(xyz))] @ inv.T)[:, :3])
            rel = ijk - np.array(lo)
            ins = np.isfinite(rel).all(1)
            rel_i = np.where(ins[:, None], rel, 0).astype(int)
            ins &= (rel_i >= 0).all(1) & (rel_i < np.array(sub.shape)).all(1)
            lab = sub[rel_i[ins, 0], rel_i[ins, 1], rel_i[ins, 2]]
            agree = float((lab == c.parcellation_index.values[ins]).mean()) if ins.any() else 0.0
            if best is None or agree > best[0]:
                best = (agree, mode, rel_i, ins)
        agree, mode, rel_i, ins = best
        d = np.full(len(c), np.nan, dtype=np.float32)
        d[ins] = depth[rel_i[ins, 0], rel_i[ins, 1], rel_i[ins, 2]]
        rows.append(pd.DataFrame({"dataset": tag, "cell_label": c.cell_label, "vmh_depth_mm": d}))
        log(f"anatomy: {tag}: coordinate convention '{mode}', label agreement {agree:.4f}, "
            f"{int(ins.sum())} of {len(c)} cells inside the sub-volume")
    pd.concat(rows, ignore_index=True).to_parquet(depth_out)


STAGES = {
    "merfish_meta": stage_merfish_meta, "merfish_expr": stage_merfish_expr,
    "zhuang1_meta": lambda: zhuang_meta("Zhuang-ABCA-1"), "zhuang1_expr": lambda: zhuang_expr("Zhuang-ABCA-1"),
    "zhuang2_meta": lambda: zhuang_meta("Zhuang-ABCA-2"), "zhuang2_expr": lambda: zhuang_expr("Zhuang-ABCA-2"),
    "tenx_meta": stage_tenx_meta, "tenx_v3": lambda: tenx_matrix("WMB-10Xv3-HY"),
    "tenx_v2": lambda: tenx_matrix("WMB-10Xv2-HY"), "tenx_merge": stage_tenx_merge,
    "clusters": stage_clusters, "anatomy": stage_anatomy,
}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stages", nargs="*", help=f"any of: {', '.join(STAGES)} (default: all, in order)")
    ap.add_argument("--force", action="store_true", help="recompute and overwrite existing outputs")
    args = ap.parse_args()
    FORCE = args.force
    unknown = [x for x in args.stages if x not in STAGES]
    if unknown:
        ap.error(f"unknown stage(s): {unknown}")
    OUT.mkdir(parents=True, exist_ok=True)
    for name in args.stages or list(STAGES):
        STAGES[name]()
    log("done")
