#!/usr/bin/env python3
"""Step 3: assign every spatial cell to the unified mouse brain atlas of Chon et al. 2019
(Franklin-Paxinos labels in Allen CCF space; Dryad doi:10.5061/dryad.t1g1jwsxw), which
divides the VMH into dorsomedial (VMHDM), central (VMHC), ventrolateral (VMHVL) and shell
(VMHSh) compartments. This replaces the hand-drawn boundaries.

Input   <data_root>/derived/reference/kim_unified_atlas/1_Anatomical_Labels.zip
        123 coronal label images, 10 um in plane, 100 um apart, Allen structure ids
        config/kim_unified_structures.csv (id -> acronym)

The atlas is a stack of plates 100 um apart, each drawn on its own, and tissue sections are
not parallel to the plates: within one section the VMH spans 50 to 220 um along the
anterior-posterior axis. Reading the nearest plate therefore stitches two or three different
drawings together in one section, with seams. To avoid that, the VMH labels are interpolated
between plates (shape-based interpolation): each plate's compartments become signed distance
fields (positive inside), lightly smoothed in the plane (SMOOTH_PX), and a monotone cubic
through the plates gives a continuous field at any position. It interpolates smoothed plate fields; smoothing can shift boundary pixels.

Three labels are stored per cell:
  kim_raw               nearest plate at the cell's released CCF coordinates
  kim_aligned           interpolated atlas after a rigid shift of the cell coordinates (default)
  kim_aligned_nearest   nearest plate after the same shift, for comparison

Why a shift: in all three brains the VMH cell types are registered 0.1 to 0.4 mm anterior
of the atlas VMH, so raw labels put most Esr1-rich neurons in the shell or outside. The
shift (three numbers per brain: AP, DV, ML) maximises F1 between "neuron of a VMH subclass"
and "inside the interpolated atlas VMH, any compartment". Compartment identity is not used in the fit, and
the resulting compartment composition describes this molecular fit
(results/tables/atlas_supertype_by_compartment.csv).

Slice position: slice k of the atlas is CCF 10 um plane 10*(k-1) + AP0, planes counted from
zero, with AP0 = 94. This is exact: the atlas ships the CCF template images its labels were
drawn on (2_AllenCCF_background.zip), and all 123 are identical, pixel for pixel, to those
planes of the Allen template. Run with --verify-ap to repeat (needs that archive beside the
labels and average_template_10.nii.gz from the Allen Brain Cell Atlas release). The file names
of the plates give bregma positions, which do not fix a CCF plane, and an earlier estimate
from structures shared with the Allen annotation had put AP0 at 93.

Outputs
  <data_root>/derived/atlas_labels.parquet          per cell: raw, aligned and aligned-nearest label
  <data_root>/derived/atlas_outlines.json           aligned compartment outlines per section
  results/tables/atlas_alignment.csv                shift and agreement per brain
  results/tables/atlas_supertype_by_compartment.csv
"""
import argparse
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from scipy.interpolate import LinearNDInterpolator, PchipInterpolator
from skimage import measure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import D, DATASETS, load_cells  # noqa: E402
from paths import CONFIG, RESULTS, abc_root  # noqa: E402

ZIP = D / "reference" / "kim_unified_atlas" / "1_Anatomical_Labels.zip"
AP0 = 94                     # CCF 10 um AP plane (from zero) of atlas slice 1; verified against the template
K0, K1 = 40, 80              # slices loaded (AP about 4.8 to 8.9 mm), ample around the hypothalamus
VMH_IDS = {777: "VMHDM", 769: "VMHC", 785: "VMHVL", 693: "VMH", 2044: "VMHSh"}
CORE = ("VMHDM", "VMHC", "VMHVL", "VMH")     # "VMH" = undivided nucleus at its anterior and posterior ends
OUTLINE = ("VMHDM", "VMHC", "VMHVL", "VMH", "VMHSh")
SUB = ("VMHDM", "VMHC", "VMHVL")             # the divided core
ROWS, COLS = (560, 780), (400, 740)          # 10 um crop of each plate around the VMH (DV rows, ML columns), both sides
KA, KB = 50, 71                              # plates held as fields; the VMH occupies 55 to 66
SMOOTH_PX = 1.5                              # in-plane Gaussian smoothing of each field, in 10 um pixels


def load_volume():
    z = zipfile.ZipFile(ZIP)
    names = {int(n.rsplit("/", 1)[-1].split("_")[0]): n for n in z.namelist() if n.endswith(".tif")}
    return np.stack([tifffile.imread(io.BytesIO(z.read(names[k]))) for k in range(K0, K1 + 1)])   # slice, DV, ML


def lookup(V, xyz, shift=(0.0, 0.0, 0.0)):
    """Atlas id at CCF coordinates xyz (mm; AP, DV, ML) after adding shift. 0 = outside the loaded block."""
    k = np.round((np.floor((xyz[:, 0] + shift[0]) / 0.01) - AP0) / 10) + 1
    r = np.floor((xyz[:, 1] + shift[1]) / 0.01)
    q = np.floor((xyz[:, 2] + shift[2]) / 0.01)
    ok = np.isfinite(k) & np.isfinite(r) & np.isfinite(q)
    ok &= (k >= K0) & (k <= K1) & (r >= 0) & (r < V.shape[1]) & (q >= 0) & (q < V.shape[2])
    lab = np.zeros(len(xyz), dtype=np.int32)
    lab[ok] = V[(k[ok] - K0).astype(int), r[ok].astype(int), q[ok].astype(int)]
    return lab


def _sdf(m):
    """Signed distance to the boundary of mask m in pixels, positive inside, zero on the pixel edge."""
    return np.where(m, ndi.distance_transform_edt(m) - 0.5, 0.5 - ndi.distance_transform_edt(~m))


class Field:
    """A stack of plate images made continuous: cubic B-spline within each plate, and between plates
    a monotone cubic (PCHIP) through the plate values at every in-plane position. The monotone cubic
    is smooth across plates and never leaves the range of the two neighbouring plates, so a position
    lies inside a structure only if at least one neighbouring plate puts it there."""

    def __init__(self, stack):
        self.pchip = PchipInterpolator(np.arange(len(stack)), stack, axis=0)
        slope = self.pchip.derivative()(np.arange(len(stack)))
        self.n = len(stack)
        self.value = ndi.spline_filter(stack, order=3)
        self.slope = ndi.spline_filter(slope, order=3)

    def on_grid(self, k):
        """Values at plate positions k on the plates' own pixel grid (no in-plane interpolation needed)."""
        return self.pchip(np.clip(k, 0, self.n - 1))

    def __call__(self, k, r, q):
        """Value at plate position k (index into the stack, clamped to its ends) and pixel (r, q)."""
        k = np.clip(k, 0, self.n - 1)
        k0 = np.minimum(np.floor(k), self.n - 2)
        t = k - k0
        at = lambda F, kk: ndi.map_coordinates(F, [kk, r, q], order=3, prefilter=False)       # noqa: E731
        f0, f1, d0, d1 = at(self.value, k0), at(self.value, k0 + 1), at(self.slope, k0), at(self.slope, k0 + 1)
        return ((2 * t ** 3 - 3 * t ** 2 + 1) * f0 + (t ** 3 - 2 * t ** 2 + t) * d0
                + (-2 * t ** 3 + 3 * t ** 2) * f1 + (t ** 3 - t ** 2) * d1)


class SmoothAtlas:
    """The VMH compartments as continuous fields (see Field), in pixels of 10 um, positive inside.

    whole   signed distance to the VMH including its shell
    core    signed distance to the VMH without the shell
    sub     one field per divided compartment (DM, C, VL), over the plates that have them
    Where a plate lacks the structure (beyond either end of the nucleus), the field continues the
    nearest plate's, lowered by twice its largest value per plate of distance, so the nucleus
    shrinks to a point halfway to the first plate without it. The plates at the two ends of the
    nucleus are undivided; the switch between undivided and divided is made halfway between the
    last undivided and the first divided plate, because that difference is one of labelling.
    Inside the core, a position belongs to the divided compartment whose field is largest."""

    def __init__(self, V):
        crop = V[KA - K0:KB - K0 + 1, ROWS[0]:ROWS[1], COLS[0]:COLS[1]]
        ks = np.arange(KA, KB + 1)
        core_ids = [i for i, n in VMH_IDS.items() if n in CORE]
        self.whole = Field(self._stack(np.isin(crop, list(VMH_IDS))))
        self.core = Field(self._stack(np.isin(crop, core_ids)))
        div = [int(k) for k, pl in zip(ks, crop) if all((pl == i).any() for i, n in VMH_IDS.items() if n in SUB)]
        assert div == list(range(div[0], div[-1] + 1)), "divided plates are not contiguous"
        self.d0, self.d1 = div[0], div[-1]
        self.sub = [Field(np.stack([ndi.gaussian_filter(_sdf(crop[k - KA] == i), SMOOTH_PX) for k in div]))
                    for i, n in VMH_IDS.items() if n in SUB]
        self.names = ["VMHSh", "VMH"] + [n for n in VMH_IDS.values() if n in SUB]

    @staticmethod
    def _stack(masks):
        have = [j for j, m in enumerate(masks) if m.any()]
        f = {j: ndi.gaussian_filter(_sdf(masks[j]), SMOOTH_PX) for j in have}
        out = []
        for j in range(len(masks)):
            if j in f:
                out.append(f[j])
            else:
                n = min(have, key=lambda h: abs(h - j))
                out.append(f[n] - 2 * f[n].max() * abs(n - j))
        return np.stack(out)

    @staticmethod
    def _index(xyz, shift):
        """Plate position and pixel position (within the crop) of CCF coordinates in mm, and which are usable."""
        kf = ((xyz[:, 0] + shift[0]) / 0.01 - 0.5 - AP0) / 10 + 1
        r = (xyz[:, 1] + shift[1]) / 0.01 - 0.5 - ROWS[0]
        q = (xyz[:, 2] + shift[2]) / 0.01 - 0.5 - COLS[0]
        ok = np.isfinite(kf) & np.isfinite(r) & np.isfinite(q)
        ok &= (kf >= KA) & (kf <= KB) & (r >= 0) & (r <= ROWS[1] - ROWS[0] - 1) & (q >= 0) & (q <= COLS[1] - COLS[0] - 1)
        return kf, r, q, ok

    def margins(self, xyz, shift=(0.0, 0.0, 0.0)):
        """Per compartment (rows as self.names), a field that is positive inside it; pixels as units."""
        kf, r, q, ok = self._index(xyz, shift)
        M = np.full((len(self.names), len(xyz)), -1e3)
        if not ok.any():
            return M
        whole = self.whole(kf[ok] - KA, r[ok], q[ok])
        M[0, ok] = whole
        ok[ok] = whole > -3                 # the other fields are only needed in and just around the VMH
        if not ok.any():
            return M
        kf, r, q, whole = kf[ok], r[ok], q[ok], whole[whole > -3]
        core = self.core(kf - KA, r, q)
        sub = np.stack([F(kf - self.d0, r, q) for F in self.sub])
        divided = np.minimum(kf - (self.d0 - 0.5), (self.d1 + 0.5) - kf) * 10      # > 0 between the divided plates
        M[0, ok] = np.minimum(whole, -core)
        M[1, ok] = np.minimum(core, -divided)
        for j in range(len(sub)):
            M[2 + j, ok] = np.minimum(np.minimum(core, divided), sub[j] - np.delete(sub, j, axis=0).max(0))
        return M

    def label(self, xyz, shift=(0.0, 0.0, 0.0)):
        M = self.margins(xyz, shift)
        return np.where(M.max(0) > 0, np.array(self.names, dtype=object)[M.argmax(0)], "outside")


def agreement(inside, vmh_type):
    tp = (inside & vmh_type).sum()
    return 2 * tp / max(inside.sum() + vmh_type.sum(), 1), tp / max(vmh_type.sum(), 1), tp / max(inside.sum(), 1)


def fit_shift(A, xyz, vmh_type):
    """Rigid shift (AP, DV, ML; mm) of the cells that maximises F1 with the interpolated atlas VMH.

    The search tries about 16,000 shifts. To keep it quick, the interpolated field is sampled once on
    the 10 um CCF grid and read linearly from there, and only neurons within 1 mm of the VMH in
    every direction are evaluated (no shift in the search range can bring any other neuron inside)."""
    a0 = 10 * (KA - 1) + AP0                                      # CCF AP voxel of plate KA
    ap = np.arange(a0, 10 * (KB - 1) + AP0 + 1)
    dense = A.whole.on_grid((ap - AP0) / 10 + 1 - KA).astype(np.float32)
    box = np.argwhere(dense > 0)
    lo, hi = box.min(0) - 100, box.max(0) + 100
    idx = np.c_[xyz[:, 0] / 0.01 - 0.5 - a0, xyz[:, 1] / 0.01 - 0.5 - ROWS[0], xyz[:, 2] / 0.01 - 0.5 - COLS[0]]
    near = ((idx >= lo) & (idx <= hi)).all(1)
    idx, vt, n_type = idx[near].T, vmh_type[near], vmh_type.sum()
    best = (-1.0, (0.0, 0.0, 0.0))
    for dap in np.round(np.arange(-0.6, 0.601, 0.05), 2):
        for ddv in np.round(np.arange(-0.3, 0.301, 0.02), 2):
            for dml in np.round(np.arange(-0.2, 0.201, 0.02), 2):
                at = [idx[0] + dap / 0.01, idx[1] + ddv / 0.01, idx[2] + dml / 0.01]
                ins = ndi.map_coordinates(dense, at, order=1, mode="constant", cval=-1e3) > 0
                f1 = 2 * (ins & vt).sum() / max(ins.sum() + n_type, 1)
                if f1 > best[0]:
                    best = (f1, (float(dap) + 0.0, float(ddv) + 0.0, float(dml) + 0.0))
    return best[1]


def section_outlines(A, cells, shift):
    """Outlines of the aligned, interpolated VMH compartments in each section's own plane coordinates.

    Section coordinates map to CCF through the cells themselves. Each compartment's outline is the
    zero contour of its continuous field, so it is traced at sub-pixel precision and agrees with the
    cell labels. The reported residual is that of a plain per-section affine, as a measure of how
    nonlinear the released registration is."""
    out, resid = {}, []
    for sec, d in cells.groupby("brain_section_label"):
        d = d[d.x_ccf.notna()]
        if len(d) < 50:
            continue
        P = np.c_[d.px.values, d.py.values, np.ones(len(d))]
        B, *_ = np.linalg.lstsq(P, d[["x_ccf", "y_ccf", "z_ccf"]].values, rcond=None)
        resid.append(np.sqrt(((P @ B - d[["x_ccf", "y_ccf", "z_ccf"]].values) ** 2).sum(1)))
        xs = np.arange(d.px.min(), d.px.max(), 0.01)
        ys = np.arange(d.py.min(), d.py.max(), 0.01)
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        # CCF position of every grid point: interpolate between the cells (registration is locally
        # nonlinear in the Zhuang brains), falling back to the affine outside their convex hull
        grid = np.c_[gx.ravel(), gy.ravel()]
        ccf = LinearNDInterpolator(P[:, :2], d[["x_ccf", "y_ccf", "z_ccf"]].values)(grid)
        miss = ~np.isfinite(ccf).all(1)
        ccf[miss] = np.c_[grid[miss], np.ones(miss.sum())] @ B
        M = A.margins(ccf, shift)
        polys = {}
        for name, m in zip(A.names, M):
            m = m.reshape(gx.shape)
            if (m > 0).sum() < 30:
                continue
            lines = []
            for c in measure.find_contours(np.pad(m, 1, constant_values=-1e3), 0.0):
                c = measure.approximate_polygon(c - 1, tolerance=0.1)
                if len(c) >= 8:
                    lines.append([[round(float(xs[0] + p[0] * 0.01), 4), round(float(ys[0] + p[1] * 0.01), 4)] for p in c])
            if lines:
                polys[name] = lines
        out[sec] = polys
    r = np.concatenate(resid)
    return out, float(np.median(r)), float(np.percentile(r, 95))


def verify_ap():
    """Check AP0 against the Allen template itself. The atlas ships the CCF template images its labels were
    drawn on (2_AllenCCF_background.zip, beside the labels); each must be identical to exactly one plane of
    image_volumes/Allen-CCF-2020/*/average_template_10.nii.gz. Needs about 3 GB of memory."""
    import nibabel as nib
    tpl = sorted((abc_root() / "image_volumes/Allen-CCF-2020").glob("*/average_template_10.nii.gz"))
    bg = ZIP.with_name("2_AllenCCF_background.zip")
    if not tpl or not bg.exists():
        raise SystemExit(f"needs {bg} and average_template_10.nii.gz under image_volumes/Allen-CCF-2020/")
    T = np.asanyarray(nib.load(tpl[0]).dataobj)                      # AP, DV, ML at 10 um
    z = zipfile.ZipFile(bg)
    names = sorted(n for n in z.namelist() if n.endswith(".tif"))    # AllenCCF_Z001 ... : image k belongs to label plate k
    good = 0
    for k, n in enumerate(names, start=1):
        img = tifffile.imread(io.BytesIO(z.read(n)))
        want = 10 * (k - 1) + AP0
        same = [a for a in range(max(want - 6, 0), min(want + 7, T.shape[0])) if np.array_equal(T[a], img)]
        good += same == [want]
        if same != [want]:
            print(f"image {k}: identical to template plane(s) {same}, expected {want}")
    print(f"{good} of {len(names)} background images are identical to template plane 10*(k-1) + {AP0} and to no plane near it")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verify-ap", action="store_true", help="check the slice position constant against the Allen template and exit")
    args = ap.parse_args()
    if args.verify_ap:
        return verify_ap()
    V = load_volume()
    A = SmoothAtlas(V)
    st = pd.read_csv(CONFIG / "kim_unified_structures.csv")
    acr = dict(zip(st.id, st.acronym))
    names_of = lambda ids: pd.Series(ids).map(acr).fillna("none").values      # noqa: E731
    (RESULTS / "tables").mkdir(parents=True, exist_ok=True)
    labels, rows, outlines, tabs = [], [], {}, []
    for ds in DATASETS:
        c = load_cells(ds)
        xyz = c[["x_ccf", "y_ccf", "z_ccf"]].values.astype(float)
        has = np.isfinite(xyz).all(1)
        fitm = (c.neuron & has).values
        vt = c.vmh_type.values
        shift = fit_shift(A, xyz[fitm], vt[fitm])
        raw, near = lookup(V, xyz), lookup(V, xyz, shift)
        smooth = np.full(len(c), "none", dtype=object)
        smooth[has] = A.label(xyz[has], shift)
        # outside the VMH the interpolated atlas has no fields; keep the nearest plate's name there
        ali = np.where(smooth == "outside", np.where(np.isin(near, list(VMH_IDS)), "none", names_of(near)), smooth)
        out = pd.DataFrame({"dataset": ds, "cell_label": c.cell_label.values, "kim_raw": names_of(raw),
                            "kim_aligned": ali, "kim_aligned_nearest": names_of(near)})
        out.loc[~has, ["kim_raw", "kim_aligned", "kim_aligned_nearest"]] = "no_ccf"
        labels.append(out)
        n = c.neuron.values & has
        row = dict(dataset=ds, sex=c.donor_sex.iloc[0], n_neurons=int(n.sum()), n_vmh_type_neurons=int((n & vt).sum()),
                   shift_ap_mm=shift[0], shift_dv_mm=shift[1], shift_ml_mm=shift[2])
        comp = {"raw": pd.Series(raw).map(VMH_IDS).fillna("outside").values,
                "aligned": np.where(np.isin(smooth, list(VMH_IDS.values())), smooth, "outside"),
                "aligned_nearest": pd.Series(near).map(VMH_IDS).fillna("outside").values}
        for name, a in comp.items():
            a = a[n]
            f1, rec, prec = agreement(a != "outside", vt[n])
            core, shell = np.isin(a, CORE), a == "VMHSh"
            row.update({f"f1_{name}": f1, f"recall_{name}": rec, f"precision_{name}": prec,
                        f"n_core_{name}": int(core.sum()), f"core_share_vmh_type_{name}": float(vt[n][core].mean()),
                        f"n_shell_{name}": int(shell.sum()), f"shell_share_vmh_type_{name}": float(vt[n][shell].mean())})
            t = (pd.crosstab(c.supertype.values[n & vt], comp[name][n & vt])
                   .reindex(columns=list(VMH_IDS.values()) + ["outside"], fill_value=0))
            t.insert(0, "labels", name)
            t.insert(0, "dataset", ds)
            tabs.append(t.rename_axis("supertype").reset_index())
        both = n & ((comp["aligned"] != "outside") | (comp["aligned_nearest"] != "outside"))
        row.update(n_in_vmh_either_rule=int(both.sum()),
                   share_relabelled_by_interpolation=float((comp["aligned"][both] != comp["aligned_nearest"][both]).mean()))
        outlines[ds], med, p95 = section_outlines(A, c, shift)
        row.update(section_affine_residual_median_mm=med, section_affine_residual_p95_mm=p95)
        rows.append(row)
        print(f"{ds}: shift AP {shift[0]:+.2f} DV {shift[1]:+.2f} ML {shift[2]:+.2f} mm; F1 {row['f1_raw']:.2f} -> {row['f1_aligned']:.2f}; "
              f"core neurons that are VMH types {row['core_share_vmh_type_raw']:.3f} -> {row['core_share_vmh_type_aligned']:.3f} "
              f"(nearest plate, same shift: {row['core_share_vmh_type_aligned_nearest']:.3f}); "
              f"{100 * row['share_relabelled_by_interpolation']:.1f}% of VMH neurons relabelled by interpolation; "
              f"section affine residual median {med * 1000:.0f} um, 95th pct {p95 * 1000:.0f} um")
    pd.concat(labels, ignore_index=True).to_parquet(D / "atlas_labels.parquet", index=False)
    pd.DataFrame(rows).to_csv(RESULTS / "tables" / "atlas_alignment.csv", index=False, float_format="%.4g")
    pd.concat(tabs, ignore_index=True).to_csv(RESULTS / "tables" / "atlas_supertype_by_compartment.csv", index=False)
    with open(D / "atlas_outlines.json", "w") as fh:
        json.dump(dict(units="mm, section coordinates of each dataset", labels="aligned, interpolated between plates",
                       datasets=outlines), fh)


if __name__ == "__main__":
    main()
