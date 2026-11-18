"""Step 10: does the unified atlas or the Allen CCF annotation fit the VMH cell types better?

Both atlases live in the same CCF space, and the Allen annotation has no VMH subdivisions, so the
comparison is of the nucleus as a whole. For each brain and each definition of the VMH
  Allen CCF VMH            the Allen annotation's VMH (structure 693)
  unified core             DM + C + VL + undivided ends of the unified atlas, nearest plate
  unified core + shell     the same plus the shell
neurons inside are compared with neurons of the two VMH subclasses, first at the released CCF
coordinates and then after the rigid shift that suits that definition best (same search as step 3),
so that neither atlas is penalised for the registration offset of the cells.

Reads   <data_root>/derived/anatomy.npz (Allen annotation around the VMH, from step 1) and the
        unified atlas plates
Writes  results/tables/atlas_comparison.csv   agreement with the cell types, per brain and definition
        results/tables/atlas_overlap.csv      how the Allen CCF VMH and the unified core overlap in CCF space
"""
import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import D, DATASETS, load_cells  # noqa: E402
from paths import RESULTS  # noqa: E402

ra = importlib.import_module("03_reference_atlas")


def masks():
    """The three VMH definitions as boolean volumes on the 10 um CCF grid, and that grid's origin (AP, DV, ML)."""
    z = np.load(D / "anatomy.npz", allow_pickle=True)
    origin = z["ccf_sub_origin"].astype(int)
    allen = z["ccf_sub"] == int(z["vmh_index"])
    V = ra.load_volume()
    ap = origin[0] + np.arange(allen.shape[0])
    plate = np.clip(np.round((ap - ra.AP0) / 10).astype(int) + 1, ra.K0, ra.K1)
    box = V[plate - ra.K0][:, origin[1]:origin[1] + allen.shape[1], origin[2]:origin[2] + allen.shape[2]]
    core = np.isin(box, [i for i, n in ra.VMH_IDS.items() if n in ra.CORE])
    whole = np.isin(box, list(ra.VMH_IDS))
    return {"Allen CCF VMH": allen, "unified core": core, "unified core + shell": whole}, origin


def inside(mask, origin, xyz, shift=(0.0, 0.0, 0.0)):
    idx = np.floor((xyz + np.asarray(shift)) / 0.01).astype(int) - origin
    ok = ((idx >= 0) & (idx < mask.shape)).all(1)
    out = np.zeros(len(xyz), dtype=bool)
    out[ok] = mask[idx[ok, 0], idx[ok, 1], idx[ok, 2]]
    return out


def score(ins, vt):
    tp = (ins & vt).sum()
    return 2 * tp / max(ins.sum() + vt.sum(), 1), tp / max(ins.sum(), 1), tp / max(vt.sum(), 1)


def best_shift(mask, origin, xyz, vt):
    best = (-1.0, (0.0, 0.0, 0.0))
    for dap in np.round(np.arange(-0.6, 0.601, 0.05), 2):
        for ddv in np.round(np.arange(-0.3, 0.301, 0.02), 2):
            for dml in np.round(np.arange(-0.2, 0.201, 0.02), 2):
                f1 = score(inside(mask, origin, xyz, (dap, ddv, dml)), vt)[0]
                if f1 > best[0]:
                    best = (f1, (float(dap) + 0.0, float(ddv) + 0.0, float(dml) + 0.0))
    return best[1]


def main():
    M, origin = masks()
    a, c = M["Allen CCF VMH"], M["unified core"]
    cen = lambda m: np.argwhere(m).mean(0)      # noqa: E731
    print(f"volume, mm3, both sides: Allen CCF VMH {a.sum() * 1e-6:.3f}, unified core {c.sum() * 1e-6:.3f}, "
          f"unified core + shell {M['unified core + shell'].sum() * 1e-6:.3f}")
    print(f"Allen CCF VMH against unified core: Dice {2 * (a & c).sum() / (a.sum() + c.sum()):.2f}; unified core centroid lies "
          f"{(cen(c) - cen(a))[0] * 10:.0f} um posterior, {(cen(c) - cen(a))[1] * 10:.0f} um ventral, of the Allen one; "
          f"{(a & M['unified core + shell']).sum() / a.sum():.0%} of the Allen VMH lies inside the unified core + shell")
    off = (cen(c) - cen(a)) * 10
    (RESULTS / "tables").mkdir(parents=True, exist_ok=True)
    pd.DataFrame([dict(allen_ccf_vmh_mm3=a.sum() * 1e-6, unified_core_mm3=c.sum() * 1e-6,
                       unified_core_and_shell_mm3=M["unified core + shell"].sum() * 1e-6,
                       dice_allen_vs_unified_core=2 * (a & c).sum() / (a.sum() + c.sum()),
                       share_of_allen_inside_unified_core_and_shell=(a & M["unified core + shell"]).sum() / a.sum(),
                       unified_core_centroid_posterior_um=off[0], unified_core_centroid_ventral_um=off[1],
                       unified_core_centroid_lateral_um=off[2])]).to_csv(RESULTS / "tables" / "atlas_overlap.csv", index=False, float_format="%.4g")
    rows = []
    for ds in DATASETS:
        cells = load_cells(ds)
        cells = cells[cells.neuron & cells.x_ccf.notna()]
        xyz = cells[["x_ccf", "y_ccf", "z_ccf"]].values.astype(float)
        vt = cells.vmh_type.values
        for name, mask in M.items():
            f1, pur, rec = score(inside(mask, origin, xyz), vt)
            sh = best_shift(mask, origin, xyz, vt)
            ins = inside(mask, origin, xyz, sh)
            g1, pur2, rec2 = score(ins, vt)
            rows.append(dict(dataset=ds, vmh_definition=name, n_vmh_type_neurons=int(vt.sum()),
                             f1_released=f1, share_vmh_type_released=pur, recall_released=rec,
                             shift_ap_mm=sh[0], shift_dv_mm=sh[1], shift_ml_mm=sh[2], n_inside_shifted=int(ins.sum()),
                             f1_shifted=g1, share_vmh_type_shifted=pur2, recall_shifted=rec2))
            print(f"{ds:10s} {name:22s} released: F1 {f1:.2f}, {pur:.0%} of neurons inside are VMH types | own best shift "
                  f"AP {sh[0]:+.2f} DV {sh[1]:+.2f} ML {sh[2]:+.2f}: F1 {g1:.2f}, {pur2:.0%}")
    pd.DataFrame(rows).to_csv(RESULTS / "tables" / "atlas_comparison.csv", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
