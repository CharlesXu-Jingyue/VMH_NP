"""Shared definitions for the spatial steps."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import data_root  # noqa: E402

D = data_root() / "derived"
NON_NEURONAL = "Astro|Oligo|OPC|Vascular|Immune|OEC"
VMH_SUBCLASSES = ("128 VMH Fezf1 Glut", "129 VMH Nr5a1 Glut")
# section-plane coordinate columns differ between the Allen and Zhuang releases
DATASETS = {
    "allen_m": dict(tag="merfish", label="Allen MERFISH, male", xy=("x_reconstructed", "y_reconstructed"), z="z_reconstructed"),
    "zhuang1_f": dict(tag="zhuang_abca_1", label="Zhuang MERFISH 1, female", xy=("x", "y"), z="z"),
    "zhuang2_m": dict(tag="zhuang_abca_2", label="Zhuang MERFISH 2, male", xy=("x", "y"), z="z"),
}


def load_cells(ds):
    """Cells of one spatial dataset with px, py (section plane, mm), neuron and vmh_type flags."""
    c = pd.read_parquet(D / f"{DATASETS[ds]['tag']}_cells.parquet")
    xk, yk = DATASETS[ds]["xy"]
    c = c.rename(columns={xk: "px", yk: "py"})
    c["neuron"] = ~c["class"].str.contains(NON_NEURONAL)
    c["vmh_type"] = c.subclass.isin(VMH_SUBCLASSES)
    c["in_ccf_vmh"] = c.parcellation_structure == "VMH"
    return c


# display colours shared by figures and the explorer (validated categorical order)
COMPARTMENT_COLOURS = {"VMHDM": "#2a78d6", "VMHC": "#eb6834", "VMHVL": "#1baf7a", "VMH": "#eda100", "VMHSh": "#e87ba4",
                       "outside": "#c9c8c2"}
COMPARTMENT_NAMES = {"VMHDM": "dorsomedial", "VMHC": "central", "VMHVL": "ventrolateral", "VMH": "undivided ends",
                     "VMHSh": "shell", "outside": "outside VMH"}
VMH_SUPERTYPES = ["0568 VMH Nr5a1 Glut_3", "0569 VMH Nr5a1 Glut_4", "0570 VMH Nr5a1 Glut_5", "0571 VMH Nr5a1 Glut_6",
                  "0566 VMH Nr5a1 Glut_1", "0567 VMH Nr5a1 Glut_2", "0564 VMH Fezf1 Glut_2", "0565 VMH Fezf1 Glut_3",
                  "0563 VMH Fezf1 Glut_1"]
