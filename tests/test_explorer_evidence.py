"""Regression checks for the published explorer. The map test uses the configured data extracts."""
import base64
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from common import COMPARTMENTS, D, DATASETS, load_cells

TABLES = ROOT / "results" / "tables"


class ExplorerEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        js = (ROOT / "docs" / "data.js").read_text()
        cls.site = json.loads(js.split(" = ", 1)[1].rstrip(";\n"))
        cls.rows = pd.read_csv(TABLES / "gene_by_compartment.csv")

    def test_all_compartments_and_paired_calls(self):
        self.assertEqual(self.site["comps"], [*COMPARTMENTS, "core"])
        for gene, rows in self.rows[self.rows.compartment.isin(COMPARTMENTS)].groupby("gene"):
            card = self.site["genes"][gene]
            positive = rows.projected_hit | rows.measured_hit.eq(True)
            if "src" in card:
                self.assertEqual(card["src"]["hit"], positive.any(), gene)
            for ds, group in rows.groupby("dataset"):
                self.assertEqual(set(group.compartment), set(COMPARTMENTS))
                disagree = group.measured & (group.projected_hit != group.measured_hit)
                expected = group.loc[disagree].set_index("compartment").evidence.to_dict()
                self.assertEqual(card["discordant"][ds], expected, gene)
                self.assertEqual(card["meas"][ds] is None, not group.measured.any(), gene)

    def test_review_cases(self):
        r = pd.read_csv(TABLES / "receptor_ligand_concordance.csv").set_index("receptor")
        for gene in ("Npffr2", "Crhr1", "Npsr1"):
            self.assertGreaterEqual(r.loc[gene, "vmh_VMH"], 0.1)
            self.assertTrue(r.loc[gene, "in_vmh"])
        self.assertEqual(r.loc["Rxfp1", "spatial_evidence"], "Measured only")
        self.assertFalse(r.loc[["Insr", "Igf1r"], "in_vmh"].any())
        self.assertIn("GLP-1 (Gcg, intestine)", r.loc["Glp1r", "endocrine_ligands"])
        self.assertIn("same-cell Ramp1", r.loc["Calcrl", "conditional_partners"])
        self.assertTrue(pd.isna(r.loc["Calcrl", "afferent_ligands"]))

    def test_measured_flags_and_registration(self):
        labels = pd.read_parquet(D / "atlas_labels.parquet")
        for ds, meta in DATASETS.items():
            brain = self.site["brains"][ds]
            cells = load_cells(ds).merge(labels[labels.dataset == ds], on="cell_label", validate="one_to_one")
            measured = np.load(D / (meta["tag"] + "_expr.npz"), allow_pickle=True)
            np.testing.assert_array_equal(cells.cell_label, measured["cell_label"])
            signal = measured["X"]
            panel = list(measured["genes"])
            inside = cells[cells.kim_aligned.isin(COMPARTMENTS)]
            in_view = (cells.neuron & cells.px.between(inside.px.min() - .35, inside.px.max() + .35)
                       & cells.py.between(inside.py.min() - .35, inside.py.max() + .35))
            ids = np.concatenate([cells.index[in_view & (cells.brain_section_label == section["id"])]
                                  for section in brain["sections"]])
            codes = np.frombuffer(base64.b64decode(brain["cp"]), dtype=np.uint8)
            np.testing.assert_array_equal(codes == self.site["compCode"]["no_ccf"], cells.iloc[ids].x_ccf.isna())
            for gene, encoded in brain["measured"].items():
                values = signal[:, panel.index(gene)].astype(float)
                cutoff = np.percentile(values[~cells.neuron], 95)
                bits = np.frombuffer(base64.b64decode(encoded["above"]), dtype=np.uint8)
                flags = np.unpackbits(bits, bitorder="little", count=len(ids))
                np.testing.assert_array_equal(flags, values[ids] > cutoff, err_msg=f"{ds} {gene}")


if __name__ == "__main__":
    unittest.main()
