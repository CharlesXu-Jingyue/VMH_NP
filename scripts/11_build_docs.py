#!/usr/bin/env python3
"""Step 11: build the documentation shown inside the explorer (its Documentation tab).

docs/doc.js          method, data, how to read the explorer, main results, limits; loaded as a script
docs/figures/*.png   the figures it shows

The text is scripts/documentation_fragment.html. Every number that depends on the pipeline is
filled in here from results/tables, so the page cannot drift from the tables. Statements in the
template that are fixed text but rest on the results are asserted below; if a re-run makes one
false, this step stops and names it.
"""
import html
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import COMPARTMENTS as COMPS, DATASETS, SCREENED, VMH_SUPERTYPES  # noqa: E402
from paths import CONFIG, REPO, RESULTS  # noqa: E402

T, OUT = RESULTS / "tables", REPO / "docs"
REPO_URL = "https://github.com/CharlesXu-Jingyue/VMH_NP"
PANEL = {"allen_m": "500 target genes and 50 blank controls", "zhuang1_f": "1,122 genes", "zhuang2_m": "1,122 genes"}
FIGURES = ["atlas_interpolation.png", "atlas_alignment_zhuang1_f.png", "supertype_by_compartment.png",
           "projected_vs_measured.png", "compartment_ligands.png", "compartment_receptors.png"]


def td(v, cls=""):
    return f'<td class="{cls}">{v}</td>' if cls else f"<td>{v}</td>"


def tr(*cells):
    return "<tr>" + "".join(cells) + "</tr>"


def pct(x):
    return f"{100 * x:.0f}"


def genes(names):
    return ", ".join(f"<i>{html.escape(g)}</i>" for g in names)


def main():
    gb = pd.read_csv(T / "gene_by_compartment.csv")
    al = pd.read_csv(T / "atlas_alignment.csv").set_index("dataset")
    cmp_ = pd.read_csv(T / "atlas_comparison.csv")
    conc = pd.read_csv(T / "receptor_ligand_concordance.csv")
    sex = pd.read_csv(T / "sex_differences.csv")
    lig = pd.read_csv(CONFIG / "neuropeptide_ligands.csv").set_index("gene")
    rec = pd.read_csv(CONFIG / "neuropeptide_receptors.csv").set_index("gene")
    regions = pd.read_csv(CONFIG / "vmh_afferent_regions.csv")
    fill = dict(repo_url=REPO_URL, n_lig=len(lig), n_rec=len(rec), n_regions=len(regions))

    # panels and datasets
    meas = gb[gb.measured & (gb.compartment == "VMHDM")].groupby(["dataset", "kind"]).gene.nunique().unstack(fill_value=0)
    short = {"allen_m": "the Allen brain", "zhuang1_f": "Zhuang brain 1", "zhuang2_m": "Zhuang brain 2"}
    fill["panel_sentence"] = "; ".join(
        f"{meas.loc[ds, 'ligand']} precursors and {meas.loc[ds, 'receptor']} receptors in {short[ds]}" for ds in DATASETS) + "."
    rows = []
    for ds, meta in DATASETS.items():
        r = al.loc[ds]
        rows.append(tr(td(html.escape(meta["label"])),
                       td(f"MERFISH, {PANEL[ds]}; {int(r.n_core_aligned):,} neurons in the VMH core and {int(r.n_shell_aligned):,} in the shell"),
                       td(f"Cell positions; in situ measurement of {meas.loc[ds, 'ligand']} precursors and {meas.loc[ds, 'receptor']} receptors")))
    fill["dataset_rows"] = "\n".join(rows)

    # alignment
    fill["alignment_rows"] = "\n".join(
        tr(td(html.escape(DATASETS[ds]["label"])), td(f"{al.loc[ds, 'shift_ap_mm']:+.2f} mm", "n"), td(f"{al.loc[ds, 'shift_dv_mm']:+.2f} mm", "n"),
           td(f"{al.loc[ds, 'shift_ml_mm']:+.2f} mm", "n"), td(pct(al.loc[ds, "core_share_vmh_type_raw"]) + "%", "n"),
           td(pct(al.loc[ds, "core_share_vmh_type_aligned"]) + "%", "n")) for ds in DATASETS)
    rel = al.share_relabelled_by_interpolation
    fill["relabel_range"] = f"{pct(rel.min())} to {pct(rel.max())}%"
    assert (al.shift_ap_mm > 0).all(), "template says the atlas VMH sits posterior of the VMH cell types in every brain"
    assert ((al.shell_share_vmh_type_aligned > 0.4) & (al.shell_share_vmh_type_aligned < 0.7)).all(), "template says about half of shell neurons are VMH types"

    # the two atlases
    rows = []
    for ds in DATASETS:
        for _, r in cmp_[cmp_.dataset == ds].iterrows():
            rows.append(tr(td(html.escape(DATASETS[ds]["label"])), td(html.escape(r.vmh_definition)), td(f"{r.f1_released:.2f}", "n"),
                           td(pct(r.share_vmh_type_released) + "%", "n"), td(f"{r.shift_ap_mm:+.2f} mm", "n"), td(f"{r.f1_shifted:.2f}", "n"),
                           td(pct(r.share_vmh_type_shifted) + "%", "n")))
    fill["comparison_rows"] = "\n".join(rows)
    piv = cmp_.pivot(index="dataset", columns="vmh_definition")
    assert (piv["f1_shifted"].idxmax(axis=1) == "Allen CCF VMH").all(), "template says the Allen CCF VMH has the highest shifted F1 in every brain"
    assert (piv["share_vmh_type_shifted"]["Allen CCF VMH"] > piv["share_vmh_type_shifted"]["unified core"]).all(), \
        "template says more of the neurons inside the Allen CCF VMH are VMH types than inside the unified core"
    ov = pd.read_csv(T / "atlas_overlap.csv").iloc[0]
    fill["overlap_sentence"] = (f"they overlap with a Dice coefficient of {ov.dice_allen_vs_unified_core:.2f}, and the unified core is centred "
                                f"{ov.unified_core_centroid_posterior_um:.0f} µm posterior of the Allen VMH. "
                                f"{100 * ov.share_of_allen_inside_unified_core_and_shell:.0f}% of the Allen VMH lies inside the unified core plus shell.")
    assert abs(ov.allen_ccf_vmh_mm3 / ov.unified_core_mm3 - 1) < 0.1, "template says the Allen CCF VMH and the unified core have the same volume"
    assert ov.unified_core_centroid_posterior_um > 0, "template sentence says the unified core lies posterior of the Allen CCF VMH"

    # projected against measured (same selection as steps 6 and 8: the five atlas compartments)
    v = gb[gb.measured & gb.compartment.isin(COMPS) & gb.kind.isin(["ligand", "receptor"])]
    rho = {ds: v[v.dataset == ds][["proj_frac", "meas_frac"]].corr(method="spearman").iloc[0, 1] for ds in DATASETS}
    fill["spearman_sentence"] = ", ".join(f"{rho[ds]:.2f} ({short[ds]})" for ds in DATASETS)

    within = {}
    for ds, group in v.groupby("dataset"):
        correlations = [g.proj_frac.corr(g.meas_frac, method="spearman") for _, g in group.groupby("gene")
                        if g.proj_frac.nunique() > 1 and g.meas_frac.nunique() > 1]
        within[ds] = f"{pd.Series(correlations).median():.2f} ({len(correlations)} genes in {short[ds]})"
    fill["within_gene"] = ", ".join(within[ds] for ds in DATASETS)
    scores = pd.read_csv(T / "ligand_calibration.csv").set_index("evaluation")
    calibration = scores.loc["calibration"]
    fill["held_out_f1"] = f"{scores.loc['held-out clusters', 'f1']:.3f}"
    fill["calibration"] = f"F1 {calibration.f1:.3f} over {int(calibration.n_pairs):,} hypothalamic cluster–gene pairs"
    fill["n_tracing"] = int((regions.evidence_type == "mouse tracing citation").sum())
    # Use the Allen projection as reference columns after screening all brains and both methods.
    src = pd.read_csv(T / "ligand_sources.csv").set_index("gene")
    evidence = pd.concat([src.spatial_evidence[src.local], conc.set_index("receptor").query("in_vmh").spatial_evidence])
    a = gb[(gb.dataset == "allen_m") & gb.compartment.isin(COMPS)].pivot(index=["gene", "kind"], columns="compartment", values="proj_frac")[COMPS]
    hit = a.reset_index().set_index("gene").loc[evidence.index].copy()
    hit["top"] = hit[SCREENED].max(axis=1)
    # the screen's controls: what blank probes reach, how much of the undivided ends is VMH, and what was left out
    fl = pd.read_csv(T / "blank_probe_floor.csv")
    fl = fl[fl.compartment.isin(SCREENED)]
    assert set(fl.dataset) == {"allen_m"}, "template says only the Allen panel carries blank probes"
    fill["n_blank"] = int(fl.n_blank_probes.iloc[0])
    fill["blank_median"] = f"{100 * fl.blank_frac_median.min():.0f} to {100 * fl.blank_frac_median.max():.0f}%"
    fill["blank_floor"] = f"{100 * fl.blank_frac_p95.min():.0f} to {100 * fl.blank_frac_p95.max():.0f}%"
    fill["blank_reaching"] = f"{fl.n_blanks_reaching_detection_min.min()} to {fl.n_blanks_reaching_detection_min.max()}"
    cc = pd.read_csv(T / "compartment_composition.csv")
    vmh_share = cc[cc.supertype.isin(VMH_SUPERTYPES)].groupby(["dataset", "compartment"]).share.sum()
    fill["ends_share"] = ", ".join(f"{100 * vmh_share.get((ds, 'VMH'), 0):.0f}% ({short[ds]})" for ds in DATASETS)
    three = vmh_share[vmh_share.index.get_level_values("compartment").isin(["VMHDM", "VMHC", "VMHVL"])]
    fill["three_share"] = f"{100 * three.min():.0f} to {100 * three.max():.0f}%"
    assert vmh_share.xs("VMH", level="compartment").min() < 0.1 < three.min(), "template contrasts the ends with the three divided compartments"
    status = pd.concat([src.spatial_evidence, conc.set_index("receptor").spatial_evidence])
    weak, ends_only = status[status.str.startswith("Measured signal not above controls")].index, status[status.str.startswith("Undivided ends only")].index
    fill["weak_genes"], fill["ends_genes"] = genes(sorted(weak)), genes(sorted(ends_only))
    marg = pd.concat([src[["borderline", "local"]].rename(columns={"local": "listed"}),
                      conc.set_index("receptor")[["borderline", "in_vmh"]].rename(columns={"in_vmh": "listed"})])
    fill["borderline_listed"] = genes(sorted(marg.index[marg.borderline & marg.listed])) or "none"
    fill["borderline_unlisted"] = genes(sorted(marg.index[marg.borderline & ~marg.listed])) or "none"
    alone = sorted(status[status.str.startswith("Measured only")].index)
    fill["measured_alone"] = (f"{genes(alone)} rest{'s' if len(alone) == 1 else ''} on a measured value alone" if alone
                              else "No listed gene rests on a measured value alone")
    # whether the brains without blank probes see the same thing
    other = gb[gb.gene.isin(alone) & (gb.dataset != "allen_m") & gb.measured & gb.compartment.isin(SCREENED)]
    fill["measured_alone_elsewhere"] = ("in the Zhuang brains the highest measured share of a compartment's neurons is "
                                        + ", ".join(f"{100 * d.meas_frac.max():.0f}% for <i>{g}</i>" for g, d in other.groupby("gene"))
                                        if len(other) else "none of them is on the Zhuang panels")
    L, R = hit[hit.kind == "ligand"].sort_values("top", ascending=False), hit[hit.kind == "receptor"].copy()
    fill["n_lig_hit"], fill["n_rec_hit"] = len(L), len(R)

    def peptide(g):
        """Named peptides of a precursor; where the gene table lists none (some granins), the protein's name from its note."""
        v = lig.loc[g, "peptides"]
        if isinstance(v, str) and v.strip():
            return v
        note = lig.loc[g, "note"]
        return note.split(";")[0].strip() if isinstance(note, str) and note.strip() else str(lig.loc[g, "family"]).lower() + " precursor"
    fill["ligand_rows"] = "\n".join(
        tr(td(f"<i>{g}</i>"), td(html.escape(peptide(g))), td(evidence[g]), *[td(pct(r[c]), "n") for c in COMPS]) for g, r in L.iterrows())

    def lean(r):
        dm, vl = r.VMHDM, r.VMHVL
        if abs(dm - vl) >= 0.05 and max(dm, vl) >= 1.5 * min(dm, vl):
            return "dorsomedial" if dm > vl else "ventrolateral"
        return ""
    R["leans"] = R.apply(lean, axis=1)
    R["order"] = R.leans.map({"": 0, "dorsomedial": 1, "ventrolateral": 2})
    R = R.sort_values(["order", "top"], ascending=[True, False])
    fill["receptor_rows"] = "\n".join(
        tr(td(f"<i>{g}</i>"), td(html.escape(str(rec.loc[g, "family"]))), td(r.leans or "no"), td(evidence[g]), *[td(pct(r[c]), "n") for c in COMPS]) for g, r in R.iterrows())

    # sex
    fill["sex_sentence"] = (f"{len(sex)} supertype-by-gene tests in {sex.supertype.nunique()} supertypes, "
                            f"lowest adjusted p-value {sex.q_value.min():.2f}")
    assert sex.q_value.min() > 0.05, "template says no sex difference survives correction"

    # ligand sources
    inv = conc[conc.in_vmh]
    fill["source_rows"] = "\n".join(
        tr(td(html.escape(verdict.capitalize())), td(str(len(group)), "n"),
           td(genes(group.sort_values("peak_fraction", ascending=False).receptor)))
        for verdict, group in inv.groupby("ligand_sources"))

    page = (HERE / "documentation_fragment.html").read_text()
    page = page[page.index("-->") + 3:].strip()                # drop the note to maintainers at the top
    for k, val in fill.items():
        page = page.replace("{{" + k + "}}", str(val))
    assert "{{" not in page, "unfilled placeholder: " + page[page.index("{{"):page.index("{{") + 40]
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    for f in FIGURES:
        shutil.copyfile(RESULTS / "figures" / f, OUT / "figures" / f)
    (OUT / "doc.js").write_text("window.__DOC__ = " + json.dumps(page).replace("</", "<\\/") + ";\n")
    print(f"wrote {OUT / 'doc.js'} ({len(page) / 1e3:.0f} kB) and {len(FIGURES)} figures; "
          f"{len(L)} precursors and {len(R)} receptors in the results tables")


if __name__ == "__main__":
    main()
