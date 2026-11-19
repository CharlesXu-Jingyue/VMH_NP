#!/usr/bin/env python3
"""Step 7 (supplemental part 3): which peptides could be released in the VMH, and do they
match the receptors found there?

Three sources are considered for each ligand gene:
  local      measured or projected transcript detection reaches 10% in any brain and compartment
  afferent   neurons in a region reported to project to VMH make it. Regions are listed in
             config/vmh_afferent_regions.csv, curated by hand from the tracing literature and
             meant to be edited.
  endocrine  the peptide has a curated peripheral or pituitary source

This is transcript-level evidence that a peptide is available, not evidence of release, and
"afferent" is not projection-specific: it asks whether a region that projects to VMH contains
neurons making the peptide, not whether those neurons are the ones that project.

A cluster "expresses" a gene when its mean log2(CPM+1) in the Allen whole-brain cluster table
is at least max(THETA, floor of the gene). The floor is the 99th percentile reached by
hypothalamic clusters that step 2 called background for that gene, and THETA is chosen to
best reproduce the broad and subset calls of step 2 in the hypothalamus. Genes whose background
in non-neuronal cells is at least AMBIENT_PRONE_CPM are held to AMBIENT_PRONE_THETA.

A region counts as a source when at least AFFERENT_MIN of its neurons belong to expressing
clusters, or, for small specialised populations (oxytocin, orexin, ...), when it holds at least
SPECIALISED_SHARE of all neurons of the expressing clusters brain-wide and at least
SPECIALISED_MIN_CELLS of them.

Outputs (results/tables)
  ligand_by_region.csv              every ligand x every atlas structure: neurons in expressing clusters, as a
                                    share of the structure and of all such neurons brain-wide
  ligand_sources.csv                per ligand: local, afferent and endocrine evidence
  receptor_ligand_concordance.csv   per receptor in VMH: which sources could supply a cognate ligand
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import COMPARTMENTS, D, DATASETS, NON_NEURONAL  # noqa: E402
from paths import CONFIG, RESULTS  # noqa: E402

T = RESULTS / "tables"
AFFERENT_MIN, MIN_REGION_NEURONS = 0.10, 200
SPECIALISED_SHARE, SPECIALISED_MIN_CELLS = 0.10, 20
AMBIENT_PRONE_CPM, AMBIENT_PRONE_THETA = 2.0, 4.0



def classifier_scores(calls, pred):
    truth = calls.call.isin(["broad", "subset"])
    tp = int((pred & truth).sum())
    return dict(n_pairs=len(calls), precision=tp / max(pred.sum(), 1),
                recall=tp / max(truth.sum(), 1), f1=2 * tp / max(pred.sum() + truth.sum(), 1))


def fit_threshold(calls, ligands, prone):
    floor = calls[calls.call == "background"].groupby("gene").val.quantile(0.99).reindex(ligands).fillna(0.0)
    floor = floor.where(~prone, np.maximum(floor, AMBIENT_PRONE_THETA))
    limits = calls.gene.map(floor)
    def score(theta):
        return classifier_scores(calls, calls.val >= np.maximum(theta, limits))["f1"]
    theta = max(np.round(np.arange(0.5, 4.01, 0.1), 1), key=score)
    return theta, np.maximum(theta, floor)


def cluster_expression(lig):
    """Boolean cluster x ligand matrix (index: cluster name) and the threshold used."""
    means = pd.read_parquet(D / "cluster_means.parquet")
    tax = pd.read_parquet(D / "cluster_taxonomy.parquet")
    tax["cluster_label"] = "CS20230722_CLUS_" + tax.cluster.str[:4]
    means = means.merge(tax[["cluster_label", "cluster", "cluster_alias", "class"]], on="cluster_label", how="inner")
    assert len(means) == len(tax), "cluster labels did not map one to one"
    M = means.set_index("cluster")[lig].astype(float)
    bg = pd.read_csv(T / "gene_background.csv")
    bg = bg[bg.chemistry == "10Xv3"].set_index("gene").background_cpm
    # calibrate against the hypothalamic per-cluster calls of step 2
    e = pd.read_parquet(T / "expr_by_type.parquet")
    for c in ("level", "type", "chemistry", "sex", "gene", "call"):
        e[c] = e[c].astype(str)
    e = e[(e.level == "cluster") & (e.chemistry == "10Xv3") & (e.sex == "all") & e.gene.isin(lig) & (e.n_cells >= 50)]
    e = e[e.type.isin(M.index)].copy()
    e["val"] = M.values[M.index.get_indexer(e.type), M.columns.get_indexer(e.gene)]
    prone = bg.reindex(lig).fillna(0) >= AMBIENT_PRONE_CPM
    theta, thr = fit_threshold(e, lig, prone)
    pred = e.val >= e.gene.map(thr)
    report = [dict(evaluation="calibration", theta=theta, **classifier_scores(e, pred))]
    # Hold out whole clusters so their calls cannot influence their gene floors or global threshold.
    clusters = np.random.default_rng(0).permutation(sorted(e.type.unique()))
    held_out = pd.Series(False, index=e.index)
    for fold in np.array_split(clusters, 5):
        test = e.type.isin(fold)
        _, fold_thr = fit_threshold(e[~test], lig, prone)
        held_out.loc[test] = e.loc[test, "val"] >= e.loc[test, "gene"].map(fold_thr)
    report.append(dict(evaluation="held-out clusters", theta=np.nan, **classifier_scores(e, held_out)))
    pd.DataFrame(report).to_csv(T / "ligand_calibration.csv", index=False)
    print(pd.DataFrame(report).to_string(index=False))
    return M.ge(thr, axis=1), means.set_index("cluster")[["cluster_alias", "class"]], thr


def spatial_evidence(d):
    """Screen each brain and compartment separately, preserving unavailable assays and paired disagreements."""
    proj = d.projected_hit.any()
    meas = d.measured_hit.fillna(False).any()
    if proj and meas:
        status = "Measured and projected"
    elif meas:
        status = "Measured only"
    elif proj:
        status = "Projected only"
    elif (d.call == "ubiquitous").any():
        status = "Unresolved at non-neuronal levels"
    else:
        status = "Below threshold"
    if not d.measured.any():
        status += " (off panel)"
    # Select the strongest observation for display, retaining its method and anatomical origin.
    peaks = d.set_index(["dataset", "compartment"])[["meas_frac", "proj_frac"]].stack().dropna()
    ds, comp, column = peaks.idxmax()
    discordant = d[d.evidence.isin(["measured only", "projected only"])]
    return dict(spatial_evidence=status, spatial_hit=bool(proj or meas), peak_dataset=ds,
                peak_compartment=comp, peak_method="measured" if column == "meas_frac" else "projected",
                peak_fraction=peaks.max(), n_discordant=len(discordant))


def main():
    ligands = pd.read_csv(CONFIG / "neuropeptide_ligands.csv").set_index("gene")
    pairs = pd.read_csv(CONFIG / "ligand_receptor_pairs.csv").fillna("")
    regions = pd.read_csv(CONFIG / "vmh_afferent_regions.csv")
    lig = list(ligands.index)
    expr, meta, thr = cluster_expression(lig)

    # whole-brain composition of each atlas structure (Allen MERFISH brain, released CCF labels)
    xs = pd.read_parquet(D / "merfish_cluster_by_structure.parquet").dropna(subset=["parcellation_structure"])
    alias = meta.reset_index().set_index("cluster_alias")
    xs = xs[xs.cluster_alias.isin(alias.index)]
    xs["cluster"] = alias.cluster.reindex(xs.cluster_alias).values
    xs["neuron"] = ~alias["class"].reindex(xs.cluster_alias).str.contains(NON_NEURONAL).values
    xs = xs[xs.neuron]
    n_region = xs.groupby("parcellation_structure").n_cells.sum()
    E = expr.reindex(xs.cluster).values.astype(float)
    num = pd.DataFrame(E * xs.n_cells.values[:, None], columns=lig).groupby(xs.parcellation_structure.values).sum()
    share = num.div(n_region, axis=0)
    long = share.stack().rename("share_of_neurons").reset_index()
    long.columns = ["structure", "gene", "share_of_neurons"]
    long["n_neurons_in_structure"] = long.structure.map(n_region).astype(int)
    long["n_expressing"] = num.stack().values.round().astype(int)
    total = num.sum(axis=0)
    long["share_of_brainwide_expressing"] = (num.div(total.where(total > 0), axis=1)).stack(future_stack=True).values
    long["is_source"] = (long.share_of_neurons >= AFFERENT_MIN) | (
        (long.share_of_brainwide_expressing >= SPECIALISED_SHARE) & (long.n_expressing >= SPECIALISED_MIN_CELLS))
    div = xs.drop_duplicates("parcellation_structure").set_index("parcellation_structure").parcellation_division
    long.insert(1, "division", long.structure.map(div))
    region_meta = regions.set_index("acronym")
    long["afferent_evidence"] = long.structure.map(region_meta.evidence_type)
    long["afferent_citation"] = long.structure.map(region_meta.citation)
    long["reported_vmh_afferent"] = long.afferent_evidence == "mouse tracing citation"
    long[long.n_expressing > 0].to_csv(T / "ligand_by_region.csv", index=False, float_format="%.4g")

    gb = pd.read_csv(T / "gene_by_compartment.csv")
    gb = gb[gb.compartment.isin(COMPARTMENTS)]
    evidence = {g: spatial_evidence(d) for g, d in gb.groupby("gene")}
    loc = {ds: gb[gb.dataset == ds].pivot(index="gene", columns="compartment", values="proj_frac")
           for ds in DATASETS}
    rows = []
    for g in lig:
        l = loc["allen_m"].loc[g].reindex(COMPARTMENTS)
        lg = long[(long.gene == g) & long.is_source]
        aff = lg[lg.reported_vmh_afferent].sort_values("n_expressing", ascending=False).set_index("structure").share_of_neurons
        top = lg[lg.n_neurons_in_structure >= MIN_REGION_NEURONS].sort_values("n_expressing", ascending=False).head(5)
        top = top.set_index("structure").share_of_neurons
        rows.append(dict(
            gene=g, peptides=ligands.loc[g, "peptides"], family=ligands.loc[g, "family"], source_class=ligands.loc[g, "source_class"],
            **{f"local_{c}": l[c] for c in COMPARTMENTS},
            local_max=l.max(), local_where=l.idxmax() if gb[(gb.gene == g) & (gb.dataset == "allen_m")].projected_hit.any() else "",
            local_max_female_brain=loc["zhuang1_f"].loc[g].reindex(COMPARTMENTS).max(),
            afferent_regions="; ".join(f"{r} {v:.0%}" for r, v in aff.items()),
            n_afferent_regions=len(aff),
            top_regions_brainwide="; ".join(f"{r} {v:.0%}" for r, v in top.items()),
            **evidence[g], local=evidence[g]["spatial_hit"], afferent=bool(len(aff)),
            endocrine=bool((pairs.loc[pairs.ligand_gene == g, "endocrine_source"] != "").any()),
            provisional_regions="; ".join(lg.loc[lg.afferent_evidence.notna() & ~lg.reported_vmh_afferent,
                                                   "structure"].sort_values())))
    src = pd.DataFrame(rows)
    src.to_csv(T / "ligand_sources.csv", index=False, float_format="%.4g")
    S = src.set_index("gene")

    # receptors in VMH against available ligands
    rec = pd.read_csv(CONFIG / "neuropeptide_receptors.csv").set_index("gene")
    out = []
    for r in rec.index:
        rr = gb[(gb.dataset == "allen_m") & (gb.gene == r) & gb.compartment.isin(COMPARTMENTS)].set_index("compartment")
        frac = rr.proj_frac.reindex(COMPARTMENTS)
        annotated = pairs[pairs.receptor_gene == r]
        eligible = annotated[annotated.qualification == ""]
        cog = list(eligible.ligand_gene.unique())
        conditional = "; ".join(f"{v.peptide} ({v.ligand_gene}): {v.qualification}"
                                for v in annotated[annotated.qualification != ""].itertuples())
        local = [g for g in cog if S.loc[g, "local"]]
        aff = [g for g in cog if S.loc[g, "afferent"]]
        endo = [f"{v.peptide} ({v.ligand_gene}, {v.endocrine_source})"
                for v in eligible[eligible.endocrine_source != ""].itertuples()]
        present = evidence[r]["spatial_hit"]
        if not present:
            local, aff, endo = [], [], []
            verdict = evidence[r]["spatial_evidence"]
        elif not cog:
            verdict = "conditional partners only" if conditional else "no ligand gene listed"
        else:
            parts = [n for n, v in (("local", local), ("afferent", aff), ("endocrine", endo)) if v]
            verdict = " + ".join(parts) if parts else "no source identified"
        out.append(dict(receptor=r, family=rec.loc[r, "family"], principal_coupling=rec.loc[r, "principal_coupling"],
                        **{f"vmh_{c}": frac[c] for c in COMPARTMENTS}, vmh_max=frac.max(),
                        vmh_where=frac.idxmax() if rr.projected_hit.any() else "",
                        **evidence[r], conditional_partners=conditional,
                        in_vmh=present, cognate_ligand_genes="; ".join(annotated.ligand_gene.unique()), local_ligands="; ".join(local),
                        afferent_ligands="; ".join(f"{g} ({S.loc[g, 'afferent_regions']})" for g in aff),
                        endocrine_ligands="; ".join(endo), ligand_sources=verdict))
    con = pd.DataFrame(out).sort_values(["in_vmh", "vmh_max"], ascending=False)
    con.to_csv(T / "receptor_ligand_concordance.csv", index=False, float_format="%.4g")
    inv = con[con.in_vmh]
    print(f"ligands: {int(src.local.sum())} made locally, {int(src.afferent.sum())} made in a reported afferent region, "
          f"{int(src.endocrine.sum())} endocrine, of {len(src)}")
    print(f"receptors in VMH: {len(inv)}; ligand sources: {inv.ligand_sources.value_counts().to_dict()}")
    # ligands made locally whose listed receptors are all absent from VMH
    present = set(inv.receptor)
    orphan = [g for g in src[src.local].gene if not (set(pairs[pairs.ligand_gene == g].receptor_gene) & present)]
    print("made locally with no listed receptor detected in VMH:", orphan)


if __name__ == "__main__":
    main()
