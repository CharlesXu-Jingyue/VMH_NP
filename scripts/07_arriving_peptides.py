#!/usr/bin/env python3
"""Step 7 (supplemental part 3): which peptides could be released in the VMH, and do they
match the receptors found there?

Three sources are considered for each ligand gene:
  local      VMH neurons make it (projected fraction by compartment, step 6)
  afferent   neurons in a region reported to project to VMH make it. Regions are listed in
             config/vmh_afferent_regions.csv, curated by hand from the tracing literature and
             meant to be edited.
  endocrine  the ligand is made mainly outside the brain (tag in neuropeptide_ligands.csv)

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
from common import D, NON_NEURONAL  # noqa: E402
from paths import CONFIG, RESULTS  # noqa: E402

T = RESULTS / "tables"
AFFERENT_MIN, LOCAL_MIN, MIN_REGION_NEURONS = 0.10, 0.10, 200
SPECIALISED_SHARE, SPECIALISED_MIN_CELLS = 0.10, 20
AMBIENT_PRONE_CPM, AMBIENT_PRONE_THETA = 2.0, 4.0
VMH_COMPARTMENTS = ["VMHDM", "VMHC", "VMHVL", "VMHSh"]


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
    # per-gene floor: what clusters carrying only background reach (ambient transcripts of abundant peptides)
    floor = e[e.call == "background"].groupby("gene").val.quantile(0.99).reindex(lig).fillna(0.0)
    truth = e.call.isin(["broad", "subset"]).values
    fl = floor.reindex(e.gene).values
    best = None
    for theta in np.arange(0.5, 4.01, 0.1):
        pred = e.val.values >= np.maximum(theta, fl)
        tp = (pred & truth).sum()
        f1 = 2 * tp / max(pred.sum() + truth.sum(), 1)
        if best is None or f1 > best[0]:
            best = (f1, theta, tp / max(truth.sum(), 1), tp / max(pred.sum(), 1))
    f1, theta, rec, prec = best
    print(f"cluster threshold: THETA = {theta:.1f} mean log2(CPM+1); against {len(e)} hypothalamic cluster x gene calls: "
          f"F1 {f1:.2f}, recall {rec:.2f}, precision {prec:.2f}")
    thr = np.maximum(theta, floor)
    # ambient-prone genes: transcripts abundant enough to contaminate droplets tissue-wide. Their true
    # sources express at very high levels, so a stricter threshold loses nothing and removes leakage.
    prone = bg.reindex(lig).fillna(0) >= AMBIENT_PRONE_CPM
    thr = thr.where(~prone, np.maximum(thr, AMBIENT_PRONE_THETA))
    print(f"ambient-prone ligand genes held to mean log2(CPM+1) >= {AMBIENT_PRONE_THETA}: {', '.join(prone[prone].index)}")
    return M.ge(thr, axis=1), means.set_index("cluster")[["cluster_alias", "class"]], thr


def main():
    ligands = pd.read_csv(CONFIG / "neuropeptide_ligands.csv").set_index("gene")
    pairs = pd.read_csv(CONFIG / "ligand_receptor_pairs.csv")
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
    long["reported_vmh_afferent"] = long.structure.isin(regions.acronym)
    long[long.n_expressing > 0].to_csv(T / "ligand_by_region.csv", index=False, float_format="%.4g")

    # local synthesis from step 6 (Allen brain weights; the female Zhuang brain alongside)
    gb = pd.read_csv(T / "gene_by_compartment.csv")
    loc = {ds: gb[(gb.dataset == ds) & gb.compartment.isin(VMH_COMPARTMENTS)].pivot(index="gene", columns="compartment", values="proj_frac")
           for ds in ("allen_m", "zhuang1_f")}
    rows = []
    for g in lig:
        l = loc["allen_m"].loc[g].reindex(VMH_COMPARTMENTS)
        lg = long[(long.gene == g) & long.is_source]
        aff = lg[lg.reported_vmh_afferent].sort_values("n_expressing", ascending=False).set_index("structure").share_of_neurons
        top = lg[lg.n_neurons_in_structure >= MIN_REGION_NEURONS].sort_values("n_expressing", ascending=False).head(5)
        top = top.set_index("structure").share_of_neurons
        rows.append(dict(
            gene=g, peptides=ligands.loc[g, "peptides"], family=ligands.loc[g, "family"], source_class=ligands.loc[g, "source_class"],
            **{f"local_{c}": l[c] for c in VMH_COMPARTMENTS},
            local_max=l.max(), local_where=l.idxmax() if l.max() >= LOCAL_MIN else "",
            local_max_female_brain=loc["zhuang1_f"].loc[g].reindex(VMH_COMPARTMENTS).max(),
            afferent_regions="; ".join(f"{r} {v:.0%}" for r, v in aff.items()),
            n_afferent_regions=len(aff),
            top_regions_brainwide="; ".join(f"{r} {v:.0%}" for r, v in top.items()),
            local=bool(l.max() >= LOCAL_MIN), afferent=bool(len(aff)), endocrine=ligands.loc[g, "source_class"] == "endocrine"))
    src = pd.DataFrame(rows)
    src.to_csv(T / "ligand_sources.csv", index=False, float_format="%.4g")
    S = src.set_index("gene")

    # receptors in VMH against available ligands
    rec = pd.read_csv(CONFIG / "neuropeptide_receptors.csv").set_index("gene")
    out = []
    for r in rec.index:
        rr = gb[(gb.dataset == "allen_m") & (gb.gene == r) & gb.compartment.isin(VMH_COMPARTMENTS)].set_index("compartment")
        frac = rr.proj_frac.reindex(VMH_COMPARTMENTS)
        ubiq = bool((rr.call == "ubiquitous").any())
        cog = [g for g in pairs[pairs.receptor_gene == r].ligand_gene.unique() if g in S.index]
        local = [g for g in cog if S.loc[g, "local"]]
        aff = [g for g in cog if S.loc[g, "afferent"]]
        endo = [g for g in cog if S.loc[g, "endocrine"]]
        present = bool(frac.max() >= LOCAL_MIN) or ubiq
        if not present:
            verdict = "receptor not detected in VMH"
        elif not cog:
            verdict = "no ligand gene listed"
        else:
            parts = [n for n, v in (("local", local), ("afferent", aff), ("endocrine", endo)) if v]
            verdict = " + ".join(parts) if parts else "no source identified"
        out.append(dict(receptor=r, family=rec.loc[r, "family"], principal_coupling=rec.loc[r, "principal_coupling"],
                        **{f"vmh_{c}": frac[c] for c in VMH_COMPARTMENTS}, vmh_max=frac.max(),
                        vmh_where=frac.idxmax() if frac.max() >= LOCAL_MIN else ("ubiquitous" if ubiq else ""),
                        in_vmh=present, cognate_ligand_genes="; ".join(cog), local_ligands="; ".join(local),
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
