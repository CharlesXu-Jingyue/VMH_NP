#!/usr/bin/env python3
"""Write the curated neuropeptide tables in config/.

Three tables:
  neuropeptide_ligands.csv     one row per precursor gene
  neuropeptide_receptors.csv   one row per receptor (or receptor component) gene
  ligand_receptor_pairs.csv    one row per (precursor gene, peptide, receptor gene)

source_class of a ligand:
  neuronal     classical neuropeptide precursor made by central neurons
  endocrine    made mainly outside the brain (or pituitary); reaches VMH by blood or CSF
  granin       granin or granin-like precursor; broadly expressed, receptors mostly unknown
  growth       insulin-like growth factors
principal_coupling is the primary transduction class from general pharmacology and should be
checked against the IUPHAR/BPS Guide to Pharmacology before any quantitative use.
"""
import csv
from pathlib import Path

CONFIG = Path(__file__).resolve().parents[1] / "config"

# gene: (family, peptides, source_class, [(peptide, [receptors])], note)
L = {}
def lig(gene, family, source, pairs, note=""):
    L[gene] = dict(family=family, source_class=source, pairs=pairs, note=note)

OP = "Opioid"
lig("Penk", OP, "neuronal", [("enkephalins", ["Oprd1", "Oprm1"])])
lig("Pdyn", OP, "neuronal", [("dynorphins", ["Oprk1"])])
lig("Pomc", "Melanocortin/opioid", "neuronal",
    [("alpha-MSH", ["Mc3r", "Mc4r"]), ("beta-endorphin", ["Oprm1", "Oprd1"]),
     ("ACTH", ["Mc2r"]), ("alpha-MSH (peripheral targets)", ["Mc1r", "Mc5r"])],
    "Arcuate and NTS neurons; also pituitary")
lig("Pnoc", OP, "neuronal", [("nociceptin/orphanin FQ", ["Oprl1"])])
lig("Agrp", "Melanocortin/opioid", "neuronal", [("AgRP (inverse agonist)", ["Mc3r", "Mc4r"])])
lig("Tac1", "Tachykinin", "neuronal", [("substance P", ["Tacr1"]), ("neurokinin A", ["Tacr2", "Tacr1"])])
lig("Tac2", "Tachykinin", "neuronal", [("neurokinin B", ["Tacr3"])])
lig("Tac4", "Tachykinin", "neuronal", [("hemokinin-1", ["Tacr1"])], "Mostly peripheral")
lig("Cck", "CCK/gastrin", "neuronal", [("CCK-8", ["Cckbr", "Cckar"])])
lig("Gast", "CCK/gastrin", "endocrine", [("gastrin", ["Cckbr"])])
lig("Sst", "Somatostatin", "neuronal", [("somatostatin", ["Sstr1", "Sstr2", "Sstr3", "Sstr4", "Sstr5"])])
lig("Cort", "Somatostatin", "neuronal", [("cortistatin", ["Sstr1", "Sstr2", "Sstr3", "Sstr4", "Sstr5"])])
Y = ["Npy1r", "Npy2r", "Npy4r", "Npy5r", "Npy6r"]
lig("Npy", "NPY", "neuronal", [("NPY", ["Npy1r", "Npy2r", "Npy5r", "Npy6r"])])
lig("Pyy", "NPY", "endocrine", [("PYY", ["Npy2r", "Npy1r", "Npy5r"])], "Gut; small brainstem population")
lig("Ppy", "NPY", "endocrine", [("pancreatic polypeptide", ["Npy4r"])])
lig("Vip", "VIP/PACAP/glucagon", "neuronal", [("VIP", ["Vipr1", "Vipr2"])])
lig("Adcyap1", "VIP/PACAP/glucagon", "neuronal", [("PACAP", ["Adcyap1r1", "Vipr1", "Vipr2"])])
lig("Ghrh", "VIP/PACAP/glucagon", "neuronal", [("GHRH", ["Ghrhr"])])
lig("Gcg", "VIP/PACAP/glucagon", "neuronal",
    [("GLP-1", ["Glp1r"]), ("GLP-2", ["Glp2r"]), ("glucagon", ["Gcgr"])],
    "Central source is the NTS; also gut and pancreas")
lig("Sct", "VIP/PACAP/glucagon", "endocrine", [("secretin", ["Sctr"])])
lig("Gip", "VIP/PACAP/glucagon", "endocrine", [("GIP", ["Gipr"])])
lig("Crh", "CRH", "neuronal", [("CRH", ["Crhr1", "Crhr2"])], "CRH has much lower affinity at Crhr2")
lig("Ucn", "CRH", "neuronal", [("urocortin 1", ["Crhr1", "Crhr2"])])
lig("Ucn2", "CRH", "neuronal", [("urocortin 2", ["Crhr2"])])
lig("Ucn3", "CRH", "neuronal", [("urocortin 3", ["Crhr2"])])
lig("Avp", "Vasopressin/oxytocin", "neuronal", [("vasopressin", ["Avpr1a", "Avpr1b", "Avpr2", "Oxtr"])])
lig("Oxt", "Vasopressin/oxytocin", "neuronal", [("oxytocin", ["Oxtr", "Avpr1a"])])
lig("Trh", "Releasing hormone", "neuronal", [("TRH", ["Trhr", "Trhr2"])])
lig("Gnrh1", "Releasing hormone", "neuronal", [("GnRH", ["Gnrhr"])])
lig("Kiss1", "RFamide", "neuronal", [("kisspeptin", ["Kiss1r"])])
lig("Npff", "RFamide", "neuronal", [("NPFF", ["Npffr2", "Npffr1"])])
lig("Npvf", "RFamide", "neuronal", [("RFRP-1/3", ["Npffr1", "Npffr2"])])
lig("Qrfp", "RFamide", "neuronal", [("QRFP", ["Qrfpr"])])
lig("Prlh", "RFamide", "neuronal", [("prolactin-releasing peptide", ["Prlhr"])])
G = ["Galr1", "Galr2", "Galr3"]
lig("Gal", "Galanin", "neuronal", [("galanin", G)])
lig("Galp", "Galanin", "neuronal", [("galanin-like peptide", G)])
lig("Spx", "Galanin", "neuronal", [("spexin", ["Galr2", "Galr3"])])
lig("Nts", "Neurotensin", "neuronal", [("neurotensin / neuromedin N", ["Ntsr1", "Ntsr2"])])
lig("Nmu", "Neuromedin", "neuronal", [("neuromedin U", ["Nmur2", "Nmur1"])])
lig("Nms", "Neuromedin", "neuronal", [("neuromedin S", ["Nmur2", "Nmur1"])])
lig("Nmb", "Bombesin-like", "neuronal", [("neuromedin B", ["Nmbr"])])
lig("Grp", "Bombesin-like", "neuronal", [("gastrin-releasing peptide", ["Grpr"])])
lig("Hcrt", "Orexin", "neuronal", [("orexin A/B", ["Hcrtr1", "Hcrtr2"])])
lig("Pmch", "MCH", "neuronal", [("MCH", ["Mchr1"])])
lig("Cartpt", "CART", "neuronal", [("CART", ["Gpr160"])], "Receptor unresolved; Gpr160 is a proposed candidate")
lig("Calca", "Calcitonin/CGRP", "neuronal", [("CGRP-alpha", ["Calcrl", "Ramp1"]), ("calcitonin", ["Calcr"])])
lig("Calcb", "Calcitonin/CGRP", "neuronal", [("CGRP-beta", ["Calcrl", "Ramp1"])])
lig("Iapp", "Calcitonin/CGRP", "endocrine", [("amylin", ["Calcr", "Ramp1", "Ramp2", "Ramp3"])], "Pancreas; some central expression reported")
lig("Adm", "Calcitonin/CGRP", "neuronal", [("adrenomedullin", ["Calcrl", "Ramp2", "Ramp3"])])
lig("Adm2", "Calcitonin/CGRP", "neuronal", [("adrenomedullin 2", ["Calcrl", "Ramp1", "Ramp2", "Ramp3"])])
lig("Nppa", "Natriuretic", "endocrine", [("ANP", ["Npr1", "Npr3"])])
lig("Nppb", "Natriuretic", "endocrine", [("BNP", ["Npr1", "Npr3"])])
lig("Nppc", "Natriuretic", "neuronal", [("CNP", ["Npr2", "Npr3"])])
lig("Edn1", "Endothelin", "neuronal", [("endothelin-1", ["Ednra", "Ednrb"])])
lig("Edn3", "Endothelin", "neuronal", [("endothelin-3", ["Ednrb"])])
lig("Agt", "Angiotensin", "neuronal", [("angiotensin II", ["Agtr1a", "Agtr1b", "Agtr2"]), ("angiotensin 1-7", ["Mas1"])], "Astrocytes are the main central source")
lig("Kng1", "Kinin", "endocrine", [("bradykinin", ["Bdkrb2", "Bdkrb1"])])
lig("Kng2", "Kinin", "endocrine", [("bradykinin", ["Bdkrb2", "Bdkrb1"])])
lig("Rln1", "Relaxin/insulin-like", "endocrine", [("relaxin", ["Rxfp1", "Rxfp2"])])
lig("Rln3", "Relaxin/insulin-like", "neuronal", [("relaxin-3", ["Rxfp3", "Rxfp4", "Rxfp1"])], "Nucleus incertus")
lig("Insl3", "Relaxin/insulin-like", "endocrine", [("INSL3", ["Rxfp2"])])
lig("Insl5", "Relaxin/insulin-like", "endocrine", [("INSL5", ["Rxfp4"])], "Colon")
lig("Ins1", "Relaxin/insulin-like", "endocrine", [("insulin", ["Insr"])])
lig("Ins2", "Relaxin/insulin-like", "endocrine", [("insulin", ["Insr"])])
lig("Igf1", "Relaxin/insulin-like", "growth", [("IGF-1", ["Igf1r", "Insr"])])
lig("Igf2", "Relaxin/insulin-like", "growth", [("IGF-2", ["Igf1r", "Insr"])])
lig("Lep", "Adipokine/cytokine-type hormone", "endocrine", [("leptin", ["Lepr"])], "Adipose")
lig("Ghrl", "Ghrelin", "endocrine", [("ghrelin", ["Ghsr"])], "Stomach")
lig("Prl", "Adipokine/cytokine-type hormone", "endocrine", [("prolactin", ["Prlr"])], "Pituitary")
lig("Gh", "Adipokine/cytokine-type hormone", "endocrine", [("growth hormone", ["Ghr"])], "Pituitary")
lig("Nps", "Neuropeptide S", "neuronal", [("neuropeptide S", ["Npsr1"])])
lig("Npb", "Neuropeptide B/W", "neuronal", [("neuropeptide B", ["Npbwr1"])])
lig("Npw", "Neuropeptide B/W", "neuronal", [("neuropeptide W", ["Npbwr1"])])
lig("Pth", "Parathyroid", "endocrine", [("PTH", ["Pth1r"])])
lig("Pthlh", "Parathyroid", "neuronal", [("PTHrP", ["Pth1r"])])
lig("Pth2", "Parathyroid", "neuronal", [("TIP39", ["Pth2r"])])
lig("Apln", "Apelin", "neuronal", [("apelin", ["Aplnr"])])
lig("Apela", "Apelin", "neuronal", [("ELABELA", ["Aplnr"])])
lig("Uts2", "Urotensin", "neuronal", [("urotensin II", ["Uts2r"])])
lig("Uts2b", "Urotensin", "neuronal", [("urotensin II-related peptide", ["Uts2r"])])
lig("Prok1", "Prokineticin", "neuronal", [("prokineticin 1", ["Prokr1", "Prokr2"])])
lig("Prok2", "Prokineticin", "neuronal", [("prokineticin 2", ["Prokr2", "Prokr1"])])
lig("Pcsk1n", "Granin-like", "granin", [("BigLEN", ["Gpr171"]), ("PEN", ["Gpr83"])], "proSAAS; broadly neuronal")
lig("Vgf", "Granin-like", "granin", [("TLQP-21", ["C3ar1"])], "Other VGF peptides have no assigned receptor")
lig("Chga", "Granin-like", "granin", [], "Chromogranin A; no assigned neuronal receptor")
lig("Chgb", "Granin-like", "granin", [], "Chromogranin B")
lig("Scg2", "Granin-like", "granin", [], "Secretogranin II (secretoneurin)")
lig("Scg3", "Granin-like", "granin", [], "Secretogranin III")
lig("Nucb2", "Granin-like", "granin", [], "Nesfatin-1; receptor unknown")

# gene: (family, principal_coupling, note)
R = {}
def rec(genes, family, coupling, note=""):
    for g in genes.split():
        R[g] = dict(family=family, principal_coupling=coupling, note=note)
rec("Oprm1 Oprd1 Oprk1 Oprl1", "Opioid", "Gi/o")
rec("Tacr1 Tacr2 Tacr3", "Tachykinin", "Gq")
rec("Cckar Cckbr", "CCK/gastrin", "Gq")
rec("Sstr1 Sstr2 Sstr3 Sstr4 Sstr5", "Somatostatin", "Gi/o")
rec("Npy1r Npy2r Npy4r Npy5r Npy6r", "NPY", "Gi/o")
rec("Vipr1 Vipr2", "VIP/PACAP/glucagon", "Gs")
rec("Adcyap1r1", "VIP/PACAP/glucagon", "Gs (also Gq)")
rec("Gcgr Glp1r Glp2r Gipr Sctr Ghrhr", "VIP/PACAP/glucagon", "Gs")
rec("Crhr1 Crhr2", "CRH", "Gs")
rec("Avpr1a Avpr1b", "Vasopressin/oxytocin", "Gq")
rec("Avpr2", "Vasopressin/oxytocin", "Gs")
rec("Oxtr", "Vasopressin/oxytocin", "Gq (also Gi/o)")
rec("Trhr Trhr2", "Releasing hormone", "Gq")
rec("Gnrhr", "Releasing hormone", "Gq")
rec("Kiss1r", "RFamide", "Gq")
rec("Npffr1 Npffr2", "RFamide", "Gi/o")
rec("Qrfpr", "RFamide", "Gq")
rec("Prlhr", "RFamide", "Gq")
rec("Galr1 Galr3", "Galanin", "Gi/o")
rec("Galr2", "Galanin", "Gq (also Gi/o)")
rec("Ntsr1 Ntsr2", "Neurotensin", "Gq")
rec("Nmur1 Nmur2", "Neuromedin", "Gq")
rec("Nmbr Grpr", "Bombesin-like", "Gq")
rec("Brs3", "Bombesin-like", "Gq", "Orphan; no endogenous peptide assigned")
rec("Hcrtr1 Hcrtr2", "Orexin", "Gq")
rec("Mchr1", "MCH", "Gi/o (also Gq)")
rec("Mc1r Mc2r Mc3r Mc4r Mc5r", "Melanocortin/opioid", "Gs")
rec("Gpr160", "CART", "unresolved", "Proposed CART receptor; unconfirmed")
rec("Calcr Calcrl", "Calcitonin/CGRP", "Gs", "Ligand selectivity set by the RAMP partner")
rec("Ramp1 Ramp2 Ramp3", "Calcitonin/CGRP", "accessory", "Receptor activity-modifying proteins; not receptors alone")
rec("Npr1 Npr2", "Natriuretic", "guanylyl cyclase")
rec("Npr3", "Natriuretic", "clearance (Gi/o)")
rec("Ednra Ednrb", "Endothelin", "Gq")
rec("Agtr1a Agtr1b", "Angiotensin", "Gq")
rec("Agtr2 Mas1", "Angiotensin", "non-canonical")
rec("Bdkrb1 Bdkrb2", "Kinin", "Gq")
rec("Rxfp1 Rxfp2", "Relaxin/insulin-like", "Gs")
rec("Rxfp3 Rxfp4", "Relaxin/insulin-like", "Gi/o")
rec("Insr Igf1r", "Relaxin/insulin-like", "receptor tyrosine kinase")
rec("Lepr Prlr Ghr", "Adipokine/cytokine-type hormone", "JAK-STAT", "Lepr has long (signalling) and short isoforms; gene-level counts cannot separate them")
rec("Ghsr", "Ghrelin", "Gq")
rec("Npsr1", "Neuropeptide S", "Gs/Gq")
rec("Npbwr1", "Neuropeptide B/W", "Gi/o")
rec("Pth1r", "Parathyroid", "Gs/Gq")
rec("Pth2r", "Parathyroid", "Gs")
rec("Aplnr", "Apelin", "Gi/o")
rec("Uts2r", "Urotensin", "Gq")
rec("Prokr1 Prokr2", "Prokineticin", "Gq")
rec("Gpr171", "Granin-like", "Gi/o", "BigLEN receptor")
rec("Gpr83", "Granin-like", "Gq", "PEN receptor (reported)")
rec("C3ar1", "Granin-like", "Gi/o", "Also the complement C3a receptor; mainly microglial")

# Context genes shown alongside (not peptides or receptors)
CONTEXT = {
    "Pcsk1": "prohormone convertase 1/3", "Pcsk2": "prohormone convertase 2", "Cpe": "carboxypeptidase E",
    "Pam": "amidating enzyme", "Scg5": "7B2, PC2 chaperone", "Crhbp": "CRH binding protein",
    "Nr5a1": "VMH marker (SF-1)", "Fezf1": "VMH marker", "Esr1": "VMHvl marker", "Pgr": "VMHvl marker",
    "Slc17a6": "glutamatergic (Vglut2)", "Slc32a1": "GABAergic (Vgat)", "Gad1": "GABAergic", "Gad2": "GABAergic",
    "Snap25": "pan-neuronal", "Ar": "androgen receptor", "Cyp19a1": "aromatase",
    "Dlk1": "VMH-adjacent marker", "Satb2": "VMH subpopulation marker", "Foxp2": "VMH subpopulation marker",
    "Nkx2-1": "ventral hypothalamic lineage", "Sox14": "VMH lineage marker",
}

def main():
    missing = sorted({r for v in L.values() for _, rs in v["pairs"] for r in rs} - set(R))
    assert not missing, f"receptors used in pairs but not defined: {missing}"
    with open(CONFIG / "neuropeptide_ligands.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["gene", "family", "source_class", "peptides", "note"])
        for g, v in L.items():
            w.writerow([g, v["family"], v["source_class"], "; ".join(p for p, _ in v["pairs"]), v["note"]])
    with open(CONFIG / "neuropeptide_receptors.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["gene", "family", "principal_coupling", "ligand_genes", "note"])
        for g, v in R.items():
            ligs = sorted({lg for lg, lv in L.items() for _, rs in lv["pairs"] if g in rs})
            w.writerow([g, v["family"], v["principal_coupling"], "; ".join(ligs), v["note"]])
    with open(CONFIG / "ligand_receptor_pairs.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(["ligand_gene", "peptide", "receptor_gene", "rank", "endocrine_source", "qualification"])
        for g, v in L.items():
            for pep, rs in v["pairs"]:
                for i, r in enumerate(rs, 1):
                    endocrine = "peripheral or pituitary" if v["source_class"] == "endocrine" else ""
                    if g == "Gcg":
                        endocrine = "pancreas" if pep == "glucagon" else "intestine"
                    if g == "Pomc" and pep == "ACTH":
                        endocrine = "pituitary"
                    qualification = ""
                    if len(rs) > 1 and rs[0] in ("Calcrl", "Calcr"):
                        partners = " or ".join(rs[1:]) if r == rs[0] else rs[0]
                        qualification = f"requires same-cell {partners}"
                    elif g == "Cartpt":
                        qualification = "proposed receptor"
                    elif pep == "calcitonin":
                        qualification = "requires calcitonin splice product"
                    w.writerow([g, pep, r, i, endocrine, qualification])
    with open(CONFIG / "context_genes.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["gene", "role"])
        w.writerows(CONTEXT.items())
    orphans = [g for g, v in R.items() if not any(g in rs for lv in L.values() for _, rs in lv["pairs"])]
    print(f"{len(L)} ligand genes, {len(R)} receptor genes, "
          f"{sum(len(rs) for v in L.values() for _, rs in v['pairs'])} pairs; receptors without a listed ligand: {orphans}")

if __name__ == "__main__":
    main()
