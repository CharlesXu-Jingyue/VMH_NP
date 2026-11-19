#!/usr/bin/env python3
"""Step 2: per-cell-type expression of curated genes in the 10x hypothalamus data,
with an explicit ambient-RNA background model.

Why a background model: abundant hypothalamic peptide transcripts (Pmch, Oxt, Avp, Agrp,
Pomc, ...) leak into every droplet, so "fraction of cells with a nonzero count" is far
above zero in cell types that do not make the peptide. Non-neuronal cells in the same
library serve as the probe: for gene g and library l the background is

    a[g, l] = min over non-neuronal classes k of mean CPM of g in class k of library l

(the minimum discards genes genuinely expressed by one glial class). A cell i with total
counts N_i then detects g from background alone with probability 1 - exp(-N_i a / 1e6).
For a cell type c:

    enrichment[c, g] = mean CPM in c / mean background of its cells
    excess[c, g]     = (f_obs - f_bg) / (1 - f_bg)      fraction expressing beyond background

A type is called expressed when excess >= EXCESS_MIN and enrichment >= ENRICH_MIN, or when
most cells exceed background (excess >= 0.5) with enrichment >= ENRICH_SOFT. Genes that
every non-neuronal class also expresses get a high background; types that detect such a
gene in most cells without exceeding the non-neuronal level are labelled "glia_level"
instead of "background" (see gene_background.csv, bg_fraction_of_bulk).

Outputs (results/tables):
    expr_by_type.parquet   long table: level, type, chemistry, sex, gene, statistics, call
    gene_background.csv    per gene and chemistry: bulk CPM, background CPM
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import RESULTS, abc_file, data_root  # noqa: E402

D = data_root() / "derived"
OUT = RESULTS / "tables"
NON_NEURONAL = "Astro|Oligo|OPC|Vascular|Immune|OEC"
ENRICH_MIN, ENRICH_SOFT, EXCESS_MIN, MIN_CELLS = 3.0, 1.5, 0.10, 20
SPARSE_SPECIFICITY = 2.0


def load():
    cells = pd.read_parquet(D / "tenx_cells.parquet").reset_index(drop=True)
    X = sp.load_npz(D / "tenx_expr.npz").tocsr()
    genes = pd.read_csv(D / "tenx_genes.csv").gene.tolist()
    C = X.copy()
    C.data = np.exp2(C.data) - 1.0                       # log2(CPM+1) -> CPM
    # Infer one UMI from the full transcriptome before selecting the curated genes.
    for matrix, group in cells.groupby("feature_matrix_label"):
        directory = matrix.rsplit("-HY", 1)[0]
        with h5py.File(abc_file(directory, matrix + "/log2"), "r") as source:
            labels = pd.Index(source["obs/cell_label"].asstr()[:])
            positions = labels.get_indexer(group.cell_label)
            if (positions < 0).any():
                raise ValueError(f"Unmatched cell IDs in {matrix}")
            ptr = source["X/indptr"][:]
            depth = np.empty(len(labels))
            for start in range(0, len(labels), 2000):
                block = ptr[start:start + 2001]
                if np.any(np.diff(block) == 0):
                    raise ValueError(f"Empty transcriptome in {matrix}")
                values = source["X/data"][block[0]:block[-1]]
                minimum = np.minimum.reduceat(values, block[:-1] - block[0])
                depth[start:start + len(minimum)] = np.rint(1e6 / np.expm1(minimum.astype(float) * np.log(2)))
            cells.loc[group.index, "n_counts"] = depth[positions]
    cells["neuron"] = ~cells["class"].str.contains(NON_NEURONAL)
    return cells, C, genes


def background(cells, C, genes):
    """a[g, library] in CPM, and a per-gene summary."""
    n_genes = len(genes)
    lib_codes, libs = pd.factorize(cells.library_label)
    A = np.zeros((len(libs), n_genes))
    nn = cells[~cells.neuron]
    glob = {}
    for chem, d in nn.groupby("library_method"):
        per_class = [np.asarray(C[group.index].mean(0)).ravel()
                     for _, group in d.groupby("class") if len(group) >= 200]
        glob[chem] = np.min(per_class, axis=0)
    chem_of_lib = cells.groupby("library_label").library_method.first()
    for li, lib in enumerate(libs):
        d = nn[nn.library_label == lib]
        per_class = [np.asarray(C[g.index.values].mean(0)).ravel() for _, g in d.groupby("class") if len(g) >= 50]
        # with fewer than two usable classes fall back to the chemistry-wide estimate
        A[li] = np.min(per_class, axis=0) if len(per_class) >= 2 else glob[chem_of_lib[lib]]
    rows = []
    for chem in cells.library_method.unique():
        m = (cells.library_method == chem).values
        bulk = np.asarray(C[m].mean(0)).ravel()
        w = pd.Series(lib_codes[m]).value_counts(normalize=True)
        bg = (A[w.index.values] * w.values[:, None]).sum(0)
        for j, g in enumerate(genes):
            rows.append(dict(gene=g, chemistry=chem, bulk_cpm=bulk[j], background_cpm=bg[j],
                             bg_fraction_of_bulk=bg[j] / bulk[j] if bulk[j] > 0 else np.nan))
    return A, lib_codes, pd.DataFrame(rows)


def summarise(cells, C, genes, A, lib_codes):
    Cc = C.tocsr()
    Xlog = Cc.copy()
    Xlog.data = np.log2(Xlog.data + 1.0)
    N = cells.n_counts.values
    out = []
    groups = []
    for level in ("subclass", "supertype", "cluster"):
        for chem in ("10Xv3", "10Xv2"):
            for sex in ("all", "F", "M"):
                m = (cells.library_method == chem)
                if sex != "all":
                    m &= cells.donor_sex == sex
                subset = cells[m]
                for name, idx in subset.groupby(level).indices.items():
                    groups.append((level, name, chem, sex, subset.index.values[idx]))
    for level, name, chem, sex, idx in groups:
        if len(idx) < 5:
            continue
        sub = Cc[idx]
        mean_cpm = np.asarray(sub.mean(0)).ravel()
        frac = np.asarray((sub > 0).mean(0)).ravel()
        mean_log = np.asarray(Xlog[idx].mean(0)).ravel()
        a = A[lib_codes[idx]]                                   # cells x genes background CPM
        bg_cpm = a.mean(0)
        f_bg = (1.0 - np.exp(-(N[idx][:, None] * a) / 1e6)).mean(0)
        with np.errstate(divide="ignore", invalid="ignore"):
            enrich = np.where(bg_cpm > 0, mean_cpm / bg_cpm, np.where(mean_cpm > 0, np.inf, 1.0))
            excess = np.clip((frac - f_bg) / (1.0 - f_bg), 0, 1)
        out.append(pd.DataFrame(dict(level=level, type=name, chemistry=chem, sex=sex, gene=genes, n_cells=len(idx),
                                     frac=frac, mean_log2cpm=mean_log, mean_cpm=mean_cpm, background_cpm=bg_cpm,
                                     frac_background=f_bg, excess_frac=excess, enrichment=enrich)))
    df = pd.concat(out, ignore_index=True)
    df = add_neuronal_median(df)
    df["call"] = assign_calls(df)
    for c in ("frac", "mean_log2cpm", "mean_cpm", "background_cpm", "frac_background", "excess_frac", "enrichment",
              "frac_median_neuronal"):
        df[c] = df[c].astype(np.float32)
    return df


def add_neuronal_median(df):
    """Median detection fraction of each gene across neuronal types (same level, chemistry, sex).

    Non-neuronal types carry the suffix "NN" in the Allen taxonomy."""
    name = df.type.astype(str)
    neuronal = ~name.str.contains(r" NN(?:_|$)") & (df.n_cells >= MIN_CELLS)
    med = (df[neuronal].groupby(["level", "chemistry", "sex", "gene"], observed=True).frac.median()
           .rename("frac_median_neuronal").reset_index())
    for c in ("level", "chemistry", "sex", "gene"):
        med[c] = med[c].astype(str)
        df[c] = df[c].astype(str)
    return df.drop(columns="frac_median_neuronal", errors="ignore").merge(med, on=["level", "chemistry", "sex", "gene"], how="left")


def assign_calls(df):
    """broad / subset / sparse: expressed above background in >=50%, 20-50%, 10-20% of cells.
    glia_level: detected in most cells but not above the non-neuronal level, so the test cannot
    separate real expression from background (for example Insr). background: neither, or a
    sparse signal that is no higher than in the typical neuronal type."""
    strong = (df.excess_frac >= EXCESS_MIN) & (df.enrichment >= ENRICH_MIN)
    clear = (df.excess_frac >= 0.5) & (df.enrichment >= ENRICH_SOFT)     # most cells beyond background
    ok = (strong | clear) & (df.n_cells >= MIN_CELLS)
    # a sparse call must also be type-specific: at least twice the median detection across neuronal
    # types. Flat low-level detection in every neuron type (for example Insl5) is not informative.
    specific = df.frac >= SPARSE_SPECIFICITY * df.frac_median_neuronal.fillna(0)
    sparse = ok & (df.excess_frac < 0.2) & specific
    return np.select(
        [df.n_cells < MIN_CELLS, ok & (df.excess_frac >= 0.5), ok & (df.excess_frac >= 0.2), sparse, df.frac >= 0.5],
        ["too_few_cells", "broad", "subset", "sparse", "glia_level"], default="background")


def save(df):
    """Groups with fewer than MIN_CELLS cells carry no usable estimate and are not stored."""
    df = df[df.n_cells >= MIN_CELLS].copy()
    for c in ("level", "type", "chemistry", "sex", "gene", "call"):
        df[c] = df[c].astype("category")
    df["n_cells"] = df.n_cells.astype(np.int32)
    df.to_parquet(OUT / "expr_by_type.parquet", compression="zstd", index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cells, C, genes = load()
    A, lib_codes, gb = background(cells, C, genes)
    gb.to_csv(OUT / "gene_background.csv", index=False, float_format="%.4g")
    df = summarise(cells, C, genes, A, lib_codes)
    save(df)
    print(f"{len(cells)} cells; {df.type.nunique()} types; {len(df)} rows; calls:\n{df.call.value_counts().to_string()}")


if __name__ == "__main__":
    main()
