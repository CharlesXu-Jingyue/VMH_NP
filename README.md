# VMH_NP: a neuropeptide atlas of the ventromedial hypothalamus

Goal: for neurons whose somata lie in the mouse ventromedial hypothalamic nucleus (VMH),
map (1) which neuropeptide precursors they express and where, (2) which neuropeptide
receptors they express and where, and (3, supplemental) which peptides are likely to be
released into the VMH from local, afferent and endocrine sources, and whether those
agree with the receptors found in (2).

Deliverables: a static web site in `docs/` (an interactive explorer with its documentation)
and a written report, all built from the tables in `results/`.

## Layout

    config/             curated gene tables, afferent-region list, pinned data manifest
    scripts/            numbered pipeline steps
    results/tables/     derived tables (CSV, one parquet)
    results/figures/    static figures (PNG)
    docs/               the web site: explorer, its data and documentation (GitHub Pages serves this folder)
    archive/            the retired hand-annotation route

Downloaded Allen data live outside the repository, in the `data_root` set in
`config/paths.json` (currently `/Users/charlesxu/local/data/VMH_NP`), under `abc_atlas/`.
Set the environment variable `VMH_NP_DATA` to override it for one run.

## Step 0: download the data

The Allen Brain Cell (ABC) Atlas bucket is not reachable from Claude's sandboxes, so this
step has to be run natively. Standard library only, Python 3.8 or later:

    cd ~/local/repo/VMH_NP
    python3 scripts/00_download_data.py --dry-run     # optional: list files and sizes
    python3 scripts/00_download_data.py               # core + Zhuang tiers for all three brains

Downloads resume if interrupted, and every file is checked against the size and MD5 in
`config/data_manifest.csv` (ABC Atlas release manifest 20250531).

| Tier | Size | Contents | Needed for |
|---|---|---|---|
| `core` | 23.8 GB | 10x hypothalamus matrices (v2, v3), whole-brain per-cluster statistics, Allen MERFISH (500 genes) with CCF parcellation, taxonomy, CCF annotation | Parts 1 to 3 |
| `zhuang` | 4.6 GB | Zhuang MERFISH brains 1 and 2 (1,122 genes) | Second measured spatial panel; replication in other animals |
| `imputed` | 50.2 GB | Imputed MERFISH (8,460 genes per cell) | Optional upgrade of the spatial maps for unmeasured genes |

## Pipeline

Run from the repository root. Steps 1 onward read `<data_root>/abc_atlas` and write small
extracts to `<data_root>/derived` and tables to `results/`.

| Step | Script | What it does |
|---|---|---|
| 0 | `00_download_data.py` | Downloads the pinned Allen files |
| 1 | `01_extract.py` | Cuts VMH-centred windows from the three spatial brains; pulls the curated genes for all 262,175 hypothalamus 10x cells; whole-brain cluster means; atlas crops |
| 2 | `02_expression.py` | Expression of every curated gene per subclass, supertype and cluster, by chemistry and sex, with an ambient-RNA background model and calls |
| 3 | `03_reference_atlas.py` | Labels every spatial cell with the unified atlas of Chon et al. 2019 (VMHDM, VMHC, VMHVL, VMH (undivided ends), VMHSh): interpolates the atlas between its plates, fits a per-brain rigid alignment, and writes labels and smooth compartment outlines per section |
| 4 | `04_atlas_figure.py` | Figures of the aligned compartments over the cell types, and of the interpolated atlas against its plates |
| 5 | `05_vmh_type_summary.py` | Readable tables per VMH supertype (`vmh_supertype_calls.csv`, `vmh_type_counts.csv`) |
| 6 | `06_compartments.py` | Parts 1 and 2 in space: every gene by compartment (projected from the cell-type mix, and measured in situ where the gene is on a panel), by anterior-posterior position, and a sex comparison |
| 7 | `07_arriving_peptides.py` | Part 3: which regions make each ligand, restricted to regions that project to the VMH, and the receptor-ligand concordance table |
| 8 | `08_figures.py` | Summary figures for the report |
| 9 | `09_build_explorer.py` | Builds the explorer (`docs/index.html`, `docs/data.js`) from the tables and `explorer_template.html` |
| 10 | `10_atlas_comparison.py` | Compares the unified atlas and the Allen CCF annotation against the VMH cell types, each with its own best shift |
| 11 | `11_build_docs.py` | Builds the text of the explorer's Documentation tab (`docs/doc.js`) from `documentation_fragment.html`, filling every number from the tables |

`common.py` and `paths.py` are shared helpers. Steps 1 and 2 are the slow ones (several
minutes each; step 1 is staged and resumable, `--force` to redo a stage). Step 3 takes about
a minute and a half and step 10 about two minutes; the others take seconds.

    export VMH_NP_DATA=/Users/charlesxu/local/data/VMH_NP     # optional override
    python3 scripts/01_extract.py
    for s in 02_expression 03_reference_atlas 04_atlas_figure 05_vmh_type_summary \
             06_compartments 07_arriving_peptides 08_figures 09_build_explorer \
             10_atlas_comparison 11_build_docs; do
        python3 scripts/$s.py || break
    done

Python packages: numpy, pandas, scipy, pyarrow, h5py, nibabel, tifffile, scikit-image, matplotlib.

The reference atlas (Dryad doi:10.5061/dryad.t1g1jwsxw, file `1_Anatomical_Labels.zip`) goes in
`<data_root>/derived/reference/kim_unified_atlas/`. Its structure names are in
`config/kim_unified_structures.csv`. The earlier hand-annotation route (seed boundaries and a
browser editor) is kept under `archive/hand_annotation/` and is no longer part of the pipeline.

## Method in brief

1. **Which cells are "VMH neurons".** VMH and its compartments are those of the unified
   atlas of Chon et al. 2019 (Franklin-Paxinos labels in Allen CCF space): dorsomedial,
   central, ventrolateral, and a surrounding shell. Two corrections are applied before a
   cell is labelled.
   *Position of the plates.* Plate k of the atlas is CCF coronal plane 10 (k - 1) + 94,
   counted from zero. This is exact: the template images shipped with the atlas
   (`2_AllenCCF_background.zip`) are identical to those planes of the Allen template
   (`python scripts/03_reference_atlas.py --verify-ap`, which needs that archive and
   `average_template_10.nii.gz`).
   *Interpolation between plates.* The atlas is a stack of coronal plates 100 um apart,
   each drawn on its own, and the tissue sections are not parallel to them: within one
   section the VMH spans 50 to 220 um along the anterior-posterior axis. Reading the
   nearest plate stitches two or three drawings together in one section, with seams. The
   VMH labels are therefore interpolated (shape-based interpolation): each plate's
   compartments become signed distance fields, smoothed in the plane with a 15 um
   Gaussian, and joined across plates by a monotone cubic. The smoothed fields approximate the published plate boundaries. The nucleus tapers to a point halfway beyond its first and last
   plates. Those end plates are undivided; the change from undivided to divided is placed
   halfway to the next plate. See `results/figures/atlas_interpolation.png`.
   *Alignment.* Labels read at the released CCF coordinates are unreliable, because in all
   three brains the VMH cell types are registered 0.1 to 0.4 mm anterior of the atlas VMH:
   only 55 to 77% of neurons in the core compartments are then of a VMH subclass. Each
   brain is therefore aligned by a rigid shift (three numbers) that maximises overlap
   between neurons of the two VMH subclasses and the atlas VMH as a whole. The same neurons determine the shift and its reported F1. Independent anatomical landmarks would be needed to validate subdivision borders.
   Three labels are kept per cell: `kim_raw` (nearest plate, released coordinates),
   `kim_aligned` (interpolated, shifted; used by every later step) and
   `kim_aligned_nearest` (nearest plate, same shift). Interpolation changes the label of
   11 to 13% of the neurons in or next to the VMH, nearly all between adjacent compartments
   or at the outer edge. See `results/tables/atlas_alignment.csv` and
   `atlas_supertype_by_compartment.csv`. Type-level results (step 2) do not depend on any
   of this.
   *Atlas choice.* The unified atlas is used for its compartments, which the Allen CCF
   annotation lacks, and not because its outer boundary fits better. Step 10 gives each
   outline its own best shift and finds that the Allen CCF VMH agrees with the VMH
   subclasses at least as well as the unified atlas does in every brain. The unified core
   is centred about 120 um posterior of the Allen VMH in the same CCF space.
2. **Expression calls.** Abundant peptide transcripts leak into every droplet, so nonzero
   detection is not evidence of expression. For each gene and library the background is the
   lowest mean among non-neuronal classes; a type is called expressing when its cells
   exceed that background in at least 10% of cells beyond chance and by a threshold factor
   (details in `02_expression.py`). Calls: broad (>=50% of cells), subset (20 to 50%),
   sparse (10 to 20%), glia_level (detected in most cells but not above the non-neuronal
   level), background.
3. **Space.** Two estimates per compartment, always labelled. *Projected*: the neurons
   of a compartment are replaced by their cell types and each type contributes its 10x
   fraction of expressing cells above background (cluster level where the cluster has at
   least 20 cells, otherwise supertype, otherwise subclass); types whose call is background
   contribute zero. *Measured*: for genes on a panel (Allen brain, male, 500 target genes and 50 blank controls; Zhuang
   brains 1 female and 2 male, 1,122 genes), the fraction of neurons above the 95th
   percentile of non-neuronal cells in the same brain. The two agree with Spearman
   correlations of 0.81 to 0.87. Compartments are DM, C, VL, shell, and the undivided
   anterior and posterior ends; "core" is everything except the shell. A parallel breakdown
   by the nine VMH supertypes is kept throughout, because the Nr5a1 and Fezf1 types are
   developmental lineages and do not respect compartment borders.
4. **Arriving peptides (supplemental).** A cluster is taken to make a ligand when its
   whole-brain mean expression clears a threshold calibrated against the hypothalamic calls
   of step 2 (the final calibration metrics are in ligand_calibration.csv), with a per-gene background floor and a stricter threshold for
   ambient-prone genes. A region is a source when at least 10% of its neurons belong to
   such clusters, or when it holds at least 10% of all neurons of those clusters in the
   brain. Regional matches require a mouse tracing citation in the curated afferent list.
   Distribution-only and uncited entries are provisional. Receptor complexes and proposed
   partners are qualified separately. Endocrine evidence follows the peptide product.
   Regional matches are candidates because ligand-producing cells may differ from the projecting cells.

For local evidence, each of DM, C, VL and shell is screened independently in each brain
(`SCREENED` in `scripts/common.py`). A gene is listed when its projected detection reaches 10%
of a compartment's neurons, or when its measured detection does and also exceeds the
blank-probe control: the 95th percentile, across the 50 blank probes of the Allen panel, of
the fraction of that compartment's neurons a blank probe passes under the same rule (16 to 20%;
`tables/blank_probe_floor.csv`). The Zhuang matrices carry no blank probes, so a measured value
from those brains confirms a projection but does not list a gene by itself. The undivided ends
are reported in every table but not screened, because in two brains under a tenth of their
neurons are VMH cell types.
Paired comparisons identify measured-only and projected-only results for that brain and compartment.
Off-panel measurements are missing. Non-neuronal-level detection alone remains unresolved.
The pooled core is a descriptive aggregate. The documentation also reports within-gene spatial
rank agreement alongside the pooled correlation.

## Results

| File | Contents |
|---|---|
| `tables/expr_by_type.parquet` | Every curated gene by subclass, supertype and cluster, chemistry and sex, with background, excess fraction, enrichment and call |
| `tables/gene_background.csv` | Ambient background per gene and chemistry |
| `tables/vmh_supertype_calls.csv`, `vmh_type_counts.csv` | Calls and cell counts for the nine VMH supertypes |
| `tables/atlas_alignment.csv`, `atlas_supertype_by_compartment.csv` | Per-brain shift, agreement with the cell types (raw, aligned, and aligned with the nearest plate), share of neurons relabelled by interpolation, and where each supertype sits |
| `tables/compartment_composition.csv` | Cell-type mix of each compartment in each brain |
| `tables/gene_by_compartment.csv`, `gene_by_ap.csv` | Projected and measured expression by compartment and by anterior-posterior bin, with the paired evidence call and whether the row counts towards the screen |
| `tables/blank_probe_floor.csv` | What the blank probes of the Allen panel reach in each compartment under the measured rule |
| `tables/sex_differences.csv` | Female against male donors per supertype and gene, pooling eligible libraries |
| `tables/ligand_by_region.csv`, `ligand_sources.csv`, `receptor_ligand_concordance.csv` | Part 3 |
| `figures/*.png` | The interpolated atlas against its plates, alignment per brain, supertype by compartment, ligands and receptors by compartment and by supertype, projected against measured |
| `tables/atlas_comparison.csv`, `atlas_overlap.csv` | The unified atlas against the Allen CCF annotation: agreement with the VMH cell types, and how the two outlines overlap |

## The web site

`docs/` is a static site: one page, two script files that carry its data and text, and six
figures. It needs no server-side code and nothing is built on the host:

| File | Contents |
|---|---|
| `docs/index.html` | The explorer. Five tabs: a map of neurons in each brain coloured by one gene (measured in situ or projected) with the compartment outlines, which can be saved as PNG, SVG or PDF; a table of all genes by compartment or supertype; the ligand sources of Part 3; notes, where a reader marks genes, writes remarks and exports them (they stay in the reader's browser); and the documentation |
| `docs/data.js` | Its data (about 4 MB), loaded as a script so the page also works when opened straight from the folder |
| `docs/doc.js`, `docs/figures/` | The Documentation tab: method, data, how to read the explorer, main results and limits. Every number is filled in from `results/tables` by step 11 |

To look at it locally, open `docs/index.html` in a browser. To host it on GitHub Pages: in the
repository settings, under Pages, choose "Deploy from a branch", branch `main`, folder `/docs`.
The site then appears at `https://<user>.github.io/VMH_NP/`. The empty `docs/.nojekyll` tells
Pages to serve the files as they are.

## Interpretation limits

Transcript presence is evidence for capacity, not for release or for functional receptor
at the soma. Precursors need processing (Pcsk1, Pcsk2, Cpe are reported alongside), one
precursor can yield several peptides, receptor protein can be trafficked to terminals
outside the VMH, and sparse receptor transcripts are under-detected by droplet sequencing.

Known gaps: Part 3 has no projection-specific evidence (the Allen connectivity service was
not reachable when this was built); the afferent list and the structure-name table are
unchecked third-party or hand-curated inputs; the sex comparison has few independent donors per
sex and supertype and is underpowered; the imputed MERFISH tier is unused.

## Gene tables

`scripts/build_gene_tables.py` writes `config/neuropeptide_ligands.csv` (89 precursor
genes), `config/neuropeptide_receptors.csv` (98 receptor or receptor-component genes) and
`config/ligand_receptor_pairs.csv` (159 pairings). They are hand-curated; edit the script
and re-run it to add or remove genes. The `principal_coupling` column should be checked
against the IUPHAR/BPS Guide to Pharmacology before quantitative use.

## Data terms

The Allen whole-mouse-brain datasets are shared under the CC BY-NC 4.0 licence. The tables in
`results/` and the data behind the site in `docs/` are derived from them and are offered on the
same terms: non-commercial use, with attribution to the Allen Institute, the Zhuang laboratory,
and Chon et al. 2019 for the unified atlas. This repository does not yet carry a licence for its
own code.

## Checking a build

Run `python -m unittest discover -s tests -v` after building the explorer. The checks compare
the published evidence against the compartment table and the map’s background flags against
the full-precision spatial extracts under the configured data root. A local fixture exercises
checksum failures without downloading any atlas data.
