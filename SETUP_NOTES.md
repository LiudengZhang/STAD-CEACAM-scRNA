# Setup notes

The detail behind the README.

## Install

```bash
bash install.sh
conda activate stad_ceacam
```

`install.sh` is two steps and has to be: `conda env create` first, then one
`pip install --no-deps` for the eight PyPI packages. They cannot live in
`environment.yml`, because conda hands a `pip:` block to pip as a *requirements
file*, and `--no-deps` is not an option a requirements file may carry — pip
answers `ERROR: Invalid requirement: --no-deps`. The flag matters: cnsplots and
cellphonedb both declare numpy, scanpy and more with no upper bound, and letting
pip resolve them moves the versions every figure in the paper was produced at.
`install.sh` reads the five pins back out of the environment afterwards and
exits non-zero if any of them moved.

Three analyses need an environment of their own. The driver selects them by
name; create the ones you intend to run.

```bash
# Figure 2's Milo differential-abundance panel
conda create -n pertpy_milo -c conda-forge python=3.11 -y
conda run -n pertpy_milo pip install pertpy scanpy matplotlib seaborn filelock
conda install -n pertpy_milo -c conda-forge -c bioconda rpy2 r-base \
  bioconductor-edger bioconductor-limma r-statmod -y

# CellTypist and pyDESeq2. pydeseq2 needs numpy>=2 and the main environment is
# pinned to 1.26.4, the version every figure in the paper was produced with;
# installing it alongside would upgrade numpy under scanpy and anndata.
conda create -n stad_numpy2 -c conda-forge python=3.11.14 -y
conda run -n stad_numpy2 python -m pip install numpy==2.4.6 scanpy==1.11.5 \
  anndata==0.12.19 celltypist==1.7.1 pydeseq2==0.5.4

# The TCGA ESTIMATE immune score (R). BayesPrism is installed in the same
# environment; the BayesPrism pipeline itself is in upstream/ as a record.
conda create -n r_bayesprism -c conda-forge r-base r-data.table r-devtools -y
conda run -n r_bayesprism R -e 'if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(c("NMF","scran")); devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")'
conda run -n r_bayesprism R -e "install.packages('estimate', repos='http://r-forge.r-project.org')"
```

## Where the data is

A clone of this repository does not contain the single-cell objects, the
spatial data, the IHC scans or the prepared intermediate results. They are a
restricted Zenodo record under
[DOI:10.5281/zenodo.18737073](https://doi.org/10.5281/zenodo.18737073), which
resolves to the current version; access is available upon reasonable request.
Unpack the three archives so that one directory holds `01_H5AD/`,
`02_IHC_Raw/` and `02_Preparation_for_Panels/`, then point the code at it:

```bash
export STAD_RAW_INPUTS=/path/to/data
export STAD_PREPARED_INPUTS=/path/to/data/02_Preparation_for_Panels
```

`01_Raw_Inputs/` and `02_Preparation_for_Panels/` in the clone hold the expected
layout and a few small tables, so a path that is wrong fails where the file was
going to be read rather than deeper in.

Five public datasets are fetched from their own sources rather than from the
record: GSE251950 (spatial), GSE183904 and GSE239676 (single-cell validation),
PRJEB25780 (bulk, TIGER), and TCGA-STAD from the NCI GDC.

Two of them are read as **single files at fixed paths**, so they are named here
rather than left to the dataset accession. Neither is in the Zenodo archives.
The TCGA one is an input and nothing here writes it; the TIGER one **is written
by this repository**, by `upstream/TIGER_PRJEB25780/`, and that is how you get
it:

| file | placed at | read by |
|---|---|---|
| TIGER metadata for PRJEB25780, with the ESTIMATE stromal/immune/purity columns computed here — **regenerate it with `upstream/TIGER_PRJEB25780/`** | `$STAD_RAW_INPUTS/02_External/Bulk/PRJEB25780_Tiger/metadata_with_purity.tsv` | `00_Config/paths.py:TIGER_META`, and through it printed Figure 2 L (`02_Figure_2/02_N`, `02_O`), printed S3 B, the two-sided sweep and the MP external validation |
| TCGA ABSOLUTE purity calls (`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`, GDC PanCanAtlas) | `$STAD_RAW_INPUTS/public_data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | `04_Revision_Analyses/05_R1.6_Spatial_Confounders/scripts/tcga_immune_exclusion.py` |

The TIGER table is the PRJEB25780 sample metadata from
<http://tiger.canceromics.org/> with four columns added — `StromalScore`,
`ImmuneScore`, `ESTIMATEScore` and `TumorPurity` — alongside `sample_id` and
the clinical fields. **Those four are computed here**, by
`upstream/TIGER_PRJEB25780/scripts/02_run_estimate.R`, from TIGER's expression
matrix through the `estimate` package. This paragraph said the metadata was
"joined to that record's ESTIMATE scores" until 2026-10-05, which reads as
though TIGER published them; it did not.

TIGER carries no redistribution licence, so neither the metadata table nor the
reformatted expression matrix is in the Zenodo archives. Both are regenerated
exactly from TIGER's own two downloads — reproduced byte for byte on
2026-10-05, 33 seconds for the pair — and
`upstream/TIGER_PRJEB25780/PROVENANCE.md` carries the md5 of each download and
of each output, so you can tell whether you hold what the paper used.

The TCGA table is the
PanCanAtlas ABSOLUTE purity/ploidy master call set, distributed by the NCI GDC.
Without them those panels and analyses stop at the read, by design: the paths
above are where they are looked for, so a missing file fails where it is needed
rather than deeper in.

SCENIC additionally needs the cisTarget motif database (~95 MB):

```bash
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  -O 02_Preparation_for_Panels/SCENIC/databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

Two sets of result tables ship **inside** this repository rather than with the
data, because the code that would rebuild them must not be re-run: the live
GSEA recompute under `04_Revision_Analyses/12_R1.8_DEG_Recompute/outputs/gsea` and
the thirteen-type sound run under
`04_Revision_Analyses/13_R1.8_Neutrophil_Rebuilt_Recompute/outputs/`.

Panel `data/` tables are not in this repository, except Figure 2F's. A panel
script computes its table from the deposited data when the table is absent.

**Supplementary Tables ST1–ST6 are published with the article**, so they are
not in this repository and neither is the script that assembled them: ST2–ST5
were carried over from the submission rather than computed, so `build_tables.py`
was a carrier for files the journal prints. One of them ships, as a declared
input: `05_Manuscript/04_Tables/ST3_gene_signatures.csv`, which
`04_Revision_Analyses/03_R1.4_MP_Direction_PrePost` reads for the metaprogram
signatures and exits without.

## How to run

Run the revision analyses **first**: the drivers for S8–S10 read tables they
write, and nothing else writes them.

```bash
bash 03_Final_Panels/_run_all_panels.sh revision   # analyses, then S1–S10
bash 03_Final_Panels/_run_all_panels.sh            # the panels of Figures 1–5, then S1–S10
```

Or one target at a time:

```bash
bash 03_Final_Panels/_run_all_panels.sh 3      # Figure 3's panels
bash 03_Final_Panels/_run_all_panels.sh supp   # Supplementary Figures S1–S10
```

**This repository supplies the panels; the journal sets the page** (the
author's ruling of 2026-10-02). Every panel of Figures 1–5 and S1–S10 is drawn
at the millimetre size it prints at, into the published slot geometry —
`03_Final_Panels/panel_rects_v2.csv`, `slot_subrects_v2.csv` and
`letter_spec_v2.csv`, one rect per printed panel,
measured off the published pages and shipped here as data. The five scripts that
place those panels onto a page and annotate it are **not** here:
`assemble_slotted.py`, `Supplementary_New/assemble_new_supplementaries.py`,
`patch_figure_annotations.py`, `build_shipped_figures.py` and
`00_Config/shared/svg_assembler.py`. So there is no `pages` target — the driver
accepts `1`–`5`, `supp` and `revision` and nothing else. A run produces every
panel at the path `03_Final_Panels/PROVENANCE.csv` declares for it, which is
also what the reproduction harness checks.

### What the driver does not run, and why

`04_Revision_Analyses/DISPATCH.csv` is the driver's manifest: one row per unit
in dependency order, with its interpreter, its conda environment, the variable
that enables an optional row, and — for every unit that does not run — the
reason. The driver runs nothing the table does not name. Two kinds of row do
not run by default:

- **steps whose results ship here and must not be recomputed**, because a
  recompute would overwrite a preserved result the figures are drawn from. Five
  rows are `status=skip`, each naming its reason in the table — four because the
  output ships here, one because its library path is absolute on the author's
  machine. The driver prints a real `SKIPPED …` line for each, so a skip is
  visible in the run log and not only recorded in a CSV.
- **one row gated on a variable**: `STAD_RUN_PREPARATION=1` enables the single
  `condition` row, the MoMac lineage recompute, which is reproducible given the
  matching source object. It is the only such variable any shipped file reads.

## Layout

```
00_Config/                    shared configuration: paths, colours, gene panels, the page grid
00_Reference/                 gene sets and the pinned Hallmark GMT
01_Raw_Inputs/                the expected input layout (the data are on Zenodo)
02_Preparation_for_Panels/    the expected prepared-intermediate layout, and small tables
03_Final_Panels/
    01_Figure_1 … 05_Figure_5 one directory per figure, one per panel inside it
    Supplementary_New/        S1–S10, one directory per panel
    _run_all_panels.sh        the driver
    PROVENANCE.csv            which printed panel each directory draws
04_Revision_Analyses/         one directory per reviewer point
    DISPATCH.csv              the driver's manifest (see above)
05_Manuscript/04_Tables/      ST3_gene_signatures.csv only, as a declared input
upstream/                     how the deposited intermediates were produced; reference only,
                              never executed by the driver
environment.yml               the conda half of the stad_ceacam environment
install.sh                    creates it, then the one pip step it needs
```

A panel directory's name is not the printed panel letter. Twelve of Figure 5's
sixteen directories and four of Figure 3's point at a different letter than
their name suggests. `03_Final_Panels/PROVENANCE.csv` is the authority; read a
panel's row there rather than inferring it from the path.
