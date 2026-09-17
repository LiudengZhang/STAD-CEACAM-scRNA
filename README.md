# STAD-CEACAM-scRNA

> Figure-generation code for:
>
> **CEACAM5/6⁺ Tumor Cells and IL-1β⁺ Myeloid Cells Mark Distinct States of Resistance to Chemo-immunotherapy in Gastric Cancer**
> Posted as a bioRxiv preprint under an earlier title: bioRxiv 2026.03.05.708917 (2026) — [Preprint](https://www.biorxiv.org/content/10.64898/2026.03.05.708917v1) · [DOI](https://doi.org/10.64898/2026.03.05.708917)

The final release checks and their scope are recorded in
[REPRODUCIBILITY_AUDIT.md](REPRODUCIBILITY_AUDIT.md).

## Study at a glance

- **542,121 cells** from **35 patients** with advanced gastric cancer treated with anti–PD-1 + chemotherapy
- **Multi-modal**: single-cell RNA-seq + spatial transcriptomics + immunohistochemistry + bulk RNA-seq
- **Two resistance programs identified**:
  - *Intrinsic* — CEACAM5/6⁺ tumor cells form immune-excluded niches with macrophage recruitment & CD8⁺ T-cell exhaustion
  - *Acquired* — IL-1β⁺ macrophages drive NF-κB activation, PD-L1 upregulation, and EMT

## Methods used

| Domain | Tools |
|---|---|
| scRNA-seq core | Scanpy, AnnData |
| Differential abundance | **Milo** (`pertpy`) |
| Gene regulatory networks | **SCENIC** (cisTarget motifs v10nr) |
| Deconvolution of bulk RNA-seq | **BayesPrism** |
| Reproducibility | Conda environments pinned per analysis |

## Repository layout

```
00_Config/                   # Shared config (paths, colors, gene panels)
01_Raw_Inputs/               # Documented input expectations (h5ad, spatial, bulk)
02_Preparation_for_Panels/   # Placeholders; the prepared intermediates are fetched from Zenodo
03_Final_Panels/             # Panel scripts — one folder per main figure
    ├── 01_Figure_1 ... 05_Figure_5
    ├── 10_Supplementaries   # the pre-submission S1–S6 scripts, kept as the record of the submitted pages
    ├── Supplementary_New/   # Supplementary Figures S1–S10 as shipped in the revision (all redrawn;
    │                        # S1 split into S1 and S2, S2–S9 renumbered S3–S10 on 2026-09-16)
    └── _run_all_panels.sh   # Orchestrator
04_Revision_Analyses/        # Analyses added in revision, one folder per reviewer point
05_Manuscript/04_Tables/     # Supplementary Tables ST1–ST8, read by the revision scripts
upstream/                    # How the deposited intermediates were produced (SCENIC, NMF, NicheNet, BayesPrism, spatial)
environment.yml              # Main conda env (stad_ceacam)
```

`upstream/` records how the deposited intermediate results were produced; it is reference material and is not executed by `./run` or `_run_all_panels.sh`.

## Environment Setup

### Main environment

```bash
conda env create -f environment.yml
conda activate stad_ceacam
```

### Milo environment (for Figure 2F)

The Milo differential abundance analysis uses a dedicated environment with `pertpy` and R dependencies:

```bash
conda create -n pertpy_milo -c conda-forge python=3.11 -y
conda run -n pertpy_milo pip install pertpy scanpy matplotlib seaborn filelock
conda install -n pertpy_milo -c conda-forge -c bioconda rpy2 r-base bioconductor-edger bioconductor-limma r-statmod -y
```

### CellTypist and pyDESeq2 environment

`pydeseq2` requires `numpy>=2`, and the main environment is pinned to numpy
1.23.5 — the version every figure in the paper was produced with. Installing it
alongside would silently upgrade numpy and break `scanpy` and `anndata`, so it
gets its own environment. `_run_all_panels.sh` selects it for
`celltypist_annotation.py`.

```bash
conda create -n stad_numpy2 -c conda-forge python=3.11.14 -y
conda run -n stad_numpy2 python -m pip install numpy==2.4.6 scanpy==1.11.5 \
  anndata==0.12.19 celltypist==1.7.1 pydeseq2==0.5.4
```

The deposited `celltypist_labels.csv` records the supporting CellTypist run.
It is not consumed by Supplementary Figures S1–S10, so the revision driver
skips recomputation by default. Set `STAD_RUN_OPTIONAL=1` to run it; CellTypist
then downloads its `Immune_All_Low.pkl` model if it is not already cached.
The same flag enables the supporting CellPhoneDB interaction analysis, which
requires `02_External/CellPhoneDB/cellphonedb.zip`; its output is also not
consumed by S1–S10. `pin_gene_sets.py` is a preparation utility and is always
skipped because the pinned Hallmark GMT already ships under `00_Reference/`.
`rescore_mast_gsea.py` is retained as an audit/preparation utility but is also
skipped: Figure S10 uses the deposited, adopted `sound13` result tables.
`momac_lineage.py` is likewise skipped by default: the approved lineage-score
tables used by Figure S10A/B are deposited with the code, while recomputation
requires the same feature-selected source matrix used for the analysis. Set
`STAD_RUN_PREPARATION=1` only when that matching source object is available.
The approved S10B dotplot assets are deposited with those tables because the
public feature-complete MoMac object gives different scaled dotplot means.
The driver also skips the internal method/audit modules 14, 16, and 17. Module
15 is retained as the complete pseudobulk reference workflow; run its own
`run_all.sh` rather than executing its Python files independently. The released
`nfkb_comparison.csv` and Supplementary Table S8 are prepared final inputs.
Module 01 is retained as a private clinical-audit reference and is skipped by
the public driver; its de-identified cached outputs are deposited.
`spatial_positive_evidence.py` is retained as an audit/reference script and is
skipped because its internal S11 output is outside the submitted S1–S10 set.

### BayesPrism and ESTIMATE (R)

The BayesPrism pipeline in `upstream/BayesPrism/` is reference only; its
results are in the Zenodo record. The revision workflow also uses the same R
environment for the TCGA ESTIMATE immune score, so both packages are required:

```bash
conda create -n r_bayesprism -c conda-forge r-base r-data.table r-devtools -y
conda run -n r_bayesprism R -e 'if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(c("NMF","scran")); devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")'
conda run -n r_bayesprism R -e "install.packages('estimate', repos='http://r-forge.r-project.org')"
```

## Usage

Run all figure panels:

```bash
bash 03_Final_Panels/_run_all_panels.sh
```

Run a single figure, or the supplementary figures as shipped:

```bash
bash 03_Final_Panels/_run_all_panels.sh 3        # Figure 3 panels only
bash 03_Final_Panels/_run_all_panels.sh revision # revision analyses and Supplementary Figures S1–S10
bash 03_Final_Panels/_run_all_panels.sh supp     # the pre-submission S1–S6 scripts (10_Supplementaries)
# Optional supporting analyses, including CellTypist model inference:
STAD_RUN_OPTIONAL=1 bash 03_Final_Panels/_run_all_panels.sh revision
```

Every panel of Figures 1–5 is drawn at the size it prints at. The revised
pages themselves (171.10 × 229.31 mm, the panels placed into the published
slot geometry recorded in `03_Final_Panels/panel_rects_v2.csv`) are assembled
by a slot assembler that is not part of this release; the `assemble_figure_*.py`
scripts beside the panel folders are the pre-submission assemblers and print
the submitted layouts, not the revised pages. Supplementary Figures S1–S10 are
assembled here in full by `Supplementary_New/assemble_new_supplementaries.py`.

## External Databases

The SCENIC pipeline requires cisTarget motif databases (~95 MB, not bundled). Download into `02_Preparation_for_Panels/SCENIC/databases/`:

```bash
# cisTarget motifs database (Aerts lab)
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  -O 02_Preparation_for_Panels/SCENIC/databases/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

## Data Availability

Processed single-cell objects and the H&E/IHC images are deposited at Zenodo,
[DOI:10.5281/zenodo.18737073](https://doi.org/10.5281/zenodo.18737073), which
always resolves to the current version of the record. The same record carries the
prepared intermediates the panels read (`02_Preparation_for_Panels/`). Set
`STAD_RAW_INPUTS` and `STAD_PREPARED_INPUTS` to point at your copy.

Public datasets are obtained from their own sources: GSE251950 (spatial),
GSE183904 and GSE239676 (single-cell validation), PRJEB25780 (bulk, TIGER) and
TCGA-STAD from the NCI GDC. Supplementary Tables ST1–ST8 accompany the paper and
are also included here, because three revision scripts read the cohort and
signature definitions out of them.

## Citation

```bibtex
@article{chen2026ceacam,
  title   = {CEACAM5/6+ Tumor Cells and IL-1β+ Myeloid Cells Mark Distinct States of Resistance to Chemo-immunotherapy in Gastric Cancer},
  author  = {Chen, Jian and Zhang, Liudeng and Luo, Yikai and Han, Xiaying and Kang, Muxing
             and Chen, Jing and Liu, Wei and Xun, Zhenzhen and Chen, Guofeng and Chen, Ke
             and Xu, Shenbin and Zhang, Chaoyang and Wu, Zhiwei and Wu, Wenxuan
             and Hao, Zhixing and Han, Yaxuan and Lin, Qiaowei and Xu, Yewei
             and Wang, Lie and Liang, Han},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.03.05.708917},
  url     = {https://www.biorxiv.org/content/10.64898/2026.03.05.708917v1}
}
```

## License

MIT
