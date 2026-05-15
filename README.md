# STAD-CEACAM-scRNA

> Figure-generation code for the bioRxiv preprint:
>
> **CEACAM5/6⁺ Tumor Cells and IL-1β⁺ Macrophages Drive Resistance to Chemo-immunotherapy in Gastric Cancer**
> Chen J, **Zhang L**, Luo Y, *et al.* (Liang H, senior author). bioRxiv 2026.03.05.708917 (2026).
> [📄 Preprint](https://www.biorxiv.org/content/10.64898/2026.03.05.708917v1) · [DOI](https://doi.org/10.64898/2026.03.05.708917)

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
| Spatial transcriptomics | (per-figure scripts in `02_Preparation_for_Panels/`) |
| Reproducibility | Conda environments pinned per analysis |

## Repository layout

```
00_Config/                   # Shared config (paths, colors, gene panels)
01_Raw_Inputs/               # Documented input expectations (h5ad, spatial, bulk)
02_Preparation_for_Panels/   # Per-method preprocessing (SCENIC, Milo, BayesPrism, ST)
03_Final_Panels/             # Figure assembly scripts — one folder per main figure
    ├── 01_Figure_1 ... 05_Figure_5
    ├── 10_Supplementaries
    └── _run_all_panels.sh   # Orchestrator
environment.yml              # Main conda env (119/120 scripts)
```

## Environment Setup

### Main environment (119/120 scripts)

```bash
conda env create -f environment.yml
conda activate stad_ceacam
```

### Milo environment (Figure 2F only)

The Milo differential abundance analysis requires a separate environment with `pertpy` and R dependencies:

```bash
conda create -n pertpy_milo -c conda-forge python=3.11 -y
conda run -n pertpy_milo pip install pertpy scanpy matplotlib seaborn filelock
conda install -n pertpy_milo -c conda-forge -c bioconda rpy2 r-base bioconductor-edger bioconductor-limma r-statmod -y
```

### BayesPrism (R, preparation scripts only)

```bash
conda create -n r_bayesprism -c conda-forge r-base r-data.table r-devtools -y
conda run -n r_bayesprism R -e 'if (!require("BiocManager")) install.packages("BiocManager"); BiocManager::install(c("NMF","scran")); devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")'
```

## Usage

Run all figure panels and assemblies:

```bash
bash 03_Final_Panels/_run_all_panels.sh
```

Run a single figure:

```bash
bash 03_Final_Panels/_run_all_panels.sh 3      # Figure 3 only
bash 03_Final_Panels/_run_all_panels.sh supp   # supplementaries only
```

## External Databases

The SCENIC pipeline requires cisTarget motif databases (~95 MB, not bundled). Download into `02_Preparation_for_Panels/SCENIC/database/`:

```bash
# cisTarget motifs database (Aerts lab)
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  -O 02_Preparation_for_Panels/SCENIC/database/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

## Data Availability

Input data (h5ad files, spatial data, external datasets) are available upon request or from the repositories described in the manuscript.

## Citation

```bibtex
@article{chen2026ceacam,
  title   = {CEACAM5/6+ Tumor Cells and IL-1β+ Macrophages Drive Resistance to Chemo-immunotherapy in Gastric Cancer},
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

## Affiliations

- Zhejiang University School of Medicine — Second Affiliated Hospital (clinical, gastrointestinal surgery)
- The University of Texas MD Anderson Cancer Center — Bioinformatics & Computational Biology
- Baylor College of Medicine — Quantitative and Computational Biosciences

## License

MIT
