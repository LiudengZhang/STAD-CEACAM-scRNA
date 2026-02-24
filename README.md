# STAD-CEACAM-scRNA

Figure generation code for: **CEACAM5/6 as immunotherapy resistance markers in gastric cancer** (single-cell RNA-seq study).

## Environment Setup

```bash
conda env create -f environment.yml
conda activate stad_ceacam
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

*Manuscript in preparation.*

## License

MIT
