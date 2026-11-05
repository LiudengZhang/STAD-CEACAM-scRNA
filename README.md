# STAD-CEACAM-scRNA

Figure-generation code for **CEACAM5/6⁺ Tumor Cells and IL-1β⁺ Myeloid Cells Mark Distinct States of Resistance to Chemo-immunotherapy in Gastric Cancer**. Posted as a preprint under an earlier title: [bioRxiv 2026.03.05.708917](https://doi.org/10.64898/2026.03.05.708917).

This repository draws the panels of Figures 1–5 and Supplementary Figures S1–S10 and runs the revision analyses. It holds code and small tables; the data are deposited separately.

## Install

```bash
bash install.sh
conda activate stad_ceacam
```

Three analyses use their own environments (`pertpy_milo`, `stad_numpy2`, `r_bayesprism`). Their recipes are in the [setup notes](SETUP_NOTES.md).

## Data

The single-cell objects, spatial data, IHC scans and prepared intermediates are a restricted Zenodo record, [10.5281/zenodo.18737073](https://doi.org/10.5281/zenodo.18737073), available upon reasonable request. Unpack the archives into one directory, then:

```bash
export STAD_RAW_INPUTS=/path/to/data
export STAD_PREPARED_INPUTS=/path/to/data/02_Preparation_for_Panels
```

Public datasets come from their own sources: GSE251950, GSE183904, GSE239676, PRJEB25780 (TIGER) and TCGA-STAD (NCI GDC). The [setup notes](SETUP_NOTES.md) cover two single-file inputs and the SCENIC motif database.

## Run

```bash
bash 03_Final_Panels/_run_all_panels.sh revision   # run first: analyses, then S1–S10
bash 03_Final_Panels/_run_all_panels.sh            # Figures 1–5, then S1–S10
```

Targets `1`–`5` and `supp` run one figure at a time. Runs produce panels, not assembled pages. Most panels rebuild their `data/` table from the deposited data on first run. `04_Revision_Analyses/DISPATCH.csv` lists every step and why any is skipped.

## Layout

- `00_Config/`, `00_Reference/`: shared configuration and gene sets
- `03_Final_Panels/`: one directory per panel, and the driver
- `04_Revision_Analyses/`: one directory per reviewer point
- `upstream/`: how the deposited intermediates were made; reference only

Directory names are not always the printed panel letters; `03_Final_Panels/PROVENANCE.csv` gives the mapping.

## Citation

Chen J, Zhang L, *et al.* bioRxiv (2026). doi:[10.64898/2026.03.05.708917](https://doi.org/10.64898/2026.03.05.708917)

## License

MIT
