"""
Refit the purity-adjusted TCGA-STAD model of tcga_immune_exclusion.py with
histological type as a covariate, and within intestinal and diffuse tumours
separately.

Histology from the GDC primary diagnosis: intestinal (intestinal-type, tubular,
papillary), diffuse (diffuse-type, signet-ring cell), mucinous (its own group,
since its Lauren type is not set by the diagnosis), NOS, and other.

Reads outputs/tcga_estimate_scores.csv (from tcga_immune_exclusion.py), the
BayesPrism epithelial expression, ABSOLUTE purity and the GDC clinical table.
Writes outputs/tcga_histology_sensitivity.csv.
Run: python tcga_histology_sensitivity.py
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "00_Config"))
from paths import (PREPARATION, RAW_INPUTS, REVIEWER_MATERIALS,  # noqa: E402
                   TCGA_BAYESPRISM_EPI)

OUT = Path(__file__).resolve().parents[1] / "outputs"
ABSOLUTE_PURITY = (REVIEWER_MATERIALS / "public_data"
                   / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
CLINICAL = (RAW_INPUTS / "02_External" / "Bulk" / "TCGA_STAD" / "02_Raw_Data"
            / "Clinical" / "TCGA_STAD_clinical_data.tsv")

HISTOLOGY = {
    "Adenocarcinoma, intestinal type": "intestinal",
    "Tubular adenocarcinoma": "intestinal",
    "Papillary adenocarcinoma, NOS": "intestinal",
    "Mucinous adenocarcinoma": "mucinous",
    "Carcinoma, diffuse type": "diffuse",
    "Signet ring cell carcinoma": "diffuse",
    "Adenocarcinoma, NOS": "NOS",
}

MODELS = {
    "purity": "{y} ~ CEACAM + purity",
    "purity + histology": "{y} ~ CEACAM + purity + C(histology, Treatment('NOS'))",
}


def load():
    est = pd.read_csv(OUT / "tcga_estimate_scores.csv", index_col=0)
    est.index = est.index.astype(str).str[:12]
    epi = pd.read_csv(TCGA_BAYESPRISM_EPI,
                      sep="\t", index_col=0)
    epi.index = epi.index.astype(str).str[:12]
    common = epi.index.intersection(est.index)
    d = pd.DataFrame(index=common)
    d["CEACAM"] = np.log2(epi.loc[common, ["CEACAM5", "CEACAM6"]].mean(axis=1) + 1)
    d["ImmuneScore"] = est.loc[common, "ImmuneScore"]
    d["StromalScore"] = est.loc[common, "StromalScore"]
    ab = pd.read_csv(ABSOLUTE_PURITY, sep="\t")
    ab["patient"] = ab["array"].astype(str).str[:12]
    d["purity"] = (ab.dropna(subset=["purity"]).groupby("patient")["purity"]
                   .mean().reindex(d.index))
    cl = pd.read_csv(CLINICAL, sep="\t").drop_duplicates("submitter_id")
    d["diagnosis"] = cl.set_index("submitter_id")["primary_diagnosis"].reindex(d.index)
    d["histology"] = d["diagnosis"].map(HISTOLOGY).fillna("other")
    return d.dropna(subset=["CEACAM", "ImmuneScore", "purity"])


def main():
    d = load()
    rows = []
    for y in ("ImmuneScore", "StromalScore"):
        fits = [(name, "all", d, f) for name, f in MODELS.items()]
        fits += [("purity", g, d[d["histology"] == g], MODELS["purity"])
                 for g in ("intestinal", "diffuse")]
        for name, subset, data, formula in fits:
            m = smf.ols(formula.format(y=y), data).fit()
            rows.append({"outcome": y, "model": name, "subset": subset,
                         "n": int(m.nobs), "beta": m.params["CEACAM"],
                         "p": m.pvalues["CEACAM"]})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "tcga_histology_sensitivity.csv", index=False)
    print(d["histology"].value_counts().to_dict())
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
