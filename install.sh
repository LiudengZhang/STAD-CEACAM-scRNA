#!/bin/bash
# Install the environment for this study. One command.
#
# Step 1 creates the conda environment from environment.yml. Step 2 installs
# the pip packages under a constraints file: pip resolves the whole graph and
# installs every dependency, but may not choose any version other than the one
# named, so nothing pinned moves. Step 3 imports every package; step 4 exits
# non-zero if a pinned version moved.
set -euo pipefail

ENV_NAME="${ENV_NAME:-stad_ceacam}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Read out of the environment itself, not copied from environment.yml, so the
# comparison is against what is installed rather than what was requested.
pins() {
    conda run -n "$ENV_NAME" python -c "
import importlib.metadata as m
for p in ['numpy', 'scanpy', 'anndata', 'matplotlib', 'seaborn']:
    try:
        print(p, m.version(p))
    except Exception:
        print(p, 'ABSENT')
"
}

echo ">>> 1/4  conda env create -n $ENV_NAME"
conda env create -n "$ENV_NAME" -f "$HERE/environment.yml"
before="$(pins)"

echo ">>> 2/4  pip install, under constraints"
# The constraints: every version this environment is pinned to. pip may install
# anything it needs and may not choose a different version of any of these.
CONSTRAINTS="$(mktemp)"
trap 'rm -f "$CONSTRAINTS"' EXIT
printf '%s
' 'PyComplexHeatmap==1.8.5' \
         'adjusttext==1.3.0' \
         'anndata2ri==1.3.2' \
         'anndata==0.10.9' \
         'arboreto==0.1.6' \
         'biopython==1.84' \
         'cairocffi==1.7.1' \
         'cairosvg==2.8.2' \
         'cellphonedb==5.0.1' \
         'cellrank==2.0.6' \
         'cnsplots==0.2.0' \
         'ctxcore==0.2.0' \
         'fbpca==1.0' \
         'geosketch==1.3' \
         'glasbey==0.3.0' \
         'gseapy==1.1.4' \
         'h5py==3.16.0' \
         'ktplotspy==0.3.2' \
         'lifelines==0.30.0' \
         'loompy==3.0.7' \
         'lxml==6.0.2' \
         'matplotlib-venn==1.1.1' \
         'matplotlib==3.10.5' \
         'mudata==0.3.8' \
         'natsort==8.4.0' \
         'networkx==3.4.2' \
         'num2tex==0.8' \
         'numba==0.65.1' \
         'numpy-groupies==0.11.2' \
         'numpy==1.26.4' \
         'openpyxl==3.1.5' \
         'palettable==3.3.3' \
         'pandas==2.3.3' \
         'patsy==1.0.1' \
         'pillow==12.2.0' \
         'plotnine==0.14.5' \
         'pycirclize==1.9.1' \
         'pymupdf==1.25.4' \
         'pyscenic==0.12.1' \
         'pytest==8.3.4' \
         'python-docx==1.2.0' \
         'pyyaml==6.0.3' \
         'requests==2.34.2' \
         'rpy2==3.5.11' \
         'scanpy==1.9.6' \
         'scikit-image==0.24.0' \
         'scikit-learn==1.7.2' \
         'scipy==1.15.3' \
         'seaborn==0.12.2' \
         'setuptools==82.0.1' \
         'statannotations==0.7.1' \
         'statsmodels==0.14.4' \
         'tqdm==4.68.3' \
         'upsetplot==0.9.0' > "$CONSTRAINTS"
conda run -n "$ENV_NAME" python -m pip install -c "$CONSTRAINTS" cnsplots==0.2.0 \
        glasbey==0.3.0 \
        num2tex==0.8 \
        PyComplexHeatmap==1.8.5 \
        cellphonedb==5.0.1 \
        geosketch==1.3 \
        ktplotspy==0.3.2 \
        fbpca==1.0 \
        biopython==1.84 \
        matplotlib-venn==1.1.1 \
        natsort==8.4.0 \
        numba==0.65.1 \
        palettable==3.3.3 \
        patsy==1.0.1 \
        statannotations==0.7.1 \
        upsetplot==0.9.0 \
        numpy-groupies==0.11.2 \
        pytest==8.3.4 \
        requests==2.34.2 \
        plotnine==0.14.5 \
        pycirclize==1.9.1 \
        colorspacious

echo ">>> 3/4  can every installed package be imported?"
# A declared package that cannot be imported is not installed in any sense that
# matters.
conda run -n "$ENV_NAME" python -c "
import importlib, sys
bad = []
for mod in ['cnsplots', 'glasbey', 'num2tex', 'PyComplexHeatmap', 'cellphonedb', 'geosketch', 'ktplotspy', 'fbpca', 'Bio', 'matplotlib_venn', 'natsort', 'numba', 'palettable', 'patsy', 'statannotations', 'upsetplot', 'numpy_groupies', 'pytest', 'requests', 'plotnine', 'pycirclize', 'colorspacious']:
    try:
        importlib.import_module(mod)
    except Exception as e:
        bad.append((mod, type(e).__name__, str(e)))
for mod, kind, msg in bad:
    print('   CANNOT IMPORT', mod, '-', kind + ':', msg, file=sys.stderr)
sys.exit(1 if bad else 0)
"

echo ">>> 4/4  did anything pinned move?"
after="$(pins)"
if [ "$before" = "$after" ]; then
    echo "$after" | sed "s/^/    unchanged: /"
    echo "environment $ENV_NAME ready.  conda activate $ENV_NAME"
else
    echo "REFUSED: the pip step moved a pinned version." >&2
    diff <(echo "$before") <(echo "$after") >&2 || true
    echo "Every figure in the paper was produced at the versions in" >&2
    echo "environment.yml. Do not use this environment." >&2
    exit 1
fi
