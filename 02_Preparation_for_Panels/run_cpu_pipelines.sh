#!/usr/bin/env bash
# The filesystem roots below were specific to the machines this pipeline
# was run on. They are replaced by /path/to/machine and /path/to/home
# placeholders for deposition; repoint them at your own storage before
# running this script.
# Paths below refer to the upstream Round_4 processing pipeline, which is
# not part of this release. This script is included as a record of how the
# input was produced; it is not called by _run_all_panels.sh.
# Re-run the three CPU upstream pipelines, pinned so they reproduce themselves.
#
# These do not go to the cluster, and the reason is not preference. The user's
# PVCs map rsrch8 home and scratch only; nothing mounts
# /path/to/machine/hliang1_group, where every input to these pipelines lives. A pod
# would start and find no data. This login node has 48 cores and 754 GB, and
# BayesPrism 2.2.2 - the version PROVENANCE.txt names - and nichenetr 2.2.1 are
# already installed here, so the cluster would cost a staging copy and buy
# nothing. GraphST is the only one that needs a GPU, and it is the only one
# worth staging for.
#
# Threads are pinned, not maximised. A BLAS that varies its thread count
# reorders its sums and the result moves in the last few digits; that is how
# the NMF step came to disagree with its own deposited output by 6.4e-12. A
# fixed count is reproducible on the same machine, so each pipeline gets a
# fixed budget and they run side by side rather than each grabbing 48.
#
# Usage: bash run_cpu_pipelines.sh [bayesprism_tcga|bayesprism_tiger|nichenet|all]

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT=/path/to/Project_4_05232025
PREP="$PROJECT/Round_5/02_Preparation_for_Panels"
# Nothing uses TEMP since the BayesPrism TCGA migration of 2026-09-03. Kept
# only so the comment below has something to name.
TEMP="$PROJECT/Round_5/98_Temp_workspace"
LOGS="$HERE/logs/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGS"

# One fixed budget each, well under the 48 available, so the three together
# leave the node usable and none of them varies its own thread count run to run.
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export R_DATATABLE_NUM_THREADS=8

# step1_prepare_reference.py finds paths.py as parents[2]/00_Config, which only
# resolves at the depth it was written at. Staged into a rerun directory it
# raises ModuleNotFoundError instead. The script is left as found and the
# environment is what moves.
export PYTHONPATH="$PROJECT/Round_7_major_revision/00_Config${PYTHONPATH:+:$PYTHONPATH}"

note() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

# Each pipeline runs in its own copy of the deposited directory. The originals
# are inputs to the panels that ship; a re-run that wrote into them would
# overwrite the very artefacts it is meant to be compared against.
#
# -L, and it is the whole point. The TCGA staging directory reaches two of its
# files by symlink into the deposited directory. A plain `cp -a` copies those as
# links, so the re-run's own output would have been written straight through
# into the artefact it was supposed to be compared against. Dereferencing turns
# every staged file into a real one and the copy becomes self-contained.
stage() {
    local src="$1" name="$2"
    local dst="$HERE/$name/rerun"
    rm -rf "$dst"
    mkdir -p "$dst"
    cp -aL "$src"/*.tsv "$dst"/ 2>/dev/null
    cp -aL "$src"/*.py "$src"/*.R "$dst"/ 2>/dev/null
    echo "$dst"
}

# Migrated 2026-09-03 on the author's ruling. This staged from
# $TEMP/02162026_BayesPrism, a temp workspace that can be cleared at any time,
# and check 11 of verify_panel_provenance.py could not see it: its ABS_PATH
# regex only matches a literal /path/to/machine/... path, and this one is built out
# of $TEMP. The nine data files were copied into the pipeline's declared home,
# every md5 checked on both sides, and the temp workspace left exactly as it
# was. It is no longer on any path this chain needs.
#
# stage() copies *.tsv, *.py and *.R out of the source. The prepared directory
# holds data only - shipping pipeline scripts under 02_Preparation_for_Panels
# would put them in the code release, which syncs *.py/*.R from there - so the
# code comes from scripts/, the imported copy that PROVENANCE.md fingerprints.
# That also stops the two pinned scripts living only inside rerun/, which
# stage() deletes on its way in.
run_bayesprism_tcga() {
    note "BayesPrism TCGA: staging"
    local d
    d=$(stage "$PREP/BayesPrism_TCGA" BayesPrism_TCGA)
    cp -aL "$HERE/BayesPrism_TCGA/scripts"/*.py "$HERE/BayesPrism_TCGA/scripts"/*.R "$d"/
    for f in sc_counts.tsv sc_cell_types.tsv; do
        [ -s "$d/$f" ] || { note "BayesPrism TCGA: $f missing, cannot run"; return 2; }
    done
    note "BayesPrism TCGA: running step 2 (8 threads)"
    (cd "$d" && conda run -n r_demo Rscript step2_run_bayesprism_tcga.R) \
        > "$LOGS/bayesprism_tcga.log" 2>&1
    note "BayesPrism TCGA: step 2 exit $?"

    # Step 4 pulls the two per-cell-type tables the panels read out of the
    # posterior. The pinned copy differs from the deposited script in one line:
    # the original setwd()s to an absolute path in the deposited directory and
    # would write its output over the files it is being compared against.
    note "BayesPrism TCGA: running step 4"
    (cd "$d" && conda run -n r_demo Rscript step4_extract_epi_cd8_pinned.R) \
        > "$LOGS/bayesprism_tcga_step4.log" 2>&1
    note "BayesPrism TCGA: step 4 exit $?"
}

# Step 1 builds the bulk matrix from the TCGA FPKM download. It runs in its own
# directory because it rewrites sc_counts.tsv in place and skips the copy that
# produced it, so it is only correct starting from nothing.
run_bayesprism_tcga_step1() {
    local d="$HERE/BayesPrism_TCGA/rerun_step1"
    mkdir -p "$d"
    # From scripts/, not rerun/: rerun/ is wiped by stage() on every run, so
    # the pinned copy that only lived there was one pipeline run from gone.
    cp -aL "$HERE/BayesPrism_TCGA/scripts/step1_prepare_tcga_bulk"*.py "$d"/ 2>/dev/null
    rm -f "$d/sc_counts.tsv" "$d/sc_cell_types.tsv"
    note "BayesPrism TCGA: running step 1"
    (cd "$d" && conda run -n Liudeng_Python_310 python step1_prepare_tcga_bulk_pinned.py) \
        > "$LOGS/bayesprism_tcga_step1.log" 2>&1
    note "BayesPrism TCGA: step 1 exit $?"
}

run_bayesprism_tiger() {
    note "BayesPrism Tiger: staging"
    local d
    d=$(stage "$PREP/BayesPrism" BayesPrism_Tiger)
    # The Tiger directory kept its bulk matrix but not the single-cell
    # reference; step1 regenerates it from the h5ads.
    if [ ! -s "$d/sc_counts.tsv" ]; then
        note "BayesPrism Tiger: sc_counts.tsv absent, running step1 first"
        (cd "$d" && conda run -n Liudeng_Python_310 python step1_prepare_reference.py) \
            > "$LOGS/bayesprism_tiger_step1.log" 2>&1
    fi
    [ -s "$d/sc_counts.tsv" ] || { note "BayesPrism Tiger: step1 produced no reference"; return 2; }
    note "BayesPrism Tiger: running (8 threads)"
    (cd "$d" && conda run -n r_demo Rscript step2_run_bayesprism.R) \
        > "$LOGS/bayesprism_tiger.log" 2>&1
    note "BayesPrism Tiger: exit $?"
}

# Two steps, both named rather than globbed. `find .../scripts -name "*.R" |
# head -1` returned the archived Round_4 prototype 01_Prototype/
# run_nichenet_analysis.R, which is a parameterised entry point: with no
# arguments it prints a usage line and exits 1, which is the entire content of
# logs/20260830_225318/nichenet.log. It is still on disk beside the real
# scripts, so the glob would still pick it.
#
# stage() above flattens into one directory. These two scripts self-locate, read
# <app_dir>/01_Config/config.yaml and write <app_dir>/prepared_data and
# <app_dir>/results, taking no arguments at all, so the layout has to survive the
# copy and this one stages the tree instead. Staging is needed for the usual
# reason: run in place, step 2 would write over the deposited results/ it is
# being compared against.
run_nichenet() {
    note "NicheNet: staging"
    local d="$HERE/NicheNet/rerun"
    rm -rf "$d"
    mkdir -p "$d"
    cp -aL "$HERE/NicheNet/scripts/01_Config" "$HERE/NicheNet/scripts/02_Scripts" "$d"/
    [ -s "$d/01_Config/config.yaml" ] || { note "NicheNet: config not imported"; return 2; }
    note "NicheNet: running step 1 (extract .raw matrices)"
    conda run -n Liudeng_Python_310 python "$d/02_Scripts/01_prepare_data.py" \
        > "$LOGS/nichenet_step1.log" 2>&1
    note "NicheNet: step 1 exit $?"
    [ -s "$d/prepared_data/mtx_paths.yaml" ] || { note "NicheNet: step 1 produced no matrices"; return 2; }
    note "NicheNet: running step 2 (ligand activities, 8 threads)"
    conda run -n r_demo Rscript "$d/02_Scripts/02_run_nichenet_c3mac.R" \
        > "$LOGS/nichenet.log" 2>&1
    note "NicheNet: exit $?"
}

target="${1:-all}"
note "logs -> ${LOGS#"$PROJECT/"}"
case "$target" in
    bayesprism_tcga)  run_bayesprism_tcga_step1; run_bayesprism_tcga ;;
    bayesprism_tiger) run_bayesprism_tiger ;;
    nichenet)         run_nichenet ;;
    all)
        run_bayesprism_tcga &
        run_bayesprism_tiger &
        run_nichenet &
        wait
        ;;
    *) echo "usage: run_cpu_pipelines.sh [bayesprism_tcga|bayesprism_tiger|nichenet|all]"; exit 1 ;;
esac
note "done; logs in $LOGS"
