#!/usr/bin/env bash
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
# Per-sample NMF, pinned so that it reproduces itself byte for byte.
#
# Two settings do the work, and both were established by measurement rather
# than read off the code:
#
#   nrun = 10   The scripts/ copy of 02_run_nmf_per_sample.R reads its config
#               from a hard-coded path that no longer resolves, and the two
#               candidate configs disagree on this one parameter - the
#               standardised pipeline says n_iter 10, the archived 3CA copy
#               says 100. Run both ways at seed 42 against the deposited basis
#               matrix (test_nrun.R):
#                   nrun 10   max |dW| = 6.4e-12, every program r = 1.0000
#                   nrun 100  max |dW| = 77.9,    r = 0.75 to 0.87
#               Ten is what made the deposited programs.
#
#   1 thread    At seed 42 the factorisation is deterministic in exact
#               arithmetic, but a multi-threaded BLAS reorders its sums and
#               leaves a few parts in 1e12. With the thread count pinned,
#               test_bitwise.R gets identical = TRUE and max |W1 - W2| = 0.
#               The residual 6.4e-12 against the deposited run is that same
#               effect in the original run, whose thread count was never
#               recorded.
#
# scripts/ holds the code exactly as it was found. Nothing there is edited;
# every choice needed to make it deterministic is made here, where it can be
# read.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_ENV="${STAD_R_ENV:-r_bayesprism}"

# Single-threaded BLAS. This is the difference between "reproduces to 1e-12"
# and "reproduces exactly", and it costs wall-clock, not correctness.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "NMF, pinned: seed 42, method brunet, nrun 10, 1 BLAS thread"
echo "R environment: $R_ENV"
echo

case "${1:-verify}" in
    verify)
        # Re-run a few samples and compare the gene lists the panels read.
        conda run -n "$R_ENV" Rscript "$HERE/verify_nmf.R" "${2:-3}"
        ;;
    settings)
        # Show why nrun is 10 and not 100.
        conda run -n "$R_ENV" Rscript "$HERE/test_nrun.R" "${2:-inLN_1221}" 4
        ;;
    bitwise)
        # Show that two runs of the pinned pipeline are identical.
        conda run -n "$R_ENV" Rscript "$HERE/test_bitwise.R" "${2:-inLN_1221}" 4
        ;;
    full)
        # The whole pipeline, all 54 samples, K = 4 to 9. Hours.
        echo "The scripts in scripts/ carry the original hard-coded paths and"
        echo "have not been repointed. Running the full sweep needs the config"
        echo "and input directory wired up first; see PROVENANCE.md."
        exit 2
        ;;
    *)
        echo "usage: run.sh [verify N | settings SAMPLE | bitwise SAMPLE | full]"
        exit 1
        ;;
esac
