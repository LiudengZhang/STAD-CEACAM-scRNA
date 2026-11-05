# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""Determinism pinning for the GraphST reproduction.

Loaded automatically by CPython at interpreter start (standard `site` mechanism)
because $GRAPHST_ROOT/tools/pysite is on PYTHONPATH.  This exists so that the
staged analysis scripts stay a byte-identical provenance copy: nothing here
edits them.

Active only when GRAPHST_PIN=1, so ordinary use of the env (pip, etc.) is
unaffected.  Every required variable is read without a default: an unset one
raises rather than silently pinning something different.
"""
import os

if os.environ.get("GRAPHST_PIN") == "1":
    n = int(os.environ["GRAPHST_NUM_THREADS"])
    seed = int(os.environ["GRAPHST_SEED"])

    import random
    import numpy as np
    import torch

    torch.set_num_threads(n)
    try:
        torch.set_num_interop_threads(n)
        interop = torch.get_num_interop_threads()
    except RuntimeError as exc:                       # already initialised
        interop = f"REFUSED: {exc}"

    # GraphST.preprocess.fix_seed re-seeds these at model construction with the
    # config's seed; this covers everything that happens before that point.
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    det_err = None
    if os.environ.get("GRAPHST_DETERMINISTIC") == "1":
        try:
            torch.use_deterministic_algorithms(True)
            det = torch.are_deterministic_algorithms_enabled()
        except Exception as exc:                      # recorded, never dropped
            det, det_err = False, f"{type(exc).__name__}: {exc}"
    else:
        det = False

    print("[pin] " + " ".join([
        f"torch={torch.__version__}",
        f"num_threads={torch.get_num_threads()}",
        f"interop_threads={interop}",
        f"deterministic_algorithms={det}",
        f"seed={seed}",
        f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}",
        f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS')}",
        f"MKL_NUM_THREADS={os.environ.get('MKL_NUM_THREADS')}",
        f"NUMEXPR_NUM_THREADS={os.environ.get('NUMEXPR_NUM_THREADS')}",
        f"PYTHONHASHSEED={os.environ.get('PYTHONHASHSEED')}",
    ]), flush=True)
    if det_err:
        print(f"[pin] use_deterministic_algorithms(True) FAILED verbatim: {det_err}", flush=True)
