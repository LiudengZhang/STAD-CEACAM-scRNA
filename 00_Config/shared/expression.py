"""The one place an AnnData is asked which of its matrices is the expression.

WHY THIS MODULE EXISTS
----------------------
Two objects in this project answer that question differently, and until now
every reader answered it for itself.

  The analysis inputs, submission-tree/01_Raw_Inputs/01_H5AD/*.h5ad
      `.X`      dense, min -15.2, max 19.6, NaN in whole cell rows - the double
                normalisation of 2025-07-30, see 00_Data_Audit/FINDINGS.md
                sections 1 and 7
      `.raw.X`  csr float32, min 0.13, max 8.7, no NaN - the real log1p CP10K
                matrix, and what every panel in the paper was drawn from

  The deposited clean set, 06_Clean_Data/01_H5AD/*.h5ad
      `.X`      the same csr log1p matrix, promoted out of `.raw` by
                build_clean_h5ad.py
      `.raw`    ABSENT. Measured 2026-09-10: twelve of the thirteen files carry
                no `.raw` group at all.

So "read `.raw`, never `.X`" is right for the first and impossible for the
second, and a script that hard-codes either one is broken against the other.
Four deposited scripts hard-coded the first and raised `SystemExit` against the
deposit - a reviewer running the capsule with the Zenodo record attached got an
exception rather than a figure.

WHAT THIS IS NOT
----------------
It is not a fallback. `.raw` absent does not mean "use `.X` and hope"; the
object has to PROVE it is a promoted deposit before `.X` is read, and an object
that cannot prove it raises. That is RULES.md rule 6: the two matrices can be
mistaken for each other, so the mistake has to raise rather than quietly
returning the damaged one. The proof is three measurements, each of which the
damaged `.X` fails:

  1. `layers['counts']` is present. build_clean_h5ad.py and
     attach_counts_layer.py put the counts there; no analysis input has it
     under that name.
  2. `.X` is sparse. The damaged matrix is dense - scaling made it so.
  3. `.X` is finite and non-negative. The damaged matrix has NaN and runs to
     -15.2; a log1p CP10K matrix cannot be either.

Failing any of the three, the caller is told which and pointed at FINDINGS.md.

USE
---
    from shared.expression import expression_source, expression_adata

    src = expression_source(adata)          # has .X and .var_names
    i = list(src.var_names).index("CEACAM5")
    x = src.X[:, i]

    ad = expression_adata(adata)            # a full AnnData, obs carried

Both are exactly what the callers did by hand when `.raw` is present, so a run
against the analysis inputs produces the identical numbers it produced before.

Run `python expression.py` for the self-test, which shows the guard refusing a
damaged matrix before it is believed.
"""

from __future__ import annotations

__all__ = ["expression_source", "expression_adata", "promotion_problems",
           "self_test"]


def _describe(adata, name):
    return name or getattr(adata, "filename", None) or "this object"


def promotion_problems(adata):
    """Why `.X` may not be read as expression. Empty means it may.

    Returned rather than raised so a caller can report all three at once, and
    so the self-test can assert on them.
    """
    import numpy as np
    from scipy import sparse

    problems = []

    layers = getattr(adata, "layers", None)
    if layers is None or "counts" not in layers:
        problems.append(
            "layers['counts'] is absent - build_clean_h5ad.py and "
            "attach_counts_layer.py put the counts there, so an object without "
            "it is not the promoted deposit")

    X = adata.X
    if not sparse.issparse(X):
        problems.append(
            "`.X` is dense - the promoted matrix is csr float32; the dense one "
            "is the scaled matrix the double normalisation damaged")
    else:
        data = X.data
        if data.size:
            if not bool(np.isfinite(data).all()):
                problems.append(
                    "`.X` contains NaN or inf - a log1p CP10K matrix cannot; "
                    "see 00_Data_Audit/FINDINGS.md sections 1 and 7")
            if float(data.min()) < 0.0:
                problems.append(
                    f"`.X` has negative values (min {float(data.min()):.4g}) - "
                    f"a log1p CP10K matrix cannot")
    return problems


def _check_promoted(adata, name):
    problems = promotion_problems(adata)
    if problems:
        raise SystemExit(
            f"{_describe(adata, name)} has no `.raw`, and its `.X` cannot be "
            f"read as the log-normalised expression matrix:\n  - "
            + "\n  - ".join(problems)
            + "\nRefusing to read it. The analysis inputs keep the expression "
              "in `.raw`; the deposited clean h5ads promote it into `.X`. "
              "Nothing else is either.")


def expression_source(adata, name=None):
    """The object whose `.X` is the log-normalised expression matrix.

    `.raw` where the object has one - which is the documented rule and what
    every published number was produced from - and the object itself where it
    is a clean deposit that has proved it.
    """
    if getattr(adata, "raw", None) is not None:
        return adata.raw
    _check_promoted(adata, name)
    return adata


def expression_adata(adata, name=None):
    """The same, as a full AnnData with `obs` carried.

    For callers that go on to subset by an obs column rather than pull one
    gene's column out. `.raw.to_adata()` re-indexed by the parent's obs names
    is what those callers wrote by hand, so this changes nothing where `.raw`
    is present.
    """
    if getattr(adata, "raw", None) is not None:
        return adata.raw.to_adata()[adata.obs_names]
    _check_promoted(adata, name)
    return adata


# ------------------------------------------------------------------ self-test
def _fixtures():
    """Three objects: an analysis input, a clean deposit, and a damaged one."""
    import anndata
    import numpy as np
    import pandas as pd
    from scipy import sparse

    genes = ["CEACAM5", "CEACAM6", "ACTB"]
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(4)])
    logn = sparse.csr_matrix(np.array([[0.0, 1.2, 3.4],
                                       [2.1, 0.0, 0.5],
                                       [0.0, 0.0, 1.1],
                                       [4.2, 0.3, 0.0]], dtype="float32"))
    counts = sparse.csr_matrix((logn.toarray() * 10).astype("float32"))
    var = pd.DataFrame(index=genes)

    # 1. an analysis input: damaged dense .X, good .raw
    damaged = np.array([[-15.2, 1.0, 2.0],
                        [np.nan, 0.5, 1.0],
                        [3.0, np.nan, 0.0],
                        [1.0, 2.0, 19.6]], dtype="float32")
    an_input = anndata.AnnData(X=damaged.copy(), obs=obs.copy(), var=var.copy())
    an_input.raw = anndata.AnnData(X=logn.copy(), obs=obs.copy(),
                                   var=var.copy())

    # 2. a clean deposit: promoted .X, counts layer, no .raw
    deposit = anndata.AnnData(X=logn.copy(), obs=obs.copy(), var=var.copy(),
                              layers={"counts": counts.copy()})

    # 3. the damaged object with .raw stripped - the mistake rule 6 is about
    stripped = anndata.AnnData(X=damaged.copy(), obs=obs.copy(),
                               var=var.copy())
    return an_input, deposit, stripped, logn


def self_test():
    """Show the guard passing, then show it failing. Rule 4."""
    import numpy as np

    an_input, deposit, stripped, logn = _fixtures()
    fails = []

    # 1. An analysis input resolves to .raw, and to the SAME numbers the
    #    callers got by hand. This is the "changes nothing" claim, measured.
    src = expression_source(an_input)
    got = src.X.toarray()
    ok = np.array_equal(got, logn.toarray())
    print(f"  analysis input      -> .raw, matrix identical to .raw.X   "
          f"{'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append("analysis input did not resolve to .raw.X")

    ad = expression_adata(an_input)
    ok = (np.array_equal(ad.X.toarray(), logn.toarray())
          and list(ad.obs_names) == list(an_input.obs_names))
    print(f"  analysis input      -> expression_adata carries obs        "
          f"{'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append("expression_adata lost obs or changed the matrix")

    # 2. A clean deposit resolves to .X, and to the same numbers again.
    src = expression_source(deposit)
    ok = src is deposit and np.array_equal(src.X.toarray(), logn.toarray())
    print(f"  clean deposit       -> .X, same matrix                     "
          f"{'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append("clean deposit did not resolve to .X")

    # 3. MUTATION. The damaged matrix with .raw stripped must be REFUSED. If
    #    this returns anything, the guard is a fallback and every number
    #    downstream of it is drawn from a matrix that is 31% NaN.
    try:
        expression_source(stripped)
        print("  MUTANT damaged .X, no .raw: returned a matrix           FAIL")
        fails.append("mutation not caught: damaged .X accepted")
    except SystemExit as exc:
        hits = promotion_problems(stripped)
        print(f"  MUTANT damaged .X, no .raw: refused, {len(hits)} reasons   ok")
        for h in hits:
            print(f"      - {h[:66]}")
        # Two, not three: a dense `.X` disqualifies the object outright, so the
        # finiteness and sign tests - which live inside the sparse branch and
        # read `.data` - are not reached. That is deliberate. There is no dense
        # `.data` array to read, and the object has already failed.
        want = ("layers['counts'] is absent", "is dense")
        missing = [w for w in want if not any(w in h for h in hits)]
        if missing:
            fails.append(f"expected reasons {missing} among {hits}: {exc}")

    # 4. MUTATION. Each of the three proofs must be able to fail ALONE. A
    #    conjunction that only ever fails all-at-once is not three checks.
    import anndata
    from scipy import sparse
    for label, mutate, phrase in (
            ("counts layer removed",
             lambda a: anndata.AnnData(X=a.X.copy(), obs=a.obs.copy(),
                                       var=a.var.copy()),
             "layers['counts'] is absent"),
            ("`.X` densified",
             lambda a: anndata.AnnData(X=a.X.toarray(), obs=a.obs.copy(),
                                       var=a.var.copy(),
                                       layers={"counts": a.layers["counts"]}),
             "is dense"),
            ("one NaN injected",
             lambda a: _with_nan(a), "NaN or inf"),
            ("one negative injected",
             lambda a: _with_negative(a), "negative values")):
        mutant = mutate(deposit)
        hits = [h for h in promotion_problems(mutant) if phrase in h]
        caught = bool(hits)
        print(f"  MUTANT {label:<22} -> refused: "
              f"{'yes' if caught else 'NO'}")
        if not caught:
            fails.append(f"mutation not caught: {label}")

    # 5. And the unmutated deposit must still pass, so the four above are
    #    failing for their own reason and not because the fixture is broken.
    ok = not promotion_problems(deposit)
    print(f"  unmutated deposit still accepted                          "
          f"{'ok' if ok else 'FAIL'}")
    if not ok:
        fails.append("the fixture itself does not pass")

    return 1 if fails else 0


def _with_nan(a):
    import anndata
    import numpy as np
    X = a.X.copy()
    X.data = X.data.copy()
    X.data[0] = np.nan
    return anndata.AnnData(X=X, obs=a.obs.copy(), var=a.var.copy(),
                           layers={"counts": a.layers["counts"]})


def _with_negative(a):
    import anndata
    X = a.X.copy()
    X.data = X.data.copy()
    X.data[0] = -0.5
    return anndata.AnnData(X=X, obs=a.obs.copy(), var=a.var.copy(),
                           layers={"counts": a.layers["counts"]})


if __name__ == "__main__":
    raise SystemExit(self_test())
