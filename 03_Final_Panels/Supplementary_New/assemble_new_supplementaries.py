"""
Assemble Supplementary Figures S7-S11 for the revision.

Layouts are expressed in mm on a 183 mm double-column canvas, matching the
existing supplementary figures. Panel aspect ratios follow the figsize of the
generating scripts so nothing is stretched.

Run: python assemble_new_supplementaries.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_Config"))
from shared.svg_assembler import VectorAssembler  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "_assembled"
OUT.mkdir(exist_ok=True)

W = 183.0   # double-column width, mm
M = 6.0     # outer margin
G = 5.0     # gutter between panels

FIGURES = {
    "S7_Cohort_Statistics": dict(
        height=258.0,
        title=None,
        panels=[
            ("A", "S7_A/S7_A_cohort_design.svg", M, 6, W - 2 * M, 118),
            ("B", "S7_B/S7_B_twosided_forest.svg", M, 130, W - 2 * M, 74),
            ("C", "S7_C/S7_C_design_defence.svg", M, 210, W - 2 * M, 42),
        ],
    ),
    # E-G were added in the third revision pass: the CEACAM evidence read across
    # cohorts rather than one at a time (R1.3). They belong here rather than in a
    # new figure, beside the per-marker IHC panel a reader is already looking at.
    "S8_CEACAM_Metaprogram": dict(
        height=268.0,
        title=None,
        panels=[
            ("A", "S8_A/S8_A_metaprogram_four_groups.svg", M, 6, W - 2 * M, 52),
            ("B", "S8_B/S8_B_ceacam_single_double_positive.svg", M, 64, W - 2 * M, 40),
            ("C", "S8_C/S8_C_ihc_per_marker.svg", M, 108, 124, 36),
            ("D", "S8_D/S8_D_metaprogram_external_validation.svg",
             M, 150, W - 2 * M, 42),
            ("E", "S8_E/S8_E_convergence_forest.svg", M, 198, W - 2 * M, 34),
            ("F", "S8_F/S8_F_leave_one_out.svg", M, 238, 100, 26),
            ("G", "S8_G/S8_G_rna_protein_concordance.svg", 112, 238, 65, 26),
        ],
    ),
    "S9_Mechanism_Specificity": dict(
        height=230.0,
        title=None,
        panels=[
            ("A", "S9_A/S9_A_spatial_adjusted_models.svg", M, 6, W - 2 * M, 44),
            ("B", "S9_B/S9_B_spatial_density_stratified.svg", M, 56, 118, 42),
            ("C", "S9_C/S9_C_momac_lineage_scores.svg", M, 104, 84, 62),
            ("D", "S9_D/S9_D_momac_lineage_dotplot.svg", 96, 104, 81, 62),
            ("E", "S9_E/S9_E_nfkb_per_celltype.svg", M, 172, 100, 52),
            ("F", "S9_F/S9_F_epithelial_vs_myeloid_cytokines.svg", 112, 172, 65, 46),
        ],
    ),
    "S10_PreTx_and_Adaptive": dict(
        height=225.0,
        title=None,
        panels=[
            ("A", "S10_A/S10_A_il1b_state_four_groups.svg", M, 6, 54, 44),
            ("B", "S10_B/S10_B_il1b_signature_four_groups.svg", 66, 6, 54, 44),
            ("C", "S10_C/S10_C_nfkb_pre_vs_post.svg", 126, 6, 51, 46),
            ("D", "S10_D/S10_D_adaptive_composition.svg", M, 58, W - 2 * M, 92),
            ("E", "S10_E/S10_E_adaptive_hallmark.svg", M, 156, W - 2 * M, 62),
        ],
    ),
    # The analyses added in the second revision pass, kept as their own figure
    # rather than crowded into S9, which already carries six panels.
    "S11_Affirmative_Analyses": dict(
        height=200.0,
        title=None,
        panels=[
            ("A", "S11_A/S11_A_spatial_positive_evidence.svg", M, 6, W - 2 * M, 42),
            ("B", "S11_B/S11_B_tcga_immune_exclusion.svg", M, 54, 118, 42),
            ("C", "S11_C/S11_C_nfkb_regulon_and_feedback.svg", M, 102, W - 2 * M, 44),
            ("D", "S11_D/S11_D_celltypist_annotation.svg", M, 152, W - 2 * M, 44),
        ],
    ),
}


def main():
    for name, spec in FIGURES.items():
        print(f"\n{name}")
        base = HERE / name
        asm = VectorAssembler(W, spec["height"], title=spec["title"])
        missing = []
        for label, rel, x, y, w, h in spec["panels"]:
            path = base / rel
            if not (path.exists() or path.with_suffix(".png").exists()):
                missing.append(rel)
                continue
            asm.place_panel(label, path, x, y, w, h)
        if missing:
            print(f"  MISSING panels, figure incomplete: {missing}")
        asm.save(name, output_dir=OUT)

    print(f"\nAssembled figures written to {OUT}")
    for f in sorted(OUT.iterdir()):
        print(f"   {f.name:<44} {f.stat().st_size / 1e6:7.2f} MB")


if __name__ == "__main__":
    main()
