"""The Figure 4 B node letters, and the legend sentence that defines them.
"""

# --- Figure 4 B: the node letters -------------------------------------------
#
# A network node is labelled by one letter and one digit so the label fits its
# circle at 6 pt. The digit is the state's cluster index; the letter comes from
# this table, keyed by the lineage token (second underscore field) of the state
# name. Letters are not always initials because several lineages share one.
NODE_LETTER = {
    "CD4": "H",       # helper T
    "CD8": "C",       # cytotoxic T
    "NKT": "K",
    "NK": "K",
    "Prolif": "K",    # the one proliferating state in IM-T/NK/DC is NK
    "Mono": "O",
    "Mac": "M",
    "MoMac": "M",
    "Mast": "S",
    "DC": "D",
    "Neu": "N",
    "B": "B",
    "Plasma": "P",    # plasmablast
}
# The key printed in the Figure 4 legend; checked against NODE_LETTER below.
NODE_LETTER_LEGEND = ("H, CD4+ T; C, CD8+ T; K, NK/NKT; O, monocyte; "
                      "M, macrophage; S, mast; D, dendritic; N, neutrophil; "
                      "B, B cell; P, plasmablast")
_LEGEND_LETTERS = {tok.split(",")[0].strip() for tok in NODE_LETTER_LEGEND.split(";")}
if _LEGEND_LETTERS != set(NODE_LETTER.values()):
    raise ValueError(f"NODE_LETTER_LEGEND names {sorted(_LEGEND_LETTERS)} but "
                     f"NODE_LETTER uses {sorted(set(NODE_LETTER.values()))}")
