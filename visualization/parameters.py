"""Shared parameters for the visualization scripts and notebooks."""

# Tools that produce standardized sample-match outputs.
TOOLS = (
    "Conpair",
    "CrosscheckFingerprints",
    "HYSYS",
    "NGSCheckmate",
    "ntsm",
    "OmicsPrint",
    "Peddy",
    "Somalier",
    "TimeAttackGenComp",
    "Vireo",
)

# Every tool selectable from the command line (BAMixChecker only has read depth 0 results).
ALL_TOOLS = ("BAMixChecker",) + TOOLS

# Tools for the pseudobulk experiments (missing/double samples, pathological pseudobulks).
PSEUDOBULK_TOOLS = ("CrosscheckFingerprints", "HYSYS", "NGSCheckmate", "Vireo")

# Tools used in the dendrogram analysis.
DENDROGRAM_TOOLS = tuple(tool for tool in TOOLS if tool != "TimeAttackGenComp")

# These result sources do not have a read-depth directory in the output layout.
READ_DEPTH_INDEPENDENT_TOOLS = frozenset({"BAMixChecker", "Conpair", "Somalier", "ntsm"})

TOOL_LABELS = {
    "Conpair": "Conpair",
    "CrosscheckFingerprints": "Crosscheck",
    "HYSYS": "HYSYS",
    "NGSCheckmate": "NGSCheckMate",
    "ntsm": "ntsm",
    "OmicsPrint": "OmicsPrint",
    "Peddy": "Peddy",
    "Somalier": "Somalier",
    "TimeAttackGenComp": "TimeAttack",
    "Vireo": "Vireo",
}

# tool -> (comparison operator, default threshold); "ge" means a score >= threshold
# is a predicted match, "le" means a score <= threshold is a predicted match.
PREDICTION_THRESHOLDS = {
    "Conpair": ("ge", 0.8),
    "OmicsPrint": ("ge", 1.7),
    "Peddy": ("le", 0.0),
    "Somalier": ("ge", 1.0),
    "TimeAttackGenComp": ("le", 0.05),
}


def format_threshold(value: float) -> str:
    """Format a threshold as used in standardized prediction file names."""
    return format(value, "g")
