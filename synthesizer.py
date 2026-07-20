#!/usr/bin/env python3
"""Standalone scorecard synthesizer.

Takes the identity JSON from identity_resolver.py and the graded-evidence
JSON from evidence_grader.py, and rolls them up into one antibody QC
scorecard: overall tier plus every piece of evidence annotated with its tier
and query_strategy, so the scorecard shows both the verdict and how each
piece of evidence was found. Not wired into RocketRide.
"""
import json
import sys

TIER_LABELS = {
    0: "Tier 0 — knockout/knockdown validated, independent of vendor",
    1: "Tier 1 — orthogonal validation (RNA correlation / HPA)",
    2: "Tier 2 — independent antibody agreement",
    3: "Tier 3 — peer-reviewed use with described control",
    4: "Tier 4 — peer-reviewed use, no control described",
    5: "Tier 5 — vendor/registry data only",
    "negative": "Negative — documented failure/non-specificity reported",
}

# Best (0) to worst (5); "negative" is tracked separately, not on this scale.
TIER_ORDER = [0, 1, 2, 3, 4, 5]


def synthesize(identity: dict, graded: dict, failure_result: dict = None) -> dict:
    evidence = graded.get("graded_evidence", [])
    negatives = [e for e in evidence if e["tier"] == "negative"]
    positives = [e for e in evidence if e["tier"] != "negative"]

    if not evidence:
        overall_tier = None
        overall_label = graded.get("note") or "no evidence found"
    elif positives:
        best = min(positives, key=lambda e: TIER_ORDER.index(e["tier"]))
        overall_tier = best["tier"]
        overall_label = TIER_LABELS[overall_tier]
    else:
        overall_tier = "negative"
        overall_label = TIER_LABELS["negative"]

    negative_evidence = None
    if failure_result is not None:
        negative_evidence = {
            "identifier": failure_result.get("identifier"),
            "headline": failure_result.get("headline"),
            "corroboration": failure_result.get("corroboration"),
            "reports": failure_result.get("reports", []),
        }

    return {
        "url": identity.get("url"),
        "rrid": identity.get("rrid"),
        "rrid_missing": identity.get("rrid_missing"),
        "target_protein": identity.get("target_protein"),
        "clone_id": identity.get("clone_id"),
        "host_species": identity.get("host_species"),
        "clonality": identity.get("clonality"),
        "lot_variable": identity.get("lot_variable"),
        "overall_tier": overall_tier,
        "overall_tier_label": overall_label,
        "flagged_negative_reports": negatives,
        "negative_evidence": negative_evidence,
        "evidence": [
            {
                "url": e["url"],
                "title": e["title"],
                "tier": e["tier"],
                "tier_label": TIER_LABELS[e["tier"]],
                "query_strategy": e["query_strategy"],
                "reason": e["reason"],
                "quote": e["quote"],
            }
            for e in evidence
        ],
    }


def _truncate(text, length=180):
    if not text:
        return ""
    text = text.replace("\n", " ").strip()
    return text if len(text) <= length else text[: length - 1].rstrip() + "…"


def _tier_sort_key(tier):
    # "negative" (documented failure) sorts first -- most important to flag --
    # then numeric tiers ascending (0 = strongest validation).
    return -1 if tier == "negative" else tier


def to_markdown(scorecard: dict) -> str:
    """Render a synthesize() result as the readable scorecard table."""
    lines = ["# Antibody QC Scorecard", ""]
    lines.append(f"**URL:** {scorecard.get('url') or '—'}")
    lines.append(f"**RRID:** {scorecard.get('rrid') or '—'}")
    lines.append(f"**Target protein:** {scorecard.get('target_protein') or '—'}")
    lines.append(f"**Clone:** {scorecard.get('clone_id') or '—'}")
    lines.append(f"**Host species:** {scorecard.get('host_species') or '—'}")
    lines.append(f"**Clonality:** {scorecard.get('clonality') or '—'}")
    lines.append(f"**Lot variable (polyclonal, no clone ID):** {'Yes' if scorecard.get('lot_variable') else 'No'}")
    lines.append("")
    lines.append(f"## Overall Tier: {scorecard.get('overall_tier_label') or '—'}")
    lines.append("")

    negative_evidence = scorecard.get("negative_evidence")
    if negative_evidence is not None:
        lines.append("## Negative Evidence Check (independent failure_hunter run)")
        lines.append("")
        lines.append(negative_evidence.get("headline") or "_Not run._")
        lines.append("")
        lines.append(f"**Corroboration:** {negative_evidence.get('corroboration') or '—'}")
        lines.append("")

    evidence = sorted(scorecard.get("evidence", []), key=lambda e: _tier_sort_key(e["tier"]))
    if not evidence:
        lines.append("_No positive evidence found._")
        return "\n".join(lines)

    lines.append("## Evidence")
    lines.append("")
    lines.append("| Tier | Source | Query Strategy | Reason | Quote |")
    lines.append("|------|--------|-----------------|--------|-------|")
    for e in evidence:
        tier_cell = str(e["tier"]).upper() if e["tier"] == "negative" else str(e["tier"])
        title = (e.get("title") or e["url"]).replace("|", "\\|")
        source_cell = f"[{_truncate(title, 70)}]({e['url']})"
        reason_cell = _truncate(e.get("reason"), 90).replace("|", "\\|")
        quote_cell = _truncate(e.get("quote"), 140).replace("|", "\\|")
        lines.append(
            f"| {tier_cell} | {source_cell} | {e.get('query_strategy') or '—'} | {reason_cell} | {quote_cell} |"
        )

    return "\n".join(lines)


# Fallback demo inputs, used when no file paths are given.
DEMO_IDENTITY = {
    "url": "https://www.cellsignal.com/products/primary-antibodies/b-actin-13e5-rabbit-mab/4970",
    "rrid": "RRID:AB_2223172",
    "rrid_missing": False,
    "clone_id": "13E5",
    "host_species": "Rabbit",
    "clonality": "monoclonal",
    "target_protein": "Beta Actin",
    "immunogen": None,
    "lot_variable": False,
}


def main():
    if len(sys.argv) > 2:
        with open(sys.argv[1]) as f:
            identity = json.load(f)
        with open(sys.argv[2]) as f:
            graded = json.load(f)
    else:
        identity = DEMO_IDENTITY
        # Run evidence_grader.py's own demo path for a self-contained default.
        import evidence_grader

        evidence_grader.load_dotenv()
        import os

        api_key = os.environ.get("LINKUP_API_KEY")
        if not api_key:
            print("ERROR: LINKUP_API_KEY not set (check .env)", file=sys.stderr)
            sys.exit(1)
        graded = evidence_grader.grade_evidence(
            evidence_grader.DEMO_HARVEST_RESULT,
            api_key,
            evidence_grader.DEMO_TARGET_PROTEIN,
            evidence_grader.DEMO_CLONE_ID,
        )

    print(json.dumps(synthesize(identity, graded, None), indent=2))


if __name__ == "__main__":
    main()
