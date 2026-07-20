#!/usr/bin/env python3
"""Standalone evidence grader.

Takes the JSON produced by evidence_harvester.py and assigns each piece of
evidence a validation tier per the Uhlen et al. (Nature Methods, 2016)
antibody-validation framework. Not wired into RocketRide.

Tiers (highest to lowest confidence):
  0        knockout/knockdown validated, independent of vendor
  1        orthogonal validation (RNA correlation, Human Protein Atlas, etc.)
  2        independent antibody agreement (concordant results across
           different antibodies against the same target)
  3        peer-reviewed use WITH an explicitly described control
  4        peer-reviewed use, no control described
  5        vendor-sourced data only (includes antibody-registry/citation-
           aggregator entries -- these restate vendor-submitted catalog data,
           not independent validation)
  negative documented failure / non-specificity report

Hard rule: usage alone never earns Tier 0-2, and nothing is inferred -- a
control must be explicitly described in the source text. Harvester snippets
are short excerpts, so the full page text is fetched via Linkup /fetch and
re-graded whenever either (a) the snippet already suggests Tier 0-2 language,
or (b) the hit is "directly relevant" -- its title/snippet names the
antibody's target protein or clone, not just a topically-adjacent paper --
since a snippet's silence on controls doesn't mean the full paper is silent
too. Every Tier 0-2 grade is backed by an actual quoted sentence from that
fetch, never a truncated snippet guess; hits that are neither snippet-flagged
nor directly relevant are graded from the snippet alone to bound cost.

Known scope gap: the source rubric grades per application + species
combination (e.g. "knockout-validated for WB in mouse" doesn't imply
validity for IHC in human). Neither identity_resolver.py nor
evidence_harvester.py currently capture which application/species a given
piece of evidence pertains to, so this grader assigns a single tier per
source without that finer-grained match -- a known limitation, not a
silent guess.
"""
import json
import os
import re
import sys
import urllib.error
import urllib.request

LINKUP_FETCH_URL = "https://api.linkup.so/v1/fetch"

REGISTRY_DOMAINS = ("antibodyregistry.org", "citeab.com")

NEGATIVE_PATTERNS = [
    r"did not work", r"does not work", r"failed to detect", r"non-specific",
    r"nonspecific", r"discontinued", r"no signal was detected",
    r"unable to validate", r"cross-react(?:ed|s|ivity)? with",
    r"recommend(?:ed)? against", r"pubpeer",
]
TIER0_PATTERNS = [
    r"knockout", r"knock-out", r"\bko\b", r"knockdown", r"knock-down",
    r"\bsirna\b", r"\bshrna\b", r"crispr", r"gene[- ]edited", r"null mutant",
    r"-/-\s*(cells|mice)",
]
TIER1_PATTERNS = [
    r"rna correlation", r"mrna expression correlat", r"human protein atlas",
    r"\bhpa\b", r"transcriptomic correlation", r"orthogonal validation",
]
TIER2_PATTERNS = [
    r"independent antibod", r"second(?:ary)? antibody (?:confirm|validat)",
    r"another antibody", r"concordant", r"multiple antibodies",
    r"cross-validated with",
]
TIER3_CONTROL_PATTERNS = [
    r"isotype control", r"negative control", r"peptide competition",
    r"blocking peptide", r"secondary-only control", r"loading control",
    r"positive control", r"validated? (?:in|with|using) ",
]


def load_dotenv(path=".env"):
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            key, value = key.strip(), value.strip()
            if key and value and key not in os.environ:
                os.environ[key] = value


def _matches_any(patterns, text):
    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            return match
    return None


def _is_registry_domain(url: str) -> bool:
    return any(domain in (url or "") for domain in REGISTRY_DOMAINS)


def _sentence_containing(text: str, match_span, context: int = 200) -> str:
    start, end = match_span
    window_start = max(0, start - context)
    window_end = min(len(text), end + context)
    return text[window_start:window_end].strip()


# NCBI/PMC intermittently serves a bot-challenge page instead of the article
# through Linkup's renderJs fetch (confirmed: ~1 in 3 attempts on the same
# URL); retrying recovers the real content, so a single transient block
# shouldn't cost an antibody its correct tier.
BOT_CHALLENGE_MARKER = "checking your browser"
FETCH_RETRIES = 2


class BotChallengeError(Exception):
    """Raised when every fetch attempt for a URL returned a bot-challenge page
    instead of real content -- this must be treated the same as a network
    failure, not as successfully-fetched (but content-free) text. Silently
    grading a challenge page as real content would default it to Tier 4
    ("no control described") instead of correctly falling back to the
    snippet-level grade -- confirmed to happen on unlucky runs where all
    FETCH_RETRIES+1 attempts hit the challenge page."""


def _fetch_once(url: str, api_key: str) -> str:
    payload = json.dumps({"url": url, "renderJs": True}).encode("utf-8")
    request = urllib.request.Request(
        LINKUP_FETCH_URL,
        data=payload,
        method="POST",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        body = json.loads(response.read().decode("utf-8"))
    return body.get("markdown", "")


def fetch_full_text(url: str, api_key: str) -> str:
    text = _fetch_once(url, api_key)
    attempts = 0
    while BOT_CHALLENGE_MARKER in text.lower() and attempts < FETCH_RETRIES:
        text = _fetch_once(url, api_key)
        attempts += 1
    if BOT_CHALLENGE_MARKER in text.lower():
        raise BotChallengeError(f"still bot-challenged after {attempts + 1} attempts: {url}")
    return text


def _mentions_this_antibody(text: str, rrid, clone_id) -> bool:
    """True only if the specific antibody (not just its target) is actually
    named in the text -- a paper can mention knockout/knockdown of a *gene*
    for an unrelated functional study without ever using this antibody to do
    it (confirmed: a TAF15-knockdown developmental-biology paper contained
    real knockout language but never mentioned this antibody's clone,
    catalog, or RRID anywhere)."""
    text = text.replace("\\_", "_").lower()
    if rrid:
        accession = rrid.split(":")[-1]
        if accession and accession.lower() in text:
            return True
    if clone_id and clone_id.lower() in text:
        return True
    return False


def _grade_text(text: str, source_label: str, rrid=None, clone_id=None) -> dict:
    """Run the tier cascade over a block of text (snippet or full page).

    If rrid/clone_id are given, any tier that would otherwise assert real
    validation/failure evidence (0-3, negative) additionally requires that
    this specific antibody -- not just its target -- is named in the same
    text; a control or knockout mention elsewhere in an unrelated paper isn't
    attributable to this reagent. Gets skipped when no identifiers are
    available to check against.
    """
    identifier_confirmed = (
        _mentions_this_antibody(text, rrid, clone_id) if (rrid or clone_id) else True
    )

    neg_match = _matches_any(NEGATIVE_PATTERNS, text)
    if neg_match and identifier_confirmed:
        return {
            "tier": "negative",
            "reason": f"documented failure/non-specificity language in {source_label}",
            "quote": _sentence_containing(text, neg_match.span()),
        }

    any_strong_match = None
    for tier, patterns, label in [
        (0, TIER0_PATTERNS, "knockout/knockdown language"),
        (1, TIER1_PATTERNS, "orthogonal validation language"),
        (2, TIER2_PATTERNS, "independent antibody agreement language"),
    ]:
        match = _matches_any(patterns, text)
        if match:
            any_strong_match = any_strong_match or match
            if identifier_confirmed:
                return {
                    "tier": tier,
                    "reason": f"{label} found in {source_label}",
                    "quote": _sentence_containing(text, match.span()),
                }

    control_match = _matches_any(TIER3_CONTROL_PATTERNS, text)
    if control_match and identifier_confirmed:
        return {
            "tier": 3,
            "reason": f"peer-reviewed use with an explicitly described control ({source_label})",
            "quote": _sentence_containing(text, control_match.span()),
        }

    if not identifier_confirmed and (neg_match or any_strong_match or control_match):
        return {
            "tier": 4,
            "reason": (
                f"{source_label} contains control/failure/knockout language, but never "
                "names this antibody's RRID/clone -- not attributable to it"
            ),
            "quote": None,
        }

    return {
        "tier": 4,
        "reason": f"peer-reviewed use, no control described ({source_label})",
        "quote": None,
    }


def grade_from_snippet(hit: dict, rrid=None, clone_id=None) -> dict:
    """First-pass tier signal from the harvester's title/snippet alone."""
    if _is_registry_domain(hit.get("url", "")):
        return {
            "tier": 5,
            "reason": "vendor/registry catalog entry, not independent validation",
            "quote": None,
        }
    text = f"{hit.get('title') or ''} {hit.get('snippet') or ''}"
    return _grade_text(text, "snippet", rrid, clone_id)


def _is_directly_relevant(hit: dict, target_protein, clone_id) -> bool:
    """True if the hit's own title/snippet names this antibody's target or
    clone -- i.e. it's plausibly about this specific antibody, not merely
    topically adjacent (e.g. a paper about a different member of the same
    protein family)."""
    text = f"{hit.get('title') or ''} {hit.get('snippet') or ''}".lower()
    if target_protein and target_protein.lower() in text:
        return True
    if clone_id and clone_id.lower() in text:
        return True
    return False


def confirm_with_full_text(hit: dict, snippet_grade: dict, api_key: str, rrid=None, clone_id=None) -> dict:
    """Re-grade a hit against its full page text so the tier is backed by an
    actual quoted sentence, not a snippet guess (whether the snippet already
    suggested Tier 0-2, or the hit was flagged as directly relevant and thus
    worth checking even though its snippet alone only supported Tier 3/4)."""
    try:
        full_text = fetch_full_text(hit["url"], api_key)
    except (urllib.error.HTTPError, urllib.error.URLError, OSError, json.JSONDecodeError, BotChallengeError):
        return snippet_grade | {
            "reason": snippet_grade["reason"]
            + " -- full text could not be fetched to confirm; kept snippet-level grade",
        }

    full_grade = _grade_text(full_text, "full text", rrid, clone_id)

    # Never let a full-text check downgrade a documented failure already
    # found in the snippet -- a negative signal that was explicit enough to
    # catch there shouldn't be erased just because the full-text scan missed it.
    if snippet_grade["tier"] == "negative" and full_grade["tier"] != "negative":
        return snippet_grade

    return full_grade


def grade_evidence(harvest_result: dict, api_key: str, target_protein=None, clone_id=None) -> dict:
    if harvest_result.get("note"):
        return {
            "rrid": harvest_result.get("rrid"),
            "graded_evidence": [],
            "note": harvest_result["note"],
        }

    hits = harvest_result.get("evidence", [])
    if not hits:
        return {
            "rrid": harvest_result.get("rrid"),
            "graded_evidence": [],
            "note": "no evidence found — search ran but returned nothing usable",
        }

    rrid = harvest_result.get("rrid")
    graded = []
    for hit in hits:
        grade = grade_from_snippet(hit, rrid, clone_id)
        is_registry = _is_registry_domain(hit.get("url", ""))
        snippet_suggests_strong_tier = grade["tier"] in (0, 1, 2)
        directly_relevant = not is_registry and _is_directly_relevant(hit, target_protein, clone_id)
        if not is_registry and (snippet_suggests_strong_tier or directly_relevant):
            grade = confirm_with_full_text(hit, grade, api_key, rrid, clone_id)
        graded.append(
            {
                "url": hit.get("url"),
                "title": hit.get("title"),
                "query_strategy": hit.get("query_strategy"),
                "tier": grade["tier"],
                "reason": grade["reason"],
                "quote": grade.get("quote"),
            }
        )

    return {"rrid": harvest_result.get("rrid"), "graded_evidence": graded, "note": None}


# Fallback demo harvest result (as if produced by evidence_harvester.py),
# used when no input is piped in.
DEMO_HARVEST_RESULT = {
    "rrid": "RRID:AB_2223172",
    "evidence": [
        {
            "title": "Antibody-Drug Conjugates Targeting the EGFR Ligand Epiregulin - PMC",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11875910/",
            "snippet": (
                "Primary antibodies used in this study include anti-β-actin "
                "(Cat# 4970, RRID:AB_2223172, 1:5000) ..."
            ),
            "query_strategy": "rrid_search",
        },
        {
            "title": "Cell Signaling Technology Cat# 4970, RRID:AB_2223172",
            "url": "https://www.antibodyregistry.org/AB_2223172",
            "snippet": "Antibodies: Search, Submit and Trace | Researchers only",
            "query_strategy": "rrid_search",
        },
    ],
    "note": None,
}


DEMO_TARGET_PROTEIN = "Beta Actin"
DEMO_CLONE_ID = "13E5"


def main():
    load_dotenv()
    api_key = os.environ.get("LINKUP_API_KEY")
    if not api_key:
        print("ERROR: LINKUP_API_KEY not set (check .env)", file=sys.stderr)
        sys.exit(1)

    target_protein, clone_id = None, None
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            harvest_result = json.load(f)
        if len(sys.argv) > 2:
            with open(sys.argv[2]) as f:
                identity = json.load(f)
            target_protein = identity.get("target_protein")
            clone_id = identity.get("clone_id")
    else:
        stdin_data = sys.stdin.read() if not sys.stdin.isatty() else ""
        if stdin_data.strip():
            harvest_result = json.loads(stdin_data)
        else:
            harvest_result = DEMO_HARVEST_RESULT
            target_protein, clone_id = DEMO_TARGET_PROTEIN, DEMO_CLONE_ID

    result = grade_evidence(harvest_result, api_key, target_protein, clone_id)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
