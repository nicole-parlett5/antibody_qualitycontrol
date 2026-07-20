#!/usr/bin/env python3
"""Standalone negative-evidence hunter.

Separate agent from evidence_grader.py by design: its only job is surfacing
documented failure/non-specificity reports for an antibody, run independently
so negative signal doesn't get diluted or out-voted by a synthesis step that
also weighs positive evidence. Not wired into RocketRide.

Searches Linkup /search two ways:
  1. Domain-scoped to pubpeer.com, reddit.com, researchgate.net.
  2. General web search (no domain restriction) for "<identifier> discontinued"
     and "<identifier> failed".

Tone rule: this script reports, it does not judge. Output is always phrased
as "N reports of failure found: [links]", never as an assertion that the
antibody doesn't work. A single matching hit is explicitly labeled an
uncorroborated anecdote; only 2+ independent hits are called "corroborated".

Every candidate hit is also required to actually name this antibody's target
or clone in its own text before counting as a report -- a "failed" or
"discontinued" mention that never names this antibody isn't evidence about
it (the same false-positive class caught earlier in evidence_grader.py: a
paper about a gene/topic overlapping the search terms is not the same as a
paper about this specific reagent).
"""
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

LINKUP_SEARCH_URL = "https://api.linkup.so/v1/search"

FORUM_DOMAINS = ["pubpeer.com", "reddit.com", "researchgate.net"]

FAILURE_PATTERNS = [
    r"does(?:n't| not) work", r"did not work", r"didn't work",
    r"no signal", r"failed to detect", r"non-specific", r"nonspecific",
    r"cross-react(?:ed|s|ivity)?", r"false positive", r"unreliable",
    r"discontinued", r"withdrawn", r"avoid this antibody", r"waste of money",
    r"not specific", r"couldn't reproduce", r"could not reproduce",
    r"flagged (?:as|for)", r"pubpeer", r"non-selective",
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


def _run_search(q: str, api_key: str, domains=None) -> list:
    body = {"q": q, "depth": "standard", "outputType": "searchResults"}
    if domains:
        body["includeDomains"] = domains
    payload = json.dumps(body).encode("utf-8")
    request = urllib.request.Request(
        LINKUP_SEARCH_URL,
        data=payload,
        method="POST",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        result = json.loads(response.read().decode("utf-8"))
    return [
        {"title": hit.get("name"), "url": hit.get("url"), "snippet": hit.get("content")}
        for hit in result.get("results", [])
    ]


def build_identifier(identity: dict) -> str:
    vendor = None
    if identity.get("url"):
        host = urllib.parse.urlparse(identity["url"]).netloc.lower().removeprefix("www.")
        vendor = host.split(".")[0] if host else None
    parts = [identity.get("target_protein"), identity.get("clone_id"), vendor]
    return " ".join(p for p in parts if p)


def _matches_any(patterns, text):
    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            return match
    return None


def _sentence_containing(text: str, match_span, context: int = 150) -> str:
    start, end = match_span
    return text[max(0, start - context): min(len(text), end + context)].strip()


def _normalize_title(title) -> str:
    """Strip indexing-site noise so the same paper mirrored across PubMed,
    ResearchGate, PMC, etc. collapses to one comparable key."""
    if not title:
        return ""
    t = title.lower()
    t = re.sub(r"^\(pdf\)\s*", "", t)
    t = re.sub(r"\s*[-|]\s*(pubmed|pmc|researchgate|nih|ncbi)\.?\s*$", "", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def _names_this_antibody(text: str, target_protein, clone_id) -> bool:
    text = text.lower()
    if target_protein and target_protein.lower() in text:
        return True
    if clone_id and clone_id.lower() in text:
        return True
    return False


def hunt_failures(identity: dict, api_key: str) -> dict:
    identifier = build_identifier(identity)
    target_protein = identity.get("target_protein")
    clone_id = identity.get("clone_id")

    raw_hits = []
    forum_hits = _run_search(f'"{identifier}"', api_key, FORUM_DOMAINS)
    for hit in forum_hits:
        hit["found_via"] = "forum_scoped"
    raw_hits += forum_hits

    for suffix, tag in [("discontinued", "general_discontinued"), ("failed", "general_failed")]:
        hits = _run_search(f'"{identifier}" {suffix}', api_key)
        for hit in hits:
            hit["found_via"] = tag
        raw_hits += hits

    seen_urls = set()
    reports = []
    for hit in raw_hits:
        url = hit.get("url")
        if not url or url in seen_urls:
            continue
        text = f"{hit.get('title') or ''} {hit.get('snippet') or ''}"
        failure_match = _matches_any(FAILURE_PATTERNS, text)
        if not failure_match:
            continue
        if not _names_this_antibody(text, target_protein, clone_id):
            continue
        seen_urls.add(url)
        reports.append(
            {
                "source": urllib.parse.urlparse(url).netloc,
                "url": url,
                "title": hit.get("title"),
                "claim": _sentence_containing(text, failure_match.span()),
                "found_via": hit.get("found_via"),
            }
        )

    n = len(reports)
    domains = {r["source"] for r in reports}
    # Two hits can both be the *same* underlying paper mirrored across
    # indexes (e.g. PubMed + ResearchGate both pointing at one publication) --
    # that's one documented incident, not two independent observations, so
    # corroboration strength is judged on distinct underlying reports, not
    # raw hit count.
    distinct_titles = {_normalize_title(r["title"]) for r in reports}

    if n == 0:
        corroboration = "no reports found"
    elif len(distinct_titles) == 1:
        corroboration = (
            "anecdote, uncorroborated"
            if n == 1
            else f"{n} hits found, but all refer to the same underlying report "
            f"(indexed on {len(domains)} platform{'s' if len(domains) != 1 else ''}) "
            "-- not independent corroboration"
        )
    elif len(domains) >= 2:
        corroboration = (
            f"corroborated by {len(distinct_titles)} independent reports "
            f"across {len(domains)} platforms ({n} total hits)"
        )
    else:
        corroboration = (
            f"corroborated by {len(distinct_titles)} reports, but all from the "
            f"same platform ({next(iter(domains))})"
        )

    for r in reports:
        r["corroboration"] = corroboration

    if n == 0:
        headline = (
            "0 reports of failure found in the searched sources. "
            "This does not confirm the antibody performs well -- only that no "
            "failure reports surfaced in this search."
        )
    else:
        links = ", ".join(r["url"] for r in reports)
        headline = f"{n} reports of failure found: {links}"

    return {
        "identifier": identifier,
        "headline": headline,
        "corroboration": corroboration,
        "reports": reports,
    }


# Fallback demo identity: a real, published case of antibody failure --
# Proteintech 16211-1-AP (anti-N6AMT1) shown to cross-react with Aurora
# Kinase A (PMC11130794), used when no input is piped in.
DEMO_IDENTITY = {
    "url": "https://www.ptglab.com/products/N6AMT1-Antibody-16211-1-AP.htm",
    "rrid": None,
    "rrid_missing": True,
    "clone_id": "16211-1-AP",
    "host_species": "Rabbit",
    "clonality": "polyclonal",
    "target_protein": "N6AMT1",
    "immunogen": None,
    "lot_variable": True,
    "raw_text": "",
}


def main():
    load_dotenv()
    api_key = os.environ.get("LINKUP_API_KEY")
    if not api_key:
        print("ERROR: LINKUP_API_KEY not set (check .env)", file=sys.stderr)
        sys.exit(1)

    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            identity = json.load(f)
    else:
        stdin_data = sys.stdin.read() if not sys.stdin.isatty() else ""
        identity = json.loads(stdin_data) if stdin_data.strip() else DEMO_IDENTITY

    try:
        result = hunt_failures(identity, api_key)
    except urllib.error.HTTPError as e:
        result = {"error": f"HTTP {e.code}: {e.reason}"}
    except urllib.error.URLError as e:
        result = {"error": str(e.reason)}
    except (OSError, json.JSONDecodeError) as e:
        result = {"error": str(e)}

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
