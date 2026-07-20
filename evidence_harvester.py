#!/usr/bin/env python3
"""Standalone evidence harvester.

Takes the JSON object produced by identity_resolver.py and, when it has a
valid RRID, calls Linkup /search scoped to RRID-registry and literature
domains to pull raw citation-lookup hits. Not wired into RocketRide.

This stage does not judge relevance -- it returns the raw hits so a later
stage can filter/score them.
"""
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

LINKUP_SEARCH_URL = "https://api.linkup.so/v1/search"

# RRID-registry / literature domains searched by default.
BASE_DOMAINS = ["ncbi.nlm.nih.gov", "biorxiv.org", "researchgate.net"]

# Vendor-validation domains: RRID registries and antibody-citation aggregators
# that independently confirm an antibody's identity, distinct from the
# vendor's own product page.
VENDOR_VALIDATION_DOMAINS = ["antibodyregistry.org", "citeab.com"]

# Below this many *confirmed* hits (accession literally present in the hit's
# title/snippet, not just domain-scoped noise) the RRID-only search is judged
# weak, and a keyword fallback is run alongside it.
MIN_CONFIRMED_HITS = 2

VENDOR_NAMES_BY_DOMAIN = {
    "abcam.com": "Abcam",
    "cellsignal.com": "Cell Signaling Technology",
    "thermofisher.com": "Thermo Fisher Scientific",
}


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


def _run_search(q: str, api_key: str, domains: list) -> list:
    payload = json.dumps(
        {
            "q": q,
            "depth": "standard",
            "outputType": "searchResults",
            "includeDomains": domains,
        }
    ).encode("utf-8")
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
        body = json.loads(response.read().decode("utf-8"))
    return [
        {
            "title": hit.get("name"),
            "url": hit.get("url"),
            "snippet": hit.get("content"),
        }
        for hit in body.get("results", [])
    ]


def search_rrid(rrid: str, api_key: str, domains: list) -> list:
    """Call Linkup POST /search for the exact RRID string, return raw hits.

    Queries the bare accession (e.g. "AB_2633323"), not the full "RRID:AB_..."
    string -- citing papers format the "RRID:" prefix inconsistently ("RRID: X"
    with a space, "(RRID:X)", etc.) but the accession itself always appears
    verbatim. The accession is quoted as an exact phrase -- an unquoted query
    returns semantically-similar-but-wrong RRIDs from antibodyregistry.org
    instead of the literal match (confirmed against a known heavily-cited RRID).
    """
    accession = rrid.split(":")[-1]
    return _run_search(f'"{accession}"', api_key, domains)


def vendor_from_url(url: str) -> str:
    host = urllib.parse.urlparse(url).netloc.lower().removeprefix("www.")
    for domain, name in VENDOR_NAMES_BY_DOMAIN.items():
        if host == domain or host.endswith("." + domain):
            return name
    return host


def search_keyword_fallback(identity: dict, api_key: str, domains: list) -> list:
    """Topical query built from identity fields, no RRID required.

    RRID-exact search only reliably surfaces literature for mega-cited
    antibodies (confirmed against a known heavily-cited control); for
    moderately-cited antibodies it returns registry noise even though
    relevant literature exists in the corpus, findable by target/clone/vendor
    keywords instead (confirmed against real TAF15 antibody candidates).
    """
    vendor = vendor_from_url(identity["url"]) if identity.get("url") else None
    terms = [identity.get("target_protein"), "antibody", identity.get("clone_id"), vendor]
    q = " ".join(t for t in terms if t)
    return _run_search(q, api_key, domains)


def accession_confirmed(hit: dict, accession: str) -> bool:
    haystack = f"{hit.get('title') or ''} {hit.get('snippet') or ''}"
    return accession in haystack


def dedupe_by_url(hits: list) -> list:
    seen = set()
    deduped = []
    for hit in hits:
        url = hit.get("url")
        if url in seen:
            continue
        seen.add(url)
        deduped.append(hit)
    return deduped


def harvest_evidence(identity: dict, api_key: str) -> dict:
    if identity.get("rrid_missing", True) or not identity.get("rrid"):
        return {
            "rrid": identity.get("rrid"),
            "evidence": [],
            "note": "no RRID — evidence harvest skipped, escalate to Tier 2 if needed.",
        }

    rrid = identity["rrid"]
    accession = rrid.split(":")[-1]
    domains = BASE_DOMAINS + VENDOR_VALIDATION_DOMAINS

    rrid_hits = search_rrid(rrid, api_key, domains)
    for hit in rrid_hits:
        hit["query_strategy"] = "rrid_search"

    confirmed_count = sum(accession_confirmed(h, accession) for h in rrid_hits)

    evidence = rrid_hits
    if confirmed_count < MIN_CONFIRMED_HITS:
        keyword_hits = search_keyword_fallback(identity, api_key, domains)
        for hit in keyword_hits:
            hit["query_strategy"] = "keyword_fallback"
        evidence = dedupe_by_url(rrid_hits + keyword_hits)

    return {"rrid": rrid, "evidence": evidence, "note": None}


# Fallback demo identity (as if produced by identity_resolver.py) for a
# heavily-cited antibody, used when no input is piped in.
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
        result = harvest_evidence(identity, api_key)
    except urllib.error.HTTPError as e:
        result = {"error": f"HTTP {e.code}: {e.reason}"}
    except urllib.error.URLError as e:
        result = {"error": str(e.reason)}
    except (OSError, json.JSONDecodeError) as e:
        result = {"error": str(e)}

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
