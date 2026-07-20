#!/usr/bin/env python3
"""Single CLI entry point chaining identity_resolver -> evidence_harvester ->
[evidence_grader, failure_hunter] -> synthesizer into one command.

Usage:
    python3 run_pipeline.py --url <vendor_product_url>
    python3 run_pipeline.py --rrid <RRID>   # skip identity resolution

Each stage is called as a plain Python function, in-memory -- no temp files
on disk in normal operation. Pass --debug to also dump each stage's raw JSON
into the gitignored debug/ folder for troubleshooting.
"""
import argparse
import datetime
import json
import os
import sys
import urllib.error

import evidence_grader
import evidence_harvester
import failure_hunter
import identity_resolver
import synthesizer

OUTPUT_DIR = "output"
DEBUG_DIR = "debug"

NETWORK_ERRORS = (urllib.error.HTTPError, urllib.error.URLError, OSError, json.JSONDecodeError)


def _fail(stage: str, error: Exception):
    print(f"\n[FAILED] {stage}: {error}", file=sys.stderr)
    print("Pipeline stopped -- see error above.", file=sys.stderr)
    sys.exit(1)


def _dump_debug(name: str, data, timestamp: str):
    if not os.path.isdir(DEBUG_DIR):
        os.makedirs(DEBUG_DIR, exist_ok=True)
    path = os.path.join(DEBUG_DIR, f"{timestamp}_{name}.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"    (debug) wrote {path}")


def normalize_rrid(raw: str) -> str:
    raw = raw.strip()
    return raw if raw.upper().startswith("RRID:") else f"RRID:{raw}"


def build_identity_from_rrid(rrid: str) -> dict:
    return {
        "url": None,
        "rrid": normalize_rrid(rrid),
        "rrid_missing": False,
        "clone_id": None,
        "host_species": None,
        "clonality": None,
        "target_protein": None,
        "immunogen": None,
        "lot_variable": False,
        "raw_text": "",
    }


def run_harvest_with_progress(identity: dict, api_key: str) -> dict:
    """Mirrors evidence_harvester.harvest_evidence's two-tier logic, but with
    progress printed between the two queries -- reuses the same tested
    building blocks rather than duplicating the search/grading logic."""
    if identity.get("rrid_missing", True) or not identity.get("rrid"):
        return {
            "rrid": identity.get("rrid"),
            "evidence": [],
            "note": "no RRID — evidence harvest skipped, escalate to Tier 2 if needed.",
        }

    rrid = identity["rrid"]
    accession = rrid.split(":")[-1]
    domains = evidence_harvester.BASE_DOMAINS + evidence_harvester.VENDOR_VALIDATION_DOMAINS

    print("[2/5] Harvesting evidence (rrid_search)...")
    rrid_hits = evidence_harvester.search_rrid(rrid, api_key, domains)
    for hit in rrid_hits:
        hit["query_strategy"] = "rrid_search"

    confirmed_count = sum(
        evidence_harvester.accession_confirmed(h, accession) for h in rrid_hits
    )
    print(f"    -> {len(rrid_hits)} hits, {confirmed_count} confirmed against the accession")

    evidence = rrid_hits
    if confirmed_count < evidence_harvester.MIN_CONFIRMED_HITS:
        print("[2/5] RRID-search weak, running keyword_fallback...")
        keyword_hits = evidence_harvester.search_keyword_fallback(identity, api_key, domains)
        for hit in keyword_hits:
            hit["query_strategy"] = "keyword_fallback"
        evidence = evidence_harvester.dedupe_by_url(rrid_hits + keyword_hits)
        print(f"    -> {len(keyword_hits)} additional hits from keyword_fallback")
    else:
        print("    -> RRID-search sufficient, keyword_fallback not needed")

    return {"rrid": rrid, "evidence": evidence, "note": None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--url", help="Vendor antibody product page URL")
    group.add_argument("--rrid", help="Known RRID (skips identity resolution)")
    parser.add_argument(
        "--debug", action="store_true", help="Dump each stage's raw JSON to debug/"
    )
    args = parser.parse_args()

    identity_resolver.load_dotenv()
    api_key = os.environ.get("LINKUP_API_KEY")
    if not api_key:
        print("ERROR: LINKUP_API_KEY not set (check .env)", file=sys.stderr)
        sys.exit(1)

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # --- Stage 1: identity ---
    if args.url:
        print("[1/5] Resolving identity from vendor page...")
        try:
            identity = identity_resolver.resolve_identity(args.url, api_key)
        except NETWORK_ERRORS as e:
            _fail("identity_resolver", e)
    else:
        print("[1/5] Using provided RRID directly (skipping identity resolution)...")
        identity = build_identity_from_rrid(args.rrid)
        print(
            "    NOTE: target/clone unknown without a product page -- if RRID-search "
            "is weak, the keyword_fallback query will be less specific."
        )

    print(
        f"    -> RRID: {identity.get('rrid') or '(none)'} | "
        f"target: {identity.get('target_protein') or '?'} | "
        f"clone: {identity.get('clone_id') or '?'}"
    )
    if args.debug:
        _dump_debug("identity", identity, timestamp)

    if identity.get("rrid_missing"):
        print(
            "\n[STOPPED] No RRID found on this product page. Evidence harvesting "
            "requires an RRID and has been skipped -- nothing to grade or "
            "synthesize. This antibody may need Tier 2 (manual) review instead."
        )
        sys.exit(2)

    # --- Stage 2: evidence harvest ---
    try:
        harvest_result = run_harvest_with_progress(identity, api_key)
    except NETWORK_ERRORS as e:
        _fail("evidence_harvester", e)
    if args.debug:
        _dump_debug("harvest", harvest_result, timestamp)

    # --- Stage 3: grading ---
    print("[3/5] Grading evidence...")
    try:
        graded = evidence_grader.grade_evidence(
            harvest_result, api_key, identity.get("target_protein"), identity.get("clone_id")
        )
    except NETWORK_ERRORS as e:
        _fail("evidence_grader", e)
    print(f"    -> {len(graded.get('graded_evidence', []))} pieces of evidence graded")
    if args.debug:
        _dump_debug("graded", graded, timestamp)

    # --- Stage 4: failure hunting (independent of grading -- runs off identity,
    # not harvest_result, so it doesn't depend on stage 3's output; conceptually
    # parallel to grading even though this script calls it sequentially) ---
    print("[4/5] Hunting for negative evidence (independent of grading)...")
    try:
        failure_result = failure_hunter.hunt_failures(identity, api_key)
    except NETWORK_ERRORS as e:
        _fail("failure_hunter", e)
    print(f"    -> {failure_result['headline']}")
    if args.debug:
        _dump_debug("failure", failure_result, timestamp)

    # --- Stage 5: synthesis ---
    print("[5/5] Synthesizing scorecard...")
    scorecard = synthesizer.synthesize(identity, graded, failure_result)
    markdown = synthesizer.to_markdown(scorecard)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_path = os.path.join(OUTPUT_DIR, f"scorecard_{timestamp}.md")
    with open(output_path, "w") as f:
        f.write(markdown + "\n")

    print()
    print(markdown)
    print()
    print(f"Saved scorecard to {output_path}")


if __name__ == "__main__":
    main()
