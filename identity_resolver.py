#!/usr/bin/env python3
"""Standalone antibody identity resolver.

Fetches a vendor antibody product page via the Linkup /fetch API (with JS
rendering) and pulls out identity fields via best-effort regex matching over
the returned markdown. Not wired into RocketRide.

Extraction is heuristic: vendor page layouts vary, so any field that can't be
confidently matched is left as null rather than guessed. RRID in particular is
never fabricated -- if the pattern isn't found, rrid is null and
rrid_missing is set.
"""
import json
import os
import re
import sys
import urllib.error
import urllib.request

LINKUP_FETCH_URL = "https://api.linkup.so/v1/fetch"

RRID_PATTERN = re.compile(r"RRID[^\n]{0,15}?(AB_\d+)", re.IGNORECASE)
CLONE_LABEL_PATTERN = re.compile(r"Clone(?:\s*(?:ID|No\.?|Number))?\b\s*(.*)", re.IGNORECASE)
HOST_PATTERN = re.compile(
    r"Host[^\n]{0,25}?\b(Rabbit|Mouse|Rat|Goat|Sheep|Chicken|Guinea\s?Pig|"
    r"Hamster|Donkey|Human|Bovine)\b",
    re.IGNORECASE,
)
IMMUNOGEN_LABEL_PATTERN = re.compile(r"Immunogen\b\s*(.*)", re.IGNORECASE)
HEADING_PATTERN = re.compile(r"^#+\s*(.+)$", re.MULTILINE)
ANTI_TARGET_PATTERN = re.compile(
    r"Anti[- ]([A-Za-z0-9][A-Za-z0-9/\-\. ]{1,40}?)\s+[Aa]ntibody"
)
TITLE_NOISE_WORDS = re.compile(
    r"\b(Antibody|Monoclonal|Polyclonal|Recombinant|Rabbit|Mouse|Rat|Goat|"
    r"Sheep|Chicken|Human|mAb|pAb|Ab)\b",
    re.IGNORECASE,
)


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


def fetch_page(url: str, api_key: str) -> str:
    """Call Linkup POST /fetch with renderJs=true, return the markdown body."""
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


def _first_heading(text: str):
    match = HEADING_PATTERN.search(text)
    return match.group(1) if match else None


def _label_value(rest: str):
    """Given the text following a field label, pull the value whether the
    source used 'Label: value', 'Label value', or a markdown table row
    'Label | value |'."""
    line = rest.split("\n", 1)[0]
    value = re.sub(r"^[\s:#\|]+", "", line)
    value = value.split("|", 1)[0].strip()
    return value


def extract_clone_id(text: str):
    label_match = CLONE_LABEL_PATTERN.search(text)
    if label_match:
        value = _label_value(label_match.group(1))
        # guard against matching into unrelated sentences/navigation text
        if value and len(value) <= 30 and re.search(r"\d", value) and not re.search(r"[.!?]", value):
            return value

    heading = _first_heading(text)
    if not heading:
        return None
    for token in re.findall(r"[\[\(]([A-Za-z0-9\-\.]{2,15})[\]\)]", heading):
        if re.fullmatch(r"ab\d+", token, re.IGNORECASE):
            continue  # Abcam catalog number, not a clone code
        if not re.search(r"\d", token) or not re.search(r"[A-Za-z]", token):
            continue  # clone codes mix letters and digits (D3A7, EPR1639, 3G0)
        return token
    return None


def extract_immunogen(text: str):
    label_match = IMMUNOGEN_LABEL_PATTERN.search(text)
    if not label_match:
        return None
    value = _label_value(label_match.group(1))
    return value or None


def extract_target_protein(text: str):
    anti_match = ANTI_TARGET_PATTERN.search(text)
    if anti_match:
        return anti_match.group(1).strip()

    heading = _first_heading(text)
    if not heading:
        return None
    title = re.sub(r"\[[^\]]*\]", "", heading)
    # strip parens repeatedly to handle nesting, e.g. "(1461 (3G0))"
    prev = None
    while prev != title:
        prev = title
        title = re.sub(r"\([^()]*\)", "", title)
    title = title.replace("(", "").replace(")", "")
    title = re.sub(r"#\S+", "", title)
    title = TITLE_NOISE_WORDS.sub("", title)
    title = re.sub(r"\s+", " ", title).strip(" -,.")
    return title or None


def extract_fields(text: str) -> dict:
    rrid_match = RRID_PATTERN.search(text)
    rrid = f"RRID:{rrid_match.group(1)}" if rrid_match else None

    clonality = None
    if re.search(r"\bmonoclonal\b", text, re.IGNORECASE):
        clonality = "monoclonal"
    elif re.search(r"\bpolyclonal\b", text, re.IGNORECASE):
        clonality = "polyclonal"

    clone_id = extract_clone_id(text)
    host_match = HOST_PATTERN.search(text)

    return {
        "rrid": rrid,
        "rrid_missing": rrid is None,
        "clone_id": clone_id,
        "host_species": host_match.group(1).title() if host_match else None,
        "clonality": clonality,
        "target_protein": extract_target_protein(text),
        "immunogen": extract_immunogen(text),
        "lot_variable": clonality == "polyclonal" and clone_id is None,
    }


def resolve_identity(url: str, api_key: str) -> dict:
    raw_text = fetch_page(url, api_key)
    fields = extract_fields(raw_text)
    return {"url": url, **fields, "raw_text": raw_text}


TEST_URLS = [
    # Abcam
    "https://www.abcam.com/en-us/products/primary-antibodies/raidd-antibody-epr1639-ab76465",
    # Cell Signaling Technology
    "https://www.cellsignal.com/products/primary-antibodies/akt-pan-c67e7-rabbit-mab/4691",
    # Thermo Fisher
    "https://www.thermofisher.com/antibody/product/Adenovirus-Antibody-clone-3G0-Monoclonal/MA1-82981",
]


def main():
    load_dotenv()
    api_key = os.environ.get("LINKUP_API_KEY")
    if not api_key:
        print("ERROR: LINKUP_API_KEY not set (check .env)", file=sys.stderr)
        sys.exit(1)

    for url in TEST_URLS:
        try:
            result = resolve_identity(url, api_key)
        except urllib.error.HTTPError as e:
            result = {"url": url, "error": f"HTTP {e.code}: {e.reason}"}
        except urllib.error.URLError as e:
            result = {"url": url, "error": str(e.reason)}
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
