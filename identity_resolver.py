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

SPECIES_WORDS = (
    r"(Rabbit|Mouse|Rat|Goat|Sheep|Chicken|Guinea\s?Pig|Hamster|Donkey|Human|Bovine)"
)

# (?<![\w/-]) keeps labels from matching inside URL slugs, e.g. a tracking-pixel
# URL containing "...Antibody-clone-3G0-Monoclonal%2F..." would otherwise match
# "Clone" and capture garbage from the URL instead of the real field.
LABEL_GUARD = r"(?<![\w/-])"
RRID_LABEL_PATTERN = re.compile(LABEL_GUARD + r"RRID\b", re.IGNORECASE)
CLONE_LABEL_PATTERN = re.compile(
    LABEL_GUARD + r"Clone(?:\s*(?:ID|No\.?|Number))?\b", re.IGNORECASE
)
HOST_LABEL_PATTERN = re.compile(
    LABEL_GUARD + r"(?:Host|Source)(?:\s*/\s*Isotype|\s+Species)?\b",
    re.IGNORECASE,
)
IMMUNOGEN_LABEL_PATTERN = re.compile(LABEL_GUARD + r"Immunogen\b", re.IGNORECASE)
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


def _value_after(text: str, label_match, max_span: int = 80):
    """Given a match on a field label, pull the value that follows -- whether
    the source used 'Label: value', 'Label value', a markdown table row
    'Label | value |', or the label and value sit on separate lines
    ('Label\\n\\nValue\\n\\nNextLabel...', common in rendered spec tables)."""
    window = text[label_match.end(): label_match.end() + max_span]
    window = re.sub(r"^[\s:#\|]+", "", window)
    value = re.split(r"\n\s*\n|\|", window, maxsplit=1)[0]
    return value.strip()


def extract_rrid(text: str):
    label_match = RRID_LABEL_PATTERN.search(text)
    if not label_match:
        return None
    value = _value_after(text, label_match, max_span=20)
    rrid_match = re.search(r"AB_\d+", value)
    return rrid_match.group(0) if rrid_match else None


def extract_clone_id(text: str):
    label_match = CLONE_LABEL_PATTERN.search(text)
    if label_match:
        value = _value_after(text, label_match)
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


def extract_host_species(text: str):
    label_match = HOST_LABEL_PATTERN.search(text)
    if not label_match:
        return None
    value = _value_after(text, label_match, max_span=40)
    species_match = re.search(SPECIES_WORDS, value, re.IGNORECASE)
    return species_match.group(1).title() if species_match else None


def extract_immunogen(text: str):
    label_match = IMMUNOGEN_LABEL_PATTERN.search(text)
    if not label_match:
        return None
    value = _value_after(text, label_match, max_span=200)
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
    # Linkup's markdown escapes underscores (e.g. "AB\_926850"); un-escape
    # before matching so patterns don't have to account for the backslash.
    text = text.replace("\\_", "_")

    rrid_id = extract_rrid(text)
    rrid = f"RRID:{rrid_id}" if rrid_id else None

    clonality = None
    if re.search(r"\bmonoclonal\b", text, re.IGNORECASE):
        clonality = "monoclonal"
    elif re.search(r"\bpolyclonal\b", text, re.IGNORECASE):
        clonality = "polyclonal"

    clone_id = extract_clone_id(text)

    return {
        "rrid": rrid,
        "rrid_missing": rrid is None,
        "clone_id": clone_id,
        "host_species": extract_host_species(text),
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
        except (OSError, json.JSONDecodeError) as e:
            result = {"url": url, "error": str(e)}
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
