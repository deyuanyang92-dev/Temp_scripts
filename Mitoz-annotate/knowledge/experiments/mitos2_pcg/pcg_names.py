"""Map a GenBank CDS feature to one of the 13 metazoan PCG names (shared by the scripts in this directory)."""
import re

PATTERNS = [
    ("atp8", r"\bATP\s*8\b|ATP8|ATPASE\s*(SUBUNIT\s*)?8\b|SYNTHASE.*SUBUNIT\s*8\b"),
    ("atp6", r"\bATP\s*6\b|ATP6|ATPASE\s*(SUBUNIT\s*)?6\b|SYNTHASE.*SUBUNIT\s*6\b"),
    ("nad4l", r"\bN(AD|D)\s*4L\b|SUBUNIT\s*4L\b"),
    ("nad4", r"\bN(AD|D)\s*4\b|SUBUNIT\s*4\b"),
    ("nad1", r"\bN(AD|D)\s*1\b|SUBUNIT\s*1\b"),
    ("nad2", r"\bN(AD|D)\s*2\b|SUBUNIT\s*2\b"),
    ("nad3", r"\bN(AD|D)\s*3\b|SUBUNIT\s*3\b"),
    ("nad5", r"\bN(AD|D)\s*5\b|SUBUNIT\s*5\b"),
    ("nad6", r"\bN(AD|D)\s*6\b|SUBUNIT\s*6\b"),
    ("cob", r"\bCYT\s*B\b|CYTB|\bCOB\b|CYTOCHROME\s*B\b"),
    ("cox1", r"\bCO(X|I)\s*(1|I)\b|COX1|COI\b|CYTOCHROME\s*C\s*OXIDASE\s*SUBUNIT\s*(1|I)\b"),
    ("cox2", r"\bCO(X|I)\s*(2|II)\b|COX2|COII\b|CYTOCHROME\s*C\s*OXIDASE\s*SUBUNIT\s*(2|II)\b"),
    ("cox3", r"\bCO(X|I)\s*(3|III)\b|COX3|COIII\b|CYTOCHROME\s*C\s*OXIDASE\s*SUBUNIT\s*(3|III)\b"),
]


def pcg_name(ft):
    """Gene name of a CDS feature, from /gene first, then /product; None if not one of the 13."""
    for qual in ("gene", "product"):
        for text in ft.qualifiers.get(qual, []):
            up = text.upper().replace("-", " ").replace("_", " ")
            if "NADH" in up or "NAD" in up or re.search(r"\bND\d", up):
                order = ["nad4l", "nad1", "nad2", "nad3", "nad4", "nad5", "nad6"]
            elif "ATP" in up:
                order = ["atp8", "atp6"]
            else:
                order = ["cob", "cox1", "cox2", "cox3"]
            for name, pat in PATTERNS:
                if name in order and re.search(pat, up):
                    return name
    return None
