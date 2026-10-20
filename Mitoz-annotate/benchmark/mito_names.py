"""Gene-name normalisation shared by the benchmark scripts (PCGs, rRNAs, tRNAs)."""
import re

ALIAS = {"coi": "cox1", "coii": "cox2", "coiii": "cox3", "co1": "cox1", "co2": "cox2", "co3": "cox3",
         "cytb": "cob", "cyb": "cob", "nd1": "nad1", "nd2": "nad2", "nd3": "nad3", "nd4": "nad4", "nd4l": "nad4l",
         "nd5": "nad5", "nd6": "nad6", "atpase6": "atp6", "atpase8": "atp8",
         "rrns": "rrnS", "srrna": "rrnS", "12s": "rrnS", "12srrna": "rrnS", "rns": "rrnS", "rrn12": "rrnS",
         "rrnl": "rrnL", "lrrna": "rrnL", "16s": "rrnL", "16srrna": "rrnL", "rnl": "rrnL", "rrn16": "rrnL"}
PRODUCT = {"12s ribosomal rna": "rrnS", "small subunit ribosomal rna": "rrnS", "s-rrna": "rrnS",
           "16s ribosomal rna": "rrnL", "large subunit ribosomal rna": "rrnL", "l-rrna": "rrnL"}
PCG = {"cox1", "cox2", "cox3", "cob", "atp6", "atp8", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6"}
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}
ROMAN = {"i": "1", "ii": "2", "iii": "3"}
PCG_PRODUCT = [(r"cytochrome c oxidase subunit (i{1,3}|[123])\b", "cox"), (r"cytochrome b\b", "cob"),
               (r"nadh dehydrogenase subunit (4l|[1-6])\b", "nad"), (r"atp synthase (?:f0 )?subunit (6|8)\b", "atp")]


def gene_name(qualifiers: dict, ftype: str) -> str:
    """Normalised PCG / rRNA name (cox1, nad4l, rrnS ...); '' when not a standard one.
    Split genes such as COX1a / COX1b map to their gene (cox1)."""
    g = re.sub(r"[^a-z0-9]", "", (qualifiers.get("gene", [""])[0] or "").lower())
    g = ALIAS.get(g, g)
    m = re.fullmatch(r"(cox[123]|nad[1-6]|nad4l|atp[68]|cob)[abc]", g)
    if m:
        g = m.group(1)
    prod = (qualifiers.get("product", [""])[0] or "").lower()
    if ftype == "rRNA" and g not in ("rrnS", "rrnL"):
        g = PRODUCT.get(prod, g)
    if ftype == "CDS" and g not in PCG:
        for pat, base in PCG_PRODUCT:
            m = re.search(pat, prod)
            if m:
                g = base + ROMAN.get(m.group(1), m.group(1)) if base != "cob" else "cob"
                break
    return g if g in PCG or g in ("rrnS", "rrnL") else ""


def trna_name(qualifiers: dict) -> tuple:
    """(amino acid one-letter, anticodon or '') of a tRNA feature; ('', '') if unknown."""
    prod = qualifiers.get("product", [""])[0] or ""
    gene = qualifiers.get("gene", [""])[0] or ""
    m = re.search(r"tRNA-([A-Z][a-z]{2})", prod) or re.search(r"tRNA-([A-Z][a-z]{2})", gene)
    aa = AA3.get(m.group(1), "") if m else ""
    if not aa:
        m = re.match(r"trn([A-Z])", gene)
        aa = m.group(1) if m else ""
    text = " ".join(qualifiers.get("anticodon", []))
    ac = re.search(r"seq:([acgtu]{3})", text, re.I) or re.fullmatch(r"([acgtu]{3})", text.strip(), re.I)
    if not ac:
        ac = re.search(r"\(([acgtu]{3})\)", gene + " " + prod, re.I) or \
            re.search(r"anticodon[:= ]+([acgtu]{3})", " ".join(qualifiers.get("note", [])), re.I)
    return aa, (ac.group(1).lower().replace("u", "t") if ac else "")
