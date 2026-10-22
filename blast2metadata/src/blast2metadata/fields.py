"""BLAST tabular field dictionary, column labels and -outfmt handling."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple


BLAST_HEADERS: List[str] = [
    "qseqid", "sseqid", "stitle", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs",
]

NT_META_COLS: List[str] = [
    "accession", "version", "locus", "length", "definition",
    "molecule_type", "topology", "division", "update_date", "create_date",
    "organism", "taxonomy",
    "reference_titles", "reference_authors", "reference_journals", "comment",
]

NUMERIC_FLOAT_COLS = {"pident", "evalue", "qcovs", "bitscore", "ppos", "qcovhsp", "qcovus", "pident_x_qcovs"}
NUMERIC_INT_COLS   = {"length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
                      "qlen", "slen", "nident", "positive", "gaps", "score", "qframe", "sframe"}

# 修改 (2026-10-11): 表头功能。BLAST 有 50 个列代码，用户记不住各格式的列含义；
# 这里内置官方列表，用于 --header_style 和 --mode add-header。
# BLAST+ tabular field specifiers (outfmt 6/7/10): code -> (outfmt 7 '# Fields' label, description).
# Generated from `blastn -help` and a real outfmt 7 run, BLAST 2.17.0+ (2026-10-11).
BLAST_FIELDS: Dict[str, Tuple[str, str]] = {
    'qseqid': ('query id', 'Query Seq-id'),
    'qgi': ('query gi', 'Query GI'),
    'qacc': ('query acc.', 'Query accession'),
    'qaccver': ('query acc.ver', 'Query accession.version'),
    'qlen': ('query length', 'Query sequence length'),
    'sseqid': ('subject id', 'Subject Seq-id'),
    'sallseqid': ('subject ids', "All subject Seq-id(s), separated by a ';'"),
    'sgi': ('subject gi', 'Subject GI'),
    'sallgi': ('subject gis', 'All subject GIs'),
    'sacc': ('subject acc.', 'Subject accession'),
    'saccver': ('subject acc.ver', 'Subject accession.version'),
    'sallacc': ('subject accs.', 'All subject accessions'),
    'slen': ('subject length', 'Subject sequence length'),
    'qstart': ('q. start', 'Start of alignment in query'),
    'qend': ('q. end', 'End of alignment in query'),
    'sstart': ('s. start', 'Start of alignment in subject'),
    'send': ('s. end', 'End of alignment in subject'),
    'qseq': ('query seq', 'Aligned part of query sequence'),
    'sseq': ('subject seq', 'Aligned part of subject sequence'),
    'evalue': ('evalue', 'Expect value'),
    'bitscore': ('bit score', 'Bit score'),
    'score': ('score', 'Raw score'),
    'length': ('alignment length', 'Alignment length'),
    'pident': ('% identity', 'Percentage of identical matches'),
    'nident': ('identical', 'Number of identical matches'),
    'mismatch': ('mismatches', 'Number of mismatches'),
    'positive': ('positives', 'Number of positive-scoring matches'),
    'gapopen': ('gap opens', 'Number of gap openings'),
    'gaps': ('gaps', 'Total number of gaps'),
    'ppos': ('% positives', 'Percentage of positive-scoring matches'),
    'frames': ('query/sbjct frames', "Query and subject frames separated by a '/'"),
    'qframe': ('query frame', 'Query frame'),
    'sframe': ('sbjct frame', 'Subject frame'),
    'btop': ('BTOP', 'Blast traceback operations (BTOP)'),
    'staxid': ('subject tax id', 'Subject Taxonomy ID'),
    'ssciname': ('subject sci name', 'Subject Scientific Name'),
    'scomname': ('subject com names', 'Subject Common Name'),
    'sblastname': ('subject blast name', 'Subject Blast Name'),
    'sskingdom': ('subject super kingdom', 'Subject Super Kingdom'),
    'staxids': ('subject tax ids', "unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)"),
    'sscinames': ('subject sci names', "unique Subject Scientific Name(s), separated by a ';'"),
    'scomnames': ('subject com names', "unique Subject Common Name(s), separated by a ';'"),
    'sblastnames': ('subject blast names', "unique Subject Blast Name(s), separated by a ';' (in alphabetical order)"),
    'sskingdoms': ('subject super kingdoms', "unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order)"),
    'stitle': ('subject title', 'Subject Title'),
    'salltitles': ('subject titles', "All Subject Title(s), separated by a '<>'"),
    'sstrand': ('subject strand', 'Subject Strand'),
    'qcovs': ('% query coverage per subject', 'Query Coverage Per Subject'),
    'qcovhsp': ('% query coverage per hsp', 'Query Coverage Per HSP'),
    'qcovus': ('% query coverage per uniq subject', 'Query Coverage Per Unique Subject (blastn only)'),
}
BLAST_STD_FIELDS: List[str] = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()

_BLAST_LABEL_TO_CODE: Dict[str, str] = {lab: code for code, (lab, _d) in BLAST_FIELDS.items()}

# Non-BLAST columns this script adds (Entrez / GenBank metadata).
OUTPUT_FIELD_DESC: Dict[str, str] = {
    "pident_x_qcovs": "pident x qcovs / 100 (ranking heuristic; not whole-query identity)",
    "accession": "GenBank accession",
    "version": "GenBank accession.version",
    "locus": "GenBank LOCUS name",
    "metadata_length": "Subject sequence length in GenBank",
    "definition": "GenBank DEFINITION line",
    "molecule_type": "Molecule type",
    "topology": "Topology (linear/circular)",
    "division": "GenBank division",
    "update_date": "GenBank update date",
    "create_date": "GenBank create date",
    "organism": "Organism",
    "taxonomy": "Taxonomic lineage",
    "reference_titles": "Reference titles",
    "reference_authors": "Reference authors",
    "reference_journals": "Reference journals",
    "comment": "GenBank COMMENT",
}


def column_label(code: str, style: str, metadata_table: bool = False) -> str:
    """Header text for one column: style code | label | both.

    In metadata-only tables 'length' is the GenBank sequence length, not the alignment length.
    """
    if metadata_table and code == "length":
        desc = OUTPUT_FIELD_DESC["metadata_length"]
    else:
        desc = BLAST_FIELDS[code][1] if code in BLAST_FIELDS else OUTPUT_FIELD_DESC.get(code, "")
    if style == "code" or not desc:
        return code
    return desc if style == "label" else f"{code} ({desc})"


def expand_outfmt_fields(tokens: Sequence[str]) -> List[str]:
    """Expand the BLAST keyword 'std' and drop 'delim=' tokens."""
    out: List[str] = []
    for t in tokens:
        if t == "std":
            out.extend(BLAST_STD_FIELDS)
        elif not t.startswith("delim="):
            out.append(t)
    return out


def columns_from_outfmt(spec: str) -> Optional[List[str]]:
    """'6 qseqid sseqid ...' / '10 std staxids' -> column codes; None if no explicit column list."""
    toks = (spec or "").split()
    if not toks or toks[0] not in {"6", "7", "10", "20"}:
        return None
    fields = expand_outfmt_fields(toks[1:])
    if not fields:
        return None
    unknown = [f for f in fields if f not in BLAST_FIELDS]
    if unknown:
        raise SystemExit(f"ERROR: unknown BLAST outfmt field(s): {unknown}. Known: {', '.join(BLAST_FIELDS)}")
    return fields

# Codex 修改: 增加集中式扩展名、模式和 accession DB 识别常量。
# 理由: 修复 metadata-only/topN-only 输入发现不完整，以及 Entrez 蛋白/核酸库误判问题。
FASTA_EXTENSIONS_FALLBACK = [".fasta", ".fas", ".faa", ".fa", ".fna"]
ACCESSION_LIST_EXTENSIONS = [".txt", ".tsv", ".csv", ".list", ".acc", ".accessions"]
BLAST_TABLE_EXTENSIONS    = [".tsv", ".tab", ".txt", ".csv"]
BLAST_REQUIRED_MODES      = {"blast-only", "blast2topN", "blast2meta", "blast2topN2metadata"}

PROTEIN_ACCESSION_PREFIXES = (
    "NP_", "YP_", "WP_", "XP_", "AP_", "ZP_", "CAA", "AAB", "AAC", "AAD",
    "AAH", "AAM", "AAN", "AAO", "AAP", "AAQ", "AAR", "AAS", "AAT",
)
NUCLEOTIDE_ACCESSION_PREFIXES = (
    "NC_", "NG_", "NM_", "NR_", "NT_", "NW_", "NZ_", "XM_", "XR_", "AC_",
    "AP", "CP", "CM", "JX", "KC", "KF", "KP", "KU", "KY", "MK", "MN",
    "MT", "MW", "MZ", "OP", "OQ", "OR", "PQ",
)

# ─────────────────────────────────────────────
# Fix-3: build_outfmt_string → returns (str, List[str])
# ─────────────────────────────────────────────
def build_outfmt_string(outfmt_value: str) -> Tuple[str, List[str]]:
    """
    Returns (blast_outfmt_arg, list_of_column_names).

    - "6" / "7" / ""  → use all BLAST_HEADERS
    - "6 col1 col2 …" → validate required headers present, use provided list
    - anything else   → pass through unchanged, assume BLAST_HEADERS order
    """
    s = (outfmt_value or "").strip()

    # Fix-3: default empty / bare digit to standard 14-column format
    if s in {"", "6", "7"}:
        fmt_num = s if s in {"6", "7"} else "6"
        fmt_str = fmt_num + " " + " ".join(BLAST_HEADERS)
        return fmt_str, list(BLAST_HEADERS)

    toks = s.split()
    if toks and toks[0] in {"6", "7"}:
        fields = expand_outfmt_fields(toks[1:])
        s = " ".join([toks[0], *fields])
        if not fields:
            # "6" or "7" with no columns → use defaults
            fmt_str = toks[0] + " " + " ".join(BLAST_HEADERS)
            return fmt_str, list(BLAST_HEADERS)
        missing = [h for h in BLAST_HEADERS if h not in fields]
        if missing:
            raise SystemExit(
                f"ERROR: --outfmt_blast '{outfmt_value}' is missing required fields: "
                f"{missing}. Use --outfmt_blast 6 to get all required columns automatically."
            )
        return s, list(fields)

    # Codex 修改: 非 6/7 的 BLAST 输出格式直接报错。
    # 理由: 这个 pipeline 后续解析依赖 tabular 输出；原代码会放行 XML/pairwise，
    #       然后在 TopN/metadata 阶段静默解析为空或错列。
    # 原代码:
    #   # Non-tabular format (0, 5, etc.) – pass through, assume BLAST_HEADERS for parsing
    #   return s, list(BLAST_HEADERS)
    raise SystemExit(
        f"ERROR: --outfmt_blast '{outfmt_value}' is not supported by this pipeline. "
        "Use tabular BLAST output format 6 or 7."
    )

_INT_FIELDS = {"qlen", "slen", "qstart", "qend", "sstart", "send", "length", "nident",
               "mismatch", "positive", "gapopen", "gaps", "score", "qframe", "sframe"}
_FLOAT_FIELDS = {"evalue", "bitscore", "pident", "ppos", "qcovs", "qcovhsp", "qcovus"}
ADD_HEADER_EXTENSIONS = BLAST_TABLE_EXTENSIONS + [".out", ".m8", ".blast", ".outfmt6", ".outfmt7", ".blastn", ".blastp"]
