#!/usr/bin/env python3
"""mito_intron.py - intron-aware re-annotation of mitochondrial protein-coding genes.

MitoZ cannot annotate introns. This script adds them to MitoZ GenBank output with
a Python port of the intron workflow of MFannot (Beck & Lang,
https://github.com/BFL-lab/Mfannot, GPL-2.0; data https://github.com/BFL-lab/MFannot_data,
GPL-3.0; Lang et al. 2023, Front. Plant Sci., doi:10.3389/fpls.2023.1222186).
This file is therefore distributed under GPL-3.0.

Ported steps (MFannot line numbers refer to BFL-lab/Mfannot/mfannot):
  0. Library: MFannot's curated protein collections (models/mfannot_proteins, from
     MFannot_data/protein_collections/identified); gene name from the header
     "<species> mt <gene>;" exactly as mfannot:1180 parses it.
  1. Pre-screen (RunFlip / BlastFlipVSGene / FillFLIPBLASTPROTArrayFromBlastResult /
     SelectProteinForExo, mfannot:962-1480): ORFs >= 40 aa starting at the first
     start codon (flip -m -l 40 -g <code>), blastp vs the library (BLOSUM62); hits
     kept at E <= 1e-10 (and E <= 1e-11 when the top HSP is < 50 % identical); an
     ORF is assigned to every gene whose hit covers a non-overlapping part of it;
     a gene with >= 2 ORFs on one strand separated by >= 142 nt (or ORFs on both
     strands) is a candidate intron-containing gene.
  2. Gene model (FindExonsInHypProtArray, SortAndGroupC4Report, DefineBestSol,
     mfannot:2455-2960): Exonerate protein2genome with the best-hit protein and
     MFannot's parameters (--minintron 142 --maxintron 5000 --intronpenalty -35,
     neutral splice PSSM, genetic code); all reports with raw score >= 200 are
     chained (query overlap < 5 % of protein length, mutually closest in the
     genome) and the group covering most of the protein is used - this joins
     exons around introns longer than --maxintron.
  3. Intron type: Infernal cmsearch with Rfam RF00028 (group I), RF00029 + the seven
     group II D1-D4 models (RF01998, RF01999, RF02001, RF02003-05, RF02012), and the 14
     group I subtype models of Selles Vidal et al. 2025 (doi:10.1093/bioadv/vbaf020;
     downloaded once into ~/.cache), replacing MFannot's RNAweasel/ERPIN step.
  4. Boundary adjustment (Adjust_all_intronic_junctions, mfannot:5637-6100): both
     intron ends shifted by the same k in -6..6; group I scores exon-terminal U and
     intron-terminal omega-G (adjust when >= 19), group II GUGYG...AY (>= 40).

Documented deviations: step 3 uses Rfam instead of RNAweasel; in step 4 each
shift is judged against the homolog protein (MFannot compares with the starting
model, which Exonerate can get wrong at the junction - Metridium senile COX1:
model HMFTSG, NCBI and homolog HMFTVG - and would then reject the true boundary);
3' ends use the incomplete-stop (T / TA) convention of animal mtDNA; trans-spliced
genes (exons on both strands) are reported but not modelled; genes spanning the
origin of a circular record are not handled.

Usage:
  python3 mito_intron.py -i final_gb_dir -o out_dir --genetic_code 5 [--bin env/bin]
Outputs: out_dir/<same relative paths> (changed records rewritten, others copied byte
for byte) and out_dir/intron_report.tsv.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation,
                            SeqFeature)

HERE = Path(__file__).resolve().parent
DEFAULT_LIB = HERE / "models" / "mfannot_proteins"
DEFAULT_CM = HERE / "models" / "rfam_group_I_II_introns.cm"
SPLICE_PSSM = HERE / "models" / "neutral.pssm"      # MFannot_data/models/Splice/neutral.pssm
GB_EXT = {".gb", ".gbk", ".gbf", ".genbank", ".gbff"}
ADJUST_RANGE = range(-6, 7)          # mfannot:5640 %adjustement, same k at both ends
TYPE_I_MIN, TYPE_II_MIN = 19, 40     # mfannot:5735 / Adjust_TypeII
# MFannot defaults (PirModels/Option.pir, mfannot:284)
MIN_ORF_AA, BLAST_EVALUE, SHORT_HSP_EVALUE = 40, 1e-10, 1e-11
MIN_INTRON, MAX_INTRON, MIN_REPORT_SCORE = 142, 5000, 200
HEADER_RE = re.compile(r"(pt|pt_cyano|mt|mt_alpha)\s*([\S]*)\s*;", re.I)   # mfannot:1180


@dataclass
class LibProt:
    species: str
    gene: str
    seq: str


@dataclass
class Orf:
    strand: int
    start: int                        # transcript coords (1-based) on its strand
    end: int
    evalue: float = 1.0
    prot_id: str = ""
    q_min: int = 0                    # ORF aa range covered by the hit's HSPs
    q_max: int = 0


@dataclass
class Model:
    query: str
    gene: str
    strand: int                       # +1 / -1
    exons: List[Tuple[int, int]]      # 1-based inclusive, ascending genome coords
    score: int
    qcov: float
    qmiss_c: int = 0                  # homolog residues left unaligned at the C-terminus
    qmiss_n: int = 0                  # ... and at the N-terminus
    introns: List[dict] = field(default_factory=list)

    @property
    def span(self) -> Tuple[int, int]:
        return self.exons[0][0], self.exons[-1][1]


TIMEOUT = 900                          # seconds per external program call (--timeout)
CM_CPU = 1                             # cmsearch threads (--threads)


def run_tool(cmd: List[str], what: str) -> str:
    """Run an external program; on timeout or failure return '' and warn, so one hard
    sequence cannot stall or kill a batch."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=TIMEOUT).stdout
    except subprocess.TimeoutExpired:
        print(f"[mito_intron] {what} timed out after {TIMEOUT}s - skipped", file=sys.stderr)
    except OSError as exc:
        print(f"[mito_intron] {what} failed: {exc}", file=sys.stderr)
    return ""


# ── step 0: library ──────────────────────────────────────────────────────────
def iter_files(paths: List[str], exts: set) -> List[Path]:
    out: List[Path] = []
    for p in map(Path, paths):
        if p.is_dir():
            out += sorted(x for x in p.rglob("*") if x.suffix.lower() in exts and x.is_file())
        elif p.is_file():
            out.append(p)
    return out


def norm_species(s: str) -> str:
    return " ".join(s.replace("_", " ").split()).lower()


def load_library(paths: List[str]) -> Dict[str, LibProt]:
    """MFannot-format protein FASTA ('>Homo sapiens mt cox1; 513 aa')."""
    lib: Dict[str, LibProt] = {}
    for f in iter_files(paths, {".faa", ".fa", ".fasta", ".pep"}):
        for r in SeqIO.parse(str(f), "fasta"):
            m = HEADER_RE.search(r.description)
            if not m or not m.group(2):
                continue
            lib[f"L{len(lib):05d}"] = LibProt(norm_species(r.description[:m.start()]), m.group(2),
                                              str(r.seq).upper().rstrip("*"))
    return lib


def load_ref_gb(paths: List[str], gencode: int) -> Dict[str, LibProt]:
    """Extra homologs from close-relative GenBank records (MitoFinder's reference idea): CDS
    translations of the standard mitochondrial protein genes, gene names normalised to the
    MFannot library spelling so they group with it; CDS with internal stops are skipped."""
    canon = {p.stem.lower(): p.stem for p in DEFAULT_LIB.glob("*.faa")}
    by_product = {v: k for k, v in PRODUCTS.items()}
    lib: Dict[str, LibProt] = {}
    for f in iter_files(paths, GB_EXT):
        for rec in SeqIO.parse(str(f), "genbank"):
            sp = norm_species(rec.annotations.get("organism", ""))
            for ft in rec.features:
                if ft.type != "CDS":
                    continue
                g = norm_gene(ft.qualifiers.get("gene", [""])[0])
                if g not in PRODUCTS:
                    g = by_product.get(ft.qualifiers.get("product", [""])[0].strip().lower(), "")
                if g not in PRODUCTS:
                    continue
                prot = ft.qualifiers.get("translation", [""])[0]
                if not prot:
                    table = int(ft.qualifiers.get("transl_table", [gencode])[0])
                    nt = ft.location.extract(rec.seq)
                    prot = str(nt[: len(nt) - len(nt) % 3].translate(table=table))
                prot = prot.rstrip("*").upper()
                if "*" in prot or len(prot) < 50:
                    continue
                lib[f"R{len(lib):05d}"] = LibProt(sp, canon.get(g.lower(), g), prot)
    return lib


# ── step 1: pre-screen (flip + blastp + SelectProteinForExo) ─────────────────
def find_orfs(seq: str, gencode: int, min_aa: int = MIN_ORF_AA) -> List[Tuple[Orf, str]]:
    """flip -m -l <min_aa>: stop-to-stop segments, ORF from the first start codon."""
    table = CodonTable.unambiguous_dna_by_id[gencode]
    starts, stops = set(table.start_codons), set(table.stop_codons)
    out: List[Tuple[Orf, str]] = []
    for strand, s in ((1, seq), (-1, str(Seq(seq).reverse_complement()))):
        L = len(s)
        for frame in range(3):
            seg = frame
            for i in range(frame, L - 2, 3):
                last = i + 6 > L
                if s[i:i + 3] in stops or last:
                    end = i if s[i:i + 3] in stops else i + 3
                    j = next((x for x in range(seg, end, 3) if s[x:x + 3] in starts), None)
                    if j is not None and (end - j) // 3 >= min_aa:
                        prot = "M" + str(Seq(s[j + 3:end]).translate(table=gencode))
                        out.append((Orf(strand, j + 1, end), prot))
                    seg = i + 3
    return out


def prescreen(seq: str, lib: Dict[str, LibProt], lib_db: Path, gencode: int, tools: Dict[str, Optional[str]],
              tmp: Path, exclude_species: str = "", threads: int = 1) -> Dict[str, dict]:
    """Candidate intron-containing genes: {gene: {'prot': lib id, 'strand': 1|-1|2}}."""
    orfs = find_orfs(seq, gencode)
    if not orfs:
        return {}
    q = tmp / "orfs.faa"
    q.write_text("".join(f">O{i}\n{p}\n" for i, (_, p) in enumerate(orfs)))
    res = run_tool([str(tools["blastp"]), "-query", str(q), "-db", str(lib_db), "-matrix", "BLOSUM62",
                    "-outfmt", "6 qseqid sseqid evalue pident qstart qend", "-num_threads", str(threads)], "blastp")
    hits: Dict[str, Dict[str, List[Tuple[float, float, int, int]]]] = defaultdict(dict)
    for line in res.splitlines():
        qid, sid, ev_s, pid_s, qs, qe = line.split("\t")
        if exclude_species and lib[sid].species == exclude_species:
            continue
        hits[qid].setdefault(sid, []).append((float(ev_s), float(pid_s), int(qs), int(qe)))
    genes: Dict[str, List[Orf]] = defaultdict(list)
    for qid, by_subject in hits.items():
        base = orfs[int(qid[1:])][0]
        ranked = sorted(by_subject.items(), key=lambda kv: min(h[0] for h in kv[1]))
        accepted: List[Tuple[int, int]] = []
        for sid, hsps in ranked:
            ev, pid0 = min(h[0] for h in hsps), hsps[0][1]
            if ev > BLAST_EVALUE or (pid0 < 50 and ev > SHORT_HSP_EVALUE):
                break                                         # mfannot:1160 / 1204-1205
            good = [h for h in hsps if h[0] <= BLAST_EVALUE] or hsps[:1]
            fr = (min(h[2] for h in good), max(h[3] for h in good))
            if any(not (fr[1] < a or fr[0] > b) for a, b in accepted):
                continue                                      # overlapping frame (mfannot:1226-1233)
            accepted.append(fr)
            genes[lib[sid].gene].append(Orf(base.strand, base.start, base.end, ev, sid, fr[0], fr[1]))
    cands: Dict[str, dict] = {}
    for gene, os_ in genes.items():
        if len(os_) < 2:
            continue
        prot = min(os_, key=lambda o: o.evalue).prot_id
        strands = {o.strand for o in os_}
        if len(strands) == 2:
            cands[gene] = {"prot": prot, "strand": 2}
            continue
        segs = sorted((o.start + 3 * o.q_min, o.start + 3 * o.q_max) for o in os_)
        if any(abs(segs[i][1] - segs[i + 1][0]) >= MIN_INTRON for i in range(len(segs) - 1)):
            cands[gene] = {"prot": prot, "strand": os_[0].strand}       # mfannot:1447
    return cands


# ── step 2: Exonerate + report grouping ──────────────────────────────────────
@dataclass
class Report:
    qab: int
    qae: int
    tab: int
    tae: int
    score: int
    strand: int
    exons: List[Tuple[int, int]]


def run_exonerate(exe: str, prot_faa: Path, genome_fna: Path, gencode: int, maxintron: int) -> List[Report]:
    cmd = [exe, "--model", "protein2genome", "--showalignment", "no", "--showvulgar", "no",
           "--showtargetgff", "yes", "--minintron", str(MIN_INTRON), "--maxintron", str(maxintron),
           "--intronpenalty", "-35", "--geneticcode", str(gencode),
           *(["--splice3", str(SPLICE_PSSM), "--splice5", str(SPLICE_PSSM)] if SPLICE_PSSM.exists() else []),
           # %tS = target strand (%g, the splice-derived gene orientation, is "." without introns)
           "--ryo", ">>RYO\t%qab\t%qae\t%tab\t%tae\t%s\t%tS\n", str(prot_faa), str(genome_fna)]
    out = run_tool(cmd, "exonerate")
    reps: List[Report] = []
    exons: List[Tuple[int, int]] = []
    for line in out.splitlines():
        if line.startswith(">>RYO"):
            _, qab, qae, tab, tae, score, strand = line.split("\t")
            if exons:
                reps.append(Report(int(qab), int(qae), int(tab), int(tae), int(score),
                                   1 if strand == "+" else -1, sorted(exons)))
            exons = []
        elif "\texon\t" in line:
            c = line.split("\t")
            exons.append((int(c[3]), int(c[4])))
    return reps


def group_reports(reps: List[Report], qlen: int) -> List[List[Report]]:
    """SortAndGroupC4Report + DefineBestSol (mfannot): chain reports into gene models."""
    ov = 5 * qlen / 100
    groups = [[r] for r in reps if r.score >= MIN_REPORT_SCORE]

    def best_sol(i: int) -> Optional[int]:
        ps, pe, gs, ge = groups[i][0].qab, groups[i][-1].qae, groups[i][0].tab, groups[i][-1].tae
        best, dmin = None, None
        for j in range(len(groups) - 1, -1, -1):
            if j == i:
                continue
            ps2, pe2, gs2, ge2 = groups[j][0].qab, groups[j][-1].qae, groups[j][0].tab, groups[j][-1].tae
            if ps < ps2:
                if not (ps2 > pe or pe - ps2 < ov):
                    continue
                d = abs(ge - gs2)
            elif ps > ps2:
                if not (ps > pe2 or pe2 - ps < ov):
                    continue
                d = abs(ge2 - gs)
            else:
                continue
            if dmin is None or d < dmin:
                best, dmin = j, d
        return best

    i = 0
    while i < len(groups):
        j = best_sol(i)
        if j is None or best_sol(j) != i:
            i += 1
            continue
        a, b = groups[i], groups[j]
        groups[i] = a + b if a[0].qab < b[0].qab else b + a
        del groups[j]
        if j < i:
            i -= 1
    groups.sort(key=lambda g: (-sum(abs(r.qae - r.qab) + 1 for r in g), len(g)))
    return groups


def build_model(prot_id: str, gene: str, reps: List[Report], qlen: int) -> Tuple[Optional[Model], str]:
    """Chained reports -> gene model, or (None, reason)."""
    strands = {r.strand for r in reps}
    if len(strands) != 1:
        return None, "exons on both strands (trans-spliced?) - not modelled"
    exons: List[Tuple[int, int]] = []
    for e in sorted(e for r in reps for e in r.exons):
        if exons and e[0] - exons[-1][1] - 1 < MIN_INTRON:   # gap between chained reports, not an intron
            exons[-1] = (exons[-1][0], max(exons[-1][1], e[1]))
        else:
            exons.append(e)
    if len(exons) < 2:
        return None, ""
    coding = sum(b - a + 1 for a, b in exons)
    if coding > 3 * qlen + 33:            # MitoFinder geneChecker_fasta_gaps.py:150 length sanity check
        return None, (f"chained exons ({coding} nt) longer than the homolog allows ({3 * qlen + 33} nt) "
                      "- likely joins pieces of another gene; not modelled")
    covered = sum(abs(r.qae - r.qab) for r in reps)
    return Model(prot_id, gene, strands.pop(), exons, sum(r.score for r in reps), covered / max(qlen, 1),
                 qlen - max(r.qae for r in reps), min(r.qab for r in reps)), ""


# ── step 2: intron type with Infernal + Rfam ────────────────────────────────
SUBTYPE_CM_URL = "https://api.github.com/repos/LaraSellesVidal/Group1IntronDatabase/contents/subtype_CMs"
SUBTYPE_CM_DIR = Path.home() / ".cache" / "mito_intron" / "selles_vidal_2025_subtype_CMs"
SUBTYPE_REF = "Selles Vidal et al. 2025, doi:10.1093/bioadv/vbaf020"


def subtype_cms() -> List[Path]:
    """Group I subtype covariance models (IA1..IE3) of Selles Vidal et al. 2025, trained on 42,459
    group I introns. Not redistributed (no licence stated): downloaded once into a local cache;
    if the download fails the bundled Rfam models are used alone."""
    if not SUBTYPE_CM_DIR.exists() or not list(SUBTYPE_CM_DIR.glob("*.cm")):
        try:
            import json
            import urllib.request
            SUBTYPE_CM_DIR.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(SUBTYPE_CM_URL, timeout=30) as r:
                items = json.load(r)
            for it in items:
                if it["name"].endswith(".cm"):
                    urllib.request.urlretrieve(it["download_url"], SUBTYPE_CM_DIR / it["name"])
        except Exception as exc:                       # offline: carry on with Rfam only
            print(f"[mito_intron] group I subtype models not available ({exc}); using Rfam only",
                  file=sys.stderr)
    return sorted(SUBTYPE_CM_DIR.glob("*.cm"))


def intron_types(cmsearch: str, cms: List[Path], seqs: Dict[str, str], evalue: float, tmp: Path
                 ) -> Dict[str, Tuple[str, float, str]]:
    """Best covariance-model hit per intron: (group 'I'/'II', E-value, model name).
    Models: Rfam Intron_gpI (RF00028), Intron_gpII (RF00029, domains V-VI), the seven group II
    D1-D4 models (RF01998, RF01999, RF02001, RF02003-RF02005, RF02012), and the group I subtype
    models of Selles Vidal et al. 2025 when available."""
    cms = [c for c in cms if c.exists()]
    if not seqs or not cms:
        return {}
    cm = tmp / "all_models.cm"
    cm.write_text("".join(c.read_text() for c in cms))
    best: Dict[str, Tuple[str, float, str]] = {}

    def search(todo: Dict[str, str], fast: bool) -> None:
        fna, tbl = tmp / "introns.fna", tmp / "introns.tbl"
        fna.write_text("".join(f">{k}\n{v}\n" for k, v in todo.items()))
        tbl.unlink(missing_ok=True)
        run_tool([cmsearch, "--noali", "--toponly", *(["--rfam"] if fast else []), "--cpu", str(CM_CPU),
                  "--tblout", str(tbl), str(cm), str(fna)], "cmsearch")
        for line in tbl.read_text().splitlines() if tbl.exists() else []:
            if line.startswith("#"):
                continue
            c = line.split()
            name, model, strand, e = c[0], c[2], c[9], float(c[15])
            if strand != "+" or e > evalue:
                continue
            kind = "II" if (model.endswith("gpII") or model.startswith("group-II")) else "I"
            if name not in best or e < best[name][1]:
                best[name] = (kind, e, model)

    # two passes: Infernal's fast --rfam filters first; full-sensitivity search only for introns still
    # untyped (benchmark: same 147/165 NCBI introns typed, same group for all, ~2.4x faster)
    search(seqs, fast=True)
    rest = {k: v for k, v in seqs.items() if k not in best}
    if rest:
        search(rest, fast=False)
    return best


def type_label(kind: str, ev: float, model: str) -> str:
    if not kind:
        return "intron, type not identified"
    if re.fullmatch(r"I[A-E]\d?", model):
        return f"group I intron, subtype {model} ({SUBTYPE_REF} model, E={ev:.2g})"
    acc = {"Intron_gpI": "RF00028", "Intron_gpII": "RF00029"}.get(model, model)
    return f"group {kind} intron (Infernal/Rfam {acc}, E={ev:.2g})"


# ── step 3: boundary adjustment (port of mfannot Adjust_all_intronic_junctions) ──
def _ch(s: str, pos1: int) -> str:
    return s[pos1 - 1] if 1 <= pos1 <= len(s) else ""


def score_type_i(s: str, start: int, end: int, k: int) -> float:
    """ScoreNtLTypeI + ScoreNtRTypeI: exon ends ...GGT | intron ends with omega-G."""
    sc: float = 5 * (_ch(s, start + k - 1) == "T") + 3 * (_ch(s, start + k - 2) == "G") + 3 * (_ch(s, start + k - 3) == "G")
    sc += 10 * (_ch(s, end + k) == "G")
    sc += sum(0.5 for i in range(1, 9) if _ch(s, end + k - i) in ("A", "T"))
    return sc


def score_type_ii(s: str, start: int, end: int, k: int) -> float:
    """ScoreNtLTypeII + ScoreNtRTypeII - |k|: intron GTGYG ... AY."""
    sc: float = 10 * (_ch(s, start + k) == "G") + 7.5 * (_ch(s, start + k + 1) == "T") + 9 * (_ch(s, start + k + 2) == "G")
    sc += 10 * (_ch(s, start + k + 3) in ("C", "T")) + 10 * (_ch(s, start + k + 4) == "G")
    sc += 9 * (_ch(s, end + k) in ("C", "T")) + 9 * (_ch(s, end + k - 1) == "A")
    return sc - abs(k)


def _translate(s: str, exons: List[Tuple[int, int]], gencode: int) -> str:
    nt = "".join(s[a - 1:b] for a, b in exons)
    return str(Seq(nt[: len(nt) - len(nt) % 3]).translate(table=gencode))


def homolog_scorer(ref: str):
    """Local alignment score of a protein against the homolog that produced the model."""
    al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1)
    return lambda prot: al.score(prot, ref)


def _aligner(open_gap: float) -> PairwiseAligner:
    return PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=open_gap, extend_gap_score=-1)


def homolog_gap(prot: str, refs: List[str], junction_aa: int, window: int = 4) -> Tuple[float, float]:
    """Residues the homologs have but the spliced protein lacks at an exon junction
    (AnnotateMiniExonsByHMM idea: a missing mini-exon shows as a gap in our protein).
    A single pairwise alignment rarely opens a gap for a few residues, so a mild gap
    penalty is used and the signal is taken over many homologs of the same gene.
    Returns (median gap, fraction of homologs with a gap >= 1)."""
    al = _aligner(-6)
    gaps = []
    for ref in refs:
        aln = al.align(prot, ref)[0]
        blocks = list(zip(aln.aligned[0], aln.aligned[1]))
        g = max([int((rb[0] - ra[1]) - (pb[0] - pa[1])) for (pa, ra), (pb, rb) in zip(blocks, blocks[1:])
                 if abs(pa[1] - junction_aa) <= window] or [0])
        gaps.append(g)
    if not gaps:
        return 0.0, 0.0
    gaps.sort()
    return float(gaps[len(gaps) // 2]), sum(g >= 1 for g in gaps) / len(gaps)


def find_mini_exon(s: str, exons: List[Tuple[int, int]], i: int, gencode: int, refs: List[str],
                   k_max: int) -> Optional[Tuple[int, int, float]]:
    """Search intron i (transcript coords) for a short in-frame exon (1..k_max codons) whose
    insertion improves the alignment to every homolog in refs without adding stop codons;
    both flanking introns must stay >= MIN_INTRON. Returns (start, end, mean gain)."""
    al = _aligner(-11)
    start, end = exons[i][1] + 1, exons[i + 1][0] - 1
    base_p = _translate(s, exons, gencode)
    base = [al.score(base_p, r) for r in refs]
    stops0 = base_p.rstrip("*").count("*")
    # stage 1 (fast): one homolog, local alignment of a window around the junction only
    j = sum(b - a + 1 for a, b in exons[:i + 1]) // 3
    w0 = max(0, j - 50)
    base_w = al.score(base_p[w0:j + 50], refs[0])
    shortlist: List[Tuple[float, int, int]] = []
    for ln in range(3, 3 * k_max + 1, 3):
        for x in range(start + MIN_INTRON, end - MIN_INTRON - ln + 2):
            new = exons[:i + 1] + [(x, x + ln - 1)] + exons[i + 1:]
            prot = _translate(s, new, gencode)
            if prot.rstrip("*").count("*") > stops0:
                continue
            g1 = al.score(prot[w0:j + 50 + ln // 3], refs[0]) - base_w
            if g1 > 0:
                shortlist.append((g1, x, x + ln - 1))
    shortlist = sorted(shortlist, reverse=True)[:5]
    # stage 2 (exact): full proteins against every homolog
    best: Optional[Tuple[float, int, int]] = None
    for _, x, y in shortlist:
        prot = _translate(s, exons[:i + 1] + [(x, y)] + exons[i + 1:], gencode)
        gains = [al.score(prot, r) - b for r, b in zip(refs, base)]
        if min(gains) <= 0:
            continue                          # must help every homolog
        g = sum(gains) / len(gains)
        if best is None or g > best[0]:
            best = (g, x, y)
    return (best[1], best[2], best[0]) if best and best[0] >= 2 * (best[2] - best[1] + 1) / 3 else None


def add_mini_exons(s: str, exons: List[Tuple[int, int]], kinds: List[Any], notes: List[str],
                   scores: List[float], gencode: int, refs: List[str]):
    """Insert homolog-supported mini-exons; returns updated exons / kinds / notes / scores."""
    i = 0
    while i < len(exons) - 1 and refs:
        junction = sum(b - a + 1 for a, b in exons[:i + 1]) // 3
        med, frac = homolog_gap(_translate(s, exons, gencode), refs, junction)
        if med >= 1 and frac >= 0.5:
            mx = find_mini_exon(s, exons, i, gencode, refs[:3], int(med) + 4)
            if mx:
                exons = exons[:i + 1] + [(mx[0], mx[1])] + exons[i + 1:]
                tag = (f"mini-exon {mx[1] - mx[0] + 1} nt inserted (gap of ~{med:g} aa in {frac:.0%} of "
                       f"{len(refs)} homologs; MFannot AnnotateMiniExons idea)")
                kinds = kinds[:i] + [kinds[i], kinds[i]] + kinds[i + 1:]
                notes = notes[:i] + [f"{notes[i]}; {tag}"] * 2 + notes[i + 1:]
                scores = scores[:i] + [scores[i], scores[i]] + scores[i + 1:]
                i += 2
                continue
        i += 1
    return exons, kinds, notes, scores


def aa_score(s: str, exons: List[Tuple[int, int]], new: List[Tuple[int, int]], gencode: int,
             to_ref) -> Optional[float]:
    """CheckForCorrespondance, judged against the homolog protein.

    MFannot compares the shifted protein with the starting model and rejects any
    change scoring < 0 in PAM70. The starting model here comes from Exonerate's
    built-in splice model and can be wrong at the junction (Metridium senile COX1:
    model HMFTSG, NCBI and homolog HMFTVG), so the shifted protein is scored
    against the homolog instead: rejected when it aligns worse than the starting
    model, otherwise the score gain is returned.
    """
    old_p, new_p = _translate(s, exons, gencode), _translate(s, new, gencode)
    if old_p == new_p:
        return 9999.0
    gain = to_ref(new_p) - to_ref(old_p)
    return None if gain < 0 else gain


def adjust_introns(s: str, exons: List[Tuple[int, int]], kinds: List[str], gencode: int,
                   to_ref) -> Tuple[List[Tuple[int, int]], List[str], List[float]]:
    """s and exons in transcript orientation. Returns adjusted exons and a note per intron."""
    exons = list(exons)
    notes: List[str] = []
    scores: List[float] = []
    for i, kind in enumerate(kinds):
        start, end = exons[i][1] + 1, exons[i + 1][0] - 1
        prev_len, next_len = exons[i][1] - exons[i][0] + 1, exons[i + 1][1] - exons[i + 1][0] + 1
        scored = []
        for k in ADJUST_RANGE:
            if (prev_len + k <= 5 or next_len - k <= 5) and not (kind != "II" and k == 0):
                continue
            nt = score_type_ii(s, start, end, k) if kind == "II" else score_type_i(s, start, end, k)
            new = exons[:i] + [(exons[i][0], exons[i][1] + k), (exons[i + 1][0] + k, exons[i + 1][1])] + exons[i + 2:]
            aa = aa_score(s, exons, new, gencode, to_ref)
            if aa is not None:
                scored.append((nt, aa, abs(k), k, new))
        scored.sort(key=lambda x: (-(x[0] - x[2]), -x[1], x[2], x[3]))
        if not scored:
            notes.append("splice boundaries uncertain")
            scores.append(0.0)
            continue
        nt, _, _, k, new = scored[0]
        scores.append(nt)
        threshold = TYPE_II_MIN if kind == "II" else TYPE_I_MIN
        if nt >= threshold or kind == "I":      # mfannot: typed group I adjusted even below 19
            exons = new
            notes.append(f"boundaries adjusted by {k:+d} nt (splice score {nt:g})" if k else f"splice score {nt:g}")
        else:
            notes.append(f"splice boundaries uncertain (best score {nt:g})")
    return exons, notes, scores


def extend_to_stop(s: str, exons: List[Tuple[int, int]], gencode: int, limit: int, qmiss_c: int,
                   max_codons: int = 100, rule: str = "first") -> Tuple[List[Tuple[int, int]], bool]:
    """Set the 3' end of the CDS.

    Alignment reached the homolog's C-terminus (qmiss_c <= 1): the protein is complete,
    so end on the next codon if it is a stop, else on an incomplete stop (T / TA,
    completed by polyadenylation) as NCBI mitochondrial records do.
    Alignment stopped short: extend in frame to the first complete stop within max_codons
    (overlapping the next feature is allowed, as mitochondrial genes often overlap);
    without one, end on an incomplete stop just before the next feature (limit).
    """
    stops = CodonTable.unambiguous_dna_by_id[gencode].stop_codons
    a, b = exons[-1]
    if s[b - 3:b] in stops:
        return (exons), False
    if rule == "homolog" and qmiss_c <= 1:
        if s[b:b + 3] in stops and b + 3 <= limit:
            return (exons[:-1] + [(a, b + 3)]), False
        if s[b:b + 2] == "TA" and b + 2 <= limit:
            return (exons[:-1] + [(a, b + 2)]), False
        return (exons[:-1] + [(a, b + 1)] if s[b:b + 1] == "T" and b + 1 <= limit else exons), False
    # mitochondrial genes often overlap their neighbours, so a complete stop codon may lie
    # inside the next feature; only when none is found end on an incomplete stop before it
    n_codons = (len(s) - b) // 3 if rule == "first" else max(max_codons, qmiss_c + 5)   # "first": up to the ORF end
    for i in range(n_codons):
        c0 = b + 3 * i
        if c0 + 3 > len(s):
            if rule == "first":            # open frame runs off the 3' end: 3'-partial CDS
                return (exons[:-1] + [(a, c0)]), True
            break
        if s[c0:c0 + 3] in stops:
            return (exons[:-1] + [(a, c0 + 3)]), False
    for c0 in range(b, min(limit, len(s)), 3):
        partial = s[c0:limit]
        if len(partial) < 3:
            return (exons[:-1] + [(a, limit)] if partial in ("T", "TA") else exons), False
        if s[c0:c0 + 3] in stops:
            break
    return (exons), False


def mitoz_ends(rec, m: Model, exons_tx: List[Tuple[int, int]], L: int, window: int = 90
               ) -> Tuple[Optional[int], Optional[int]]:
    """5' start / 3' end (transcript coords) of existing CDS features (MitoZ's model) on the
    same strand that overlap the model and are in frame with its first / last exon."""
    a, b = exons_tx[0][0], exons_tx[-1][1]
    g_span = to_genome([(a, b)], m.strand, L)[0]
    start = end = None
    for ft in rec.features:
        if ft.type != "CDS" or ft.location.strand != m.strand:
            continue
        fs, fe = int(ft.location.start) + 1, int(ft.location.end)
        if not overlap((fs, fe), g_span):
            continue
        s_tx, e_tx = (fs, fe) if m.strand == 1 else (L + 1 - fe, L + 1 - fs)
        if abs(s_tx - a) <= window and (a - s_tx) % 3 == 0 and s_tx <= exons_tx[0][1]:
            start = s_tx if start is None or abs(s_tx - a) < abs(start - a) else start
        if abs(e_tx - b) <= window and e_tx >= exons_tx[-1][0]:
            end = e_tx if end is None or abs(e_tx - b) < abs(end - b) else end
    return start, end


def valid_exons(exons: List[Tuple[int, int]]) -> bool:
    """Exons non-empty, ordered and separated by at least MIN_INTRON nt."""
    return all(a <= b for a, b in exons) and all(exons[i][1] + MIN_INTRON <= exons[i + 1][0]
                                                for i in range(len(exons) - 1))


def extend_to_start(s: str, exons: List[Tuple[int, int]], gencode: int, lo: int, qmiss_n: int
                    ) -> Tuple[List[Tuple[int, int]], bool]:
    """Find the start codon. Used when the alignment stopped short of the homolog's N-terminus
    or the first codon is not a start codon: upstream in frame (no stop codon, not before the
    previous feature, lo) the start codon closest to where the homolog's first residue would be;
    if none, the nearest start codon a few codons downstream. If the open frame runs off the
    5' end of the sequence the CDS is 5'-partial (returns partial=True)."""
    table = CodonTable.unambiguous_dna_by_id[gencode]
    a, b = exons[0]
    if qmiss_n <= 1 and s[a - 1:a + 2] in table.start_codons:
        return exons, False
    target = a - 3 * max(qmiss_n, 0)
    best, ran_off = None, False
    for i in range(1, max(qmiss_n + 6, 30)):
        c0 = a - 1 - 3 * i
        if c0 < 0:
            ran_off = lo <= 1
            break
        if c0 < lo - 1 or s[c0:c0 + 3] in table.stop_codons:
            break
        if s[c0:c0 + 3] in table.start_codons and (best is None or abs(c0 + 1 - target) < abs(best - target)):
            best = c0 + 1
    if best is None and s[a - 1:a + 2] in table.start_codons:
        return exons, False
    if best is None and ran_off:
        c0 = (a - 1) % 3                          # first in-frame base of the sequence
        return [(c0 + 1, b)] + exons[1:], True
    if best is None:
        for i in range(1, 11):                    # downstream, inside the first exon
            c0 = a - 1 + 3 * i
            if c0 + 3 > b or s[c0:c0 + 3] in table.stop_codons:
                break
            if s[c0:c0 + 3] in table.start_codons:
                best = c0 + 1
                break
    return ([(best, b)] + exons[1:] if best else exons), False


def prev_feature_end(rec, m: Model, exons_tx: List[Tuple[int, int]], L: int) -> int:
    """Transcript position just after the nearest upstream annotated feature (other genes)."""
    start_tx = exons_tx[0][0]
    g_span = to_genome([(start_tx, exons_tx[-1][1])], m.strand, L)[0]
    lo = 1
    for ft in rec.features:
        if ft.type not in ("gene", "CDS", "tRNA", "rRNA"):
            continue
        fs, fe = int(ft.location.start) + 1, int(ft.location.end)
        if overlap((fs, fe), g_span) >= 30:
            continue
        end_tx = fe if m.strand == 1 else L + 1 - fs
        if end_tx < start_tx:
            lo = max(lo, end_tx + 1)
    return lo


def next_feature_start(rec, m: Model, exons_tx: List[Tuple[int, int]], L: int) -> int:
    """Transcript position just before the nearest downstream annotated feature (other genes)."""
    end_tx = exons_tx[-1][1]
    g_span = to_genome([(exons_tx[0][0], end_tx)], m.strand, L)[0]
    limit = L
    for ft in rec.features:
        if ft.type not in ("gene", "CDS", "tRNA", "rRNA"):
            continue
        fs, fe = int(ft.location.start) + 1, int(ft.location.end)
        if overlap((fs, fe), g_span) >= 30:          # the gene being re-modelled (or a nested ORF)
            continue
        start_tx = fs if m.strand == 1 else L + 1 - fe
        if start_tx > end_tx:
            limit = min(limit, start_tx - 1)
    return limit


# ── write-back ───────────────────────────────────────────────────────────────
def to_transcript(m: Model, L: int) -> List[Tuple[int, int]]:
    if m.strand == 1:
        return list(m.exons)
    return sorted((L + 1 - b, L + 1 - a) for a, b in m.exons)


def to_genome(exons: List[Tuple[int, int]], strand: int, L: int) -> List[Tuple[int, int]]:
    return list(exons) if strand == 1 else sorted((L + 1 - b, L + 1 - a) for a, b in exons)


def overlap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)


# Standard products of the mitochondrial protein genes (NCBI mitogenome submissions)
PRODUCTS = {
    "cox1": "cytochrome c oxidase subunit I", "cox2": "cytochrome c oxidase subunit II",
    "cox3": "cytochrome c oxidase subunit III", "cob": "cytochrome b", "cytb": "cytochrome b",
    "atp6": "ATP synthase F0 subunit 6", "atp8": "ATP synthase F0 subunit 8", "atp9": "ATP synthase F0 subunit 9",
    "nad1": "NADH dehydrogenase subunit 1", "nad2": "NADH dehydrogenase subunit 2",
    "nad3": "NADH dehydrogenase subunit 3", "nad4": "NADH dehydrogenase subunit 4",
    "nad4l": "NADH dehydrogenase subunit 4L", "nad5": "NADH dehydrogenase subunit 5",
    "nad6": "NADH dehydrogenase subunit 6",
}
ALIASES = {"coi": "cox1", "coii": "cox2", "coiii": "cox3", "cytb": "cob", "nd1": "nad1", "nd2": "nad2",
           "nd3": "nad3", "nd4": "nad4", "nd4l": "nad4l", "nd5": "nad5", "nd6": "nad6"}


def norm_gene(gene: str) -> str:
    g = gene.lower().replace("-", "").replace("_", "")
    return ALIASES.get(g, g)


def product_for(gene: str) -> str:
    return PRODUCTS.get(norm_gene(gene), f"{gene} protein")


def conflicting_genes(rec, gene: str, strand: int, exons: List[Tuple[int, int]]) -> List[str]:
    """Existing CDS (e.g. MitoZ's) of a *different* standard mitochondrial protein gene that the
    model's exons overlap by >= 30 nt (intron-encoded ORFs such as heg are not conflicts)."""
    out = []
    for f in rec.features:
        if f.type != "CDS" or f.location.strand != strand or not f.qualifiers.get("gene"):
            continue
        other = f.qualifiers["gene"][0]
        if norm_gene(other) == norm_gene(gene) or norm_gene(other) not in PRODUCTS:
            continue                      # same gene, or intron-encoded ORF (heg/LAGLIDADG/orf...) that may overlap exons
        parts = [(int(p.start) + 1, int(p.end)) for p in f.location.parts]
        if sum(overlap(x, y) for x in exons for y in parts) >= 30:
            out.append(other)
    return out


def repair_chain(rec, prot_id: str, gene: str, groups: List[List[Report]], qlen: int
                 ) -> Optional[Tuple[Model, List[Tuple[int, int]]]]:
    """Best chain minus its exons on another gene's CDS, chained with the next report that lies on the
    side the dropped exons were on. NC_006902.1 ND5: exonerate aligns nad5 residues 135-600 as an exon on
    ND4 (nad4/nad5 are homologous) + the true exon 2, and residues 1-239 as the true exon 1 in a separate
    report; their query overlap (> 5 %) keeps group_reports from chaining them. Returns (model, dropped exons)."""
    best = groups[0]
    strand = best[0].strand
    kept = [Report(r.qab, r.qae, r.tab, r.tae, r.score, r.strand,
                   [e for e in r.exons if not conflicting_genes(rec, gene, strand, [e])]) for r in best]
    lo = min(e[0] for r in best for e in r.exons)
    k_ex = [e for r in kept for e in r.exons]
    if not k_ex or len(k_ex) == sum(len(r.exons) for r in best):
        return None
    k_lo, k_hi = min(a for a, _ in k_ex), max(b for _, b in k_ex)
    dropped_low = lo < k_lo                            # dropped exon(s) before the kept ones in genome order
    dropped = [e for r in best for e in r.exons if e not in k_ex]
    for g in groups[1:6]:
        g_ex = [e for r in g for e in r.exons]
        if g[0].strand != strand or conflicting_genes(rec, gene, strand, g_ex):
            continue
        if (dropped_low and max(b for _, b in g_ex) < k_lo) or (not dropped_low and min(a for a, _ in g_ex) > k_hi):
            m, _ = build_model(prot_id, gene, g + [r for r in kept if r.exons], qlen)
            if m is not None:
                return m, dropped
    return None


def realign_masked(rec, seq: str, prot_id: str, gene: str, groups: List[List[Report]], qlen: int, pf: Path,
                   tmp: Path, a, tools: Dict[str, Optional[str]]) -> Optional[Model]:
    """repair_chain, then exonerate again with the dropped exons masked (N) and --maxintron long enough for
    the repaired model, so both splice sites of the new junction come from one spliced alignment (the
    repaired junction keeps the acceptor exonerate chose for the wrong donor: 8 nt off in NC_006902.1)."""
    rep = repair_chain(rec, prot_id, gene, groups, qlen)
    if rep is None:
        return None
    m, dropped = rep
    masked = list(seq)
    for x, y in dropped:
        masked[x - 1:y] = "N" * (y - x + 1)
    fna = tmp / "masked.fna"
    fna.write_text(">g\n" + "".join(masked) + "\n")
    longest = max(m.exons[i + 1][0] - m.exons[i][1] for i in range(len(m.exons) - 1))
    reps = run_exonerate(str(tools["exonerate"]), pf, fna, a.genetic_code, max(a.maxintron, longest + 200))
    for g in group_reports(reps, qlen)[:3]:
        alt, _ = build_model(prot_id, gene, g, qlen)
        if alt is not None and alt.strand == m.strand and not conflicting_genes(rec, gene, alt.strand, alt.exons) \
                and overlap(alt.span, m.span) > 0:
            return alt
    return m


def stop_qualifiers(rec, cds_loc, gencode: int) -> Dict[str, List[str]]:
    """Incomplete stop codon (T / TA completed by polyadenylation): NCBI wants
    /transl_except=(pos:..,aa:TERM) and the standard note."""
    nt = str(cds_loc.extract(rec.seq)).upper()
    tail = len(nt) % 3
    if tail == 0 or nt[-tail:] not in ("T", "TA"):
        return {}                         # complete codon, or not a polyadenylation-completed stop
    last = cds_loc.parts[-1]
    if last.strand == -1:
        a = int(last.start) + 1
        pos = f"complement({a}..{a + tail - 1})" if tail > 1 else f"complement({a})"
    else:
        b = int(last.end)
        pos = f"{b - tail + 1}..{b}" if tail > 1 else f"{b}"
    return {"transl_except": [f"(pos:{pos},aa:TERM)"],
            "note": ["TAA stop codon is completed by the addition of 3' A residues to the mRNA"]}


def apply_model(rec, m: Model, exons: List[Tuple[int, int]], intron_info: List[dict], gencode: int,
                partial5: bool = False, partial3: bool = False) -> str:
    """Replace MitoZ gene/CDS features that overlap the exons; add intron features. Returns gene name."""
    strand = m.strand
    hit, keep = [], []
    for ft in rec.features:
        if ft.type in ("gene", "CDS", "mRNA") and ft.location.strand == strand:
            s, e = int(ft.location.start) + 1, int(ft.location.end)
            if sum(overlap((s, e), x) for x in exons) >= 30:
                hit.append(ft)
                continue
        keep.append(ft)
    names = [f.qualifiers.get("gene", [""])[0] for f in hit if f.type == "CDS" and f.qualifiers.get("gene")]
    name = names[0] if names else m.gene
    products = [f.qualifiers["product"][0] for f in hit if f.type == "CDS" and f.qualifiers.get("product")]
    product = products[0] if products else product_for(name)
    # partial ends (gene runs off the sequence): 5' end is the low coordinate on +, the high one on -
    lo_open = partial5 if strand == 1 else partial3
    hi_open = partial3 if strand == 1 else partial5

    def pos(i: int, j: int, a: int, b: int):
        st = BeforePosition(a - 1) if (i == 0 and lo_open) else ExactPosition(a - 1)
        en = AfterPosition(b) if (j == len(exons) - 1 and hi_open) else ExactPosition(b)
        return FeatureLocation(st, en, strand=strand)
    locs = [pos(i, i, a, b) for i, (a, b) in enumerate(exons)]
    if strand == -1:
        locs = locs[::-1]
    cds_loc = CompoundLocation(locs) if len(locs) > 1 else locs[0]
    gene_loc = FeatureLocation(BeforePosition(exons[0][0] - 1) if lo_open else ExactPosition(exons[0][0] - 1),
                               AfterPosition(exons[-1][1]) if hi_open else ExactPosition(exons[-1][1]), strand=strand)
    note = "exon-intron structure: Exonerate protein2genome + MFannot boundary rules (mito_intron.py)"
    nt = cds_loc.extract(rec.seq)
    prot = str(nt[: len(nt) - len(nt) % 3].translate(table=gencode, to_stop=True))
    quals: Dict[str, List[str]] = {"gene": [name], "codon_start": ["1"], "transl_table": [str(gencode)],
                                   "product": [product], "note": [note]}
    for k, v in ({} if partial3 else stop_qualifiers(rec, cds_loc, gencode)).items():
        quals[k] = quals.get(k, []) + v
    quals["translation"] = [prot]
    new = [SeqFeature(gene_loc, type="gene", qualifiers={"gene": [name]}),
           SeqFeature(cds_loc, type="CDS", qualifiers=quals)]
    for n, info in enumerate(intron_info, 1):
        new.append(SeqFeature(FeatureLocation(info["start"] - 1, info["end"], strand=strand), type="intron",
                              qualifiers={"gene": [name], "number": [str(n)], "note": [info["note"]]}))
    rec.features = sorted(keep + new, key=lambda f: (int(f.location.start), f.type != "source", f.type != "gene"))
    return name


UNTYPED_MIN_SPLICE = 14     # untyped introns: type-I splice score needed to be written (benchmark: true >=14.5, false <=12.5)


def evidence(s: str, exons: List[Tuple[int, int]], kinds: List[str], splice: List[float], gencode: int,
             partial5: bool = False) -> Tuple[str, str]:
    """Confidence of a gene model, from evidence independent of existing annotations.
    high: every intron typed by Rfam and the spliced CDS has no internal stop codon;
    medium: untyped introns with splice score >= UNTYPED_MIN_SPLICE; low (not written):
    internal stop codon, invalid start codon, or untyped intron with a weak splice site."""
    table = CodonTable.unambiguous_dna_by_id[gencode]
    prot = _translate(s, exons, gencode)
    if "*" in prot.rstrip("*"):
        return "low", "internal stop codon in spliced CDS"
    if not partial5 and s[exons[0][0] - 1:exons[0][0] + 2] not in table.start_codons:
        return "low", f"no start codon ({s[exons[0][0] - 1:exons[0][0] + 2]})"
    weak = [i + 1 for i, (k, sc) in enumerate(zip(kinds, splice)) if not k and sc < UNTYPED_MIN_SPLICE]
    if weak:
        return "low", f"intron {','.join(map(str, weak))} untyped with weak splice site (< {UNTYPED_MIN_SPLICE})"
    if all(kinds):
        return "high", "all introns typed (Rfam), spliced CDS without internal stop"
    return "medium", f"untyped intron(s) with splice score >= {UNTYPED_MIN_SPLICE}, spliced CDS without internal stop"


# ── driver ───────────────────────────────────────────────────────────────────
def process_record(rec, lib: Dict[str, LibProt], lib_db: Path, a, tools: Dict[str, Optional[str]],
                   tmp: Path) -> List[dict]:
    seq = str(rec.seq).upper()
    L = len(seq)
    fna = tmp / "genome.fna"
    fna.write_text(f">g\n{seq}\n")
    own = norm_species(rec.annotations.get("organism", "")) if a.exclude_own_species else ""
    cands = prescreen(seq, lib, lib_db, a.genetic_code, tools, tmp, own, a.threads)
    rows: List[dict] = []
    models: List[Model] = []
    for gene, c in sorted(cands.items()):
        prot = lib[c["prot"]].seq
        pf = tmp / "prot.faa"
        pf.write_text(f">{c['prot']}\n{prot}\n")
        groups = group_reports(run_exonerate(str(tools["exonerate"]), pf, fna, a.genetic_code, a.maxintron), len(prot))
        m, why = build_model(c["prot"], gene, groups[0], len(prot)) if groups else (None, "")
        # nad4/nad5 (and other paralogues) are homologous: when the best chain overlaps another gene's CDS,
        # the next chain that does not is used (NC_006902.1 ND5); otherwise the best one is reported as low.
        if m is not None and conflicting_genes(rec, m.gene, m.strand, m.exons):
            m = next((alt for alt in (build_model(c["prot"], gene, g, len(prot))[0] for g in groups[1:6])
                      if alt is not None and not conflicting_genes(rec, alt.gene, alt.strand, alt.exons)),
                     None) or realign_masked(rec, seq, c["prot"], gene, groups, len(prot), pf, tmp, a, tools) or m
        if m is None:
            if why:
                rows.append({"record": rec.id, "gene": gene, "written": "no", "note": why})
            continue
        if any(overlap(m.span, k.span) for k in models if k.strand == m.strand):
            continue
        if rec.annotations.get("topology") == "circular" and m.exons[0][0] <= 50 and m.exons[-1][1] >= L - 50:
            rows.append({"record": rec.id, "gene": gene, "note": "gene spans the origin of the circular record - not modelled"})
            continue
        models.append(m)
    rc_seq = str(Seq(seq).reverse_complement())
    raw = {}                                   # every intron of every model, typed in one cmsearch call
    for k, m in enumerate(models):
        use = seq if m.strand == 1 else rc_seq
        tx = to_transcript(m, L)
        raw.update({f"m{k}i{j}": use[tx[j][1]:tx[j + 1][0] - 1] for j in range(len(tx) - 1)})
    typed = intron_types(tools["cmsearch"], a.cm_files, raw, a.evalue, tmp) if tools["cmsearch"] and raw else {}
    for k, m in enumerate(models):
        use = seq if m.strand == 1 else rc_seq
        tx = to_transcript(m, L)
        tinfo: List[Any] = [typed.get(f"m{k}i{j}", ("", None, "")) for j in range(len(tx) - 1)]
        kinds = [t[0] for t in tinfo]
        adj, notes, sscores = adjust_introns(use, tx, kinds, a.genetic_code, homolog_scorer(lib[m.query].seq))
        refs = [lib[m.query].seq] + [v.seq for k_, v in lib.items() if v.gene == m.gene and k_ != m.query
                                     and not (a.exclude_own_species and v.species == own)][:9]
        adj, tinfo, notes, sscores = add_mini_exons(use, adj, tinfo, notes, sscores, a.genetic_code, refs)
        kinds = [t[0] for t in tinfo]
        z_start, z_end = mitoz_ends(rec, m, adj, L)
        started, partial5 = extend_to_start(use, adj, a.genetic_code, prev_feature_end(rec, m, adj, L), m.qmiss_n)
        fallback, partial3 = extend_to_stop(use, started, a.genetic_code, next_feature_start(rec, m, adj, L),
                                            m.qmiss_c, rule=a.stop_rule)
        cand = list(adj)
        cand[0] = (z_start if z_start is not None else fallback[0][0], cand[0][1])
        cand[-1] = (cand[-1][0], z_end if z_end is not None else fallback[-1][1])
        partial5 = partial5 and z_start is None
        partial3 = partial3 and z_end is None
        if valid_exons(cand):
            adj = cand
        else:
            z_start = z_end = None
            adj = fallback
        if not valid_exons(adj):
            rows.append({"record": rec.id, "gene": m.gene, "note": "invalid exon order (gene spans origin?) - not modelled"})
            continue
        ends = ("5'/3' end from existing CDS: " + ", ".join(x for x, v in (("start", z_start), ("end", z_end)) if v)
                if z_start or z_end else "")
        g_exons = to_genome(adj, m.strand, L)
        info: List[Dict[str, Any]] = []
        for j in range(len(adj) - 1):
            s_tx, e_tx = adj[j][1] + 1, adj[j + 1][0] - 1
            gs, ge = (s_tx, e_tx) if m.strand == 1 else (L + 1 - e_tx, L + 1 - s_tx)
            kind, ev, model = tinfo[j]
            label = type_label(kind, ev, model)
            if not kind and ge - gs + 1 > a.maxintron:
                label += " - UNCONFIRMED: long intron from chained alignments, check manually"
            info.append({"start": gs, "end": ge, "note": f"{label}; {notes[j]}; mito_intron.py (MFannot method)"})
        if ends:
            info[0]["note"] += f"; {ends}"
        confidence, reason = evidence(use, adj, kinds, sscores, a.genetic_code, partial5)
        clash = conflicting_genes(rec, m.gene, m.strand, g_exons)
        if clash:
            confidence, reason = "low", f"model overlaps existing {','.join(sorted(set(clash)))} CDS (different gene)"
        hom = lib[m.query]
        name = (apply_model(rec, m, g_exons, info, a.genetic_code, partial5, partial3)
                if confidence != "low" else m.gene)
        for n, x in enumerate(info, 1):
            rows.append({"record": rec.id, "gene": name, "strand": "+" if m.strand == 1 else "-",
                         "cds": ",".join(f"{s}..{e}" for s, e in g_exons), "intron": n,
                         "intron_start": x["start"], "intron_end": x["end"], "intron_len": x["end"] - x["start"] + 1,
                         "confidence": confidence, "written": "no" if confidence == "low" else "yes",
                         "evidence": reason, "note": x["note"], "homolog": f"{hom.species} {hom.gene}",
                         "homolog_coverage": f"{m.qcov:.2f}", "exonerate_score": m.score})
    return rows


_JOB: Dict[str, Any] = {}               # shared state for process_file (set in main; inherited by fork workers)


def process_file(f: Path) -> List[dict]:
    """Annotate every record of one GenBank file; a failing record is reported and left unchanged."""
    lib, lib_db, a, tools, out = _JOB["lib"], _JOB["lib_db"], _JOB["a"], _JOB["tools"], _JOB["out"]
    rows: List[dict] = []
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        try:
            recs = list(SeqIO.parse(str(f), "genbank"))
        except Exception as exc:
            return [{"record": f.name, "written": "no", "note": f"cannot read GenBank file: {exc}"}]
        for rec in recs:
            try:
                rows += process_record(rec, lib, lib_db, a, tools, tmp)
            except Exception as exc:              # keep going: record is written back unchanged
                rows.append({"record": rec.id, "written": "no", "note": f"error, record left unchanged: {exc!r}"})
        src_root = Path(a.input) if Path(a.input).is_dir() else f.parent
        dest = out / f.relative_to(src_root)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if any(r.get("written") == "yes" for r in rows):
            SeqIO.write(recs, str(dest), "genbank")
        elif dest.resolve() != f.resolve():
            shutil.copy2(f, dest)
    return rows


def find_tool(name: str, bin_dir: Optional[str]) -> Optional[str]:
    return mito_env.find(name, bin_dir)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank file or directory (e.g. MitoZ 04.final_gb)")
    p.add_argument("-o", "--out", required=True, help="Output directory (may equal input dir to update in place)")
    p.add_argument("--lib", action="append", default=[],
                   help=f"Protein library in MFannot header format (default {DEFAULT_LIB}); repeatable")
    p.add_argument("--ref_gb", action="append", default=[],
                   help="GenBank file/dir of close relatives: their CDS translations are added to the "
                        "homolog library (helps lineages poorly covered by MFannot's library, e.g. annelids)")
    p.add_argument("--genetic_code", type=int, default=5)
    p.add_argument("--maxintron", type=int, default=MAX_INTRON, help=f"Exonerate --maxintron (MFannot default {MAX_INTRON})")
    p.add_argument("--evalue", type=float, default=0.01, help="cmsearch E-value cut-off for intron type")
    p.add_argument("--cm", default=str(DEFAULT_CM), help="Rfam group I/II intron covariance models")
    p.add_argument("--no_subtype_cms", action="store_true",
                   help="Do not use the group I subtype models of Selles Vidal et al. 2025 (Rfam only)")
    p.add_argument("--bin", default=None, help="Directory with exonerate, blastp, makeblastdb, cmsearch")
    p.add_argument("--threads", type=int, default=1, help="Threads for blastp and cmsearch")
    p.add_argument("--jobs", type=int, default=1, help="GenBank files processed in parallel")
    p.add_argument("--timeout", type=int, default=900, help="Seconds allowed per external program call")
    p.add_argument("--stop_rule", choices=["homolog", "first"], default="first",
                   help="3' end when no MitoZ end is available: 'homolog' = stop where the homolog ends "
                        "(complete or T/TA stop there); 'first' = first in-frame complete stop codon")
    p.add_argument("--exclude_own_species", action="store_true",
                   help="Ignore library proteins of the record's own species (benchmarking)")
    a = p.parse_args(argv)

    tools = {t: find_tool(t, a.bin) for t in ("exonerate", "blastp", "makeblastdb", "cmsearch")}
    missing = [t for t in ("exonerate", "blastp", "makeblastdb") if not tools[t]]
    if missing:
        print(f"[mito_intron] missing {', '.join(missing)} "
              "(conda install -c bioconda exonerate blast infernal; or --bin)", file=sys.stderr)
        return 2
    if not tools["cmsearch"]:
        print("[mito_intron] cmsearch not found: intron type will not be identified", file=sys.stderr)
    a.cm_files = [Path(a.cm)] + ([] if a.no_subtype_cms else subtype_cms())
    lib = load_library(a.lib or [str(DEFAULT_LIB)])
    if a.ref_gb:
        ref = load_ref_gb(a.ref_gb, a.genetic_code)
        lib.update(ref)
        print(f"[mito_intron] + {len(ref)} reference proteins from --ref_gb", file=sys.stderr)
    files = iter_files([a.input], GB_EXT)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    global TIMEOUT, CM_CPU
    TIMEOUT, CM_CPU = a.timeout, max(1, a.threads)
    with tempfile.TemporaryDirectory() as td:
        lib_faa = Path(td) / "lib.faa"
        lib_faa.write_text("".join(f">{k}\n{v.seq}\n" for k, v in lib.items()))
        subprocess.run([str(tools["makeblastdb"]), "-in", str(lib_faa), "-dbtype", "prot", "-out", str(Path(td) / "lib")],
                       capture_output=True, check=True)
        _JOB.update(lib=lib, lib_db=Path(td) / "lib", a=a, tools=tools, out=out)
        if a.jobs > 1 and len(files) > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(a.jobs) as pool:      # workers inherit _JOB (no pickling)
                for r in pool.imap(process_file, files):
                    rows += r
        else:
            for f in files:
                rows += process_file(f)
    cols = ["record", "gene", "strand", "cds", "intron", "intron_start", "intron_end", "intron_len",
            "confidence", "written", "evidence", "note", "homolog", "homolog_coverage", "exonerate_score"]
    with open(out / "intron_report.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    w_ = [r for r in rows if r.get("written") == "yes"]
    low = [r for r in rows if r.get("written") == "no"]
    print(f"[mito_intron] {len(files)} file(s), library {len(lib)} proteins; "
          f"{len(w_)} intron(s) written in {len({(r['record'], r['gene']) for r in w_})} gene(s); "
          f"{len(low)} low-confidence not written -> {out / 'intron_report.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
