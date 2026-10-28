"""Minimal HTTP clients (standard library only): NCBI E-utilities, NCBI Datasets v2, WoRMS REST.

- Identity: NCBI_EMAIL / NCBI_API_KEY environment variables (optional; a key makes E-utilities ~3x faster).
- Rate limit: 0.11 s/request with a key, 0.35 s without (NCBI policy).
- Failures: HTTP errors, HTML error pages, JSON "error", XML <ERROR> all count as failures and are
  retried with exponential back-off; the last error is raised as NcbiError.
- `transport` is injectable for offline tests: callable(url, data=None) -> str.

E-utilities part adapted from yzz_taxonomy_workflow/modules/sra/eutils.py (same author).
"""
import csv
import io
import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/"
WORMS = "https://www.marinespecies.org/rest/"


class NcbiError(RuntimeError):
    pass


def xml_root(text):
    """Parse NCBI XML with the standard library; refuse <!ENTITY (the only entity-expansion entry point)."""
    if "<!ENTITY" in text:
        raise ValueError("XML contains <!ENTITY declaration; refused")
    return ET.fromstring(text)


def _urlopen(url, data=None):
    req = urllib.request.Request(url, data=data, headers={"Accept": "application/json, text/xml, */*"})
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read().decode("utf-8", errors="replace")


def _error_in(body, expect=None):
    s = body.lstrip()
    if expect == "json":
        try:
            j = json.loads(s)
        except ValueError:
            return "empty or invalid JSON"
        if isinstance(j, dict) and j.get("error"):
            return str(j["error"])[:200]
        if isinstance(j, dict) and isinstance(j.get("esearchresult"), dict) and j["esearchresult"].get("ERROR"):
            return str(j["esearchresult"]["ERROR"])[:200]
        return None
    if expect == "runinfo":
        rows = [r for r in csv.reader(io.StringIO(s.lstrip("﻿"))) if r]
        if not rows or rows[0][0] != "Run":
            return "empty or invalid RunInfo page"
        if any(len(r) != len(rows[0]) for r in rows):
            return "truncated RunInfo page"
        return None
    if s[:15].upper().startswith("<!DOCTYPE HTML") or s[:5].lower() == "<html":
        return "HTML error page"
    if "<ERROR>" in s[:3000]:
        return s[s.find("<ERROR>"):][:200]
    if expect == "xml":
        try:
            xml_root(s)
        except Exception as e:  # truncated / empty XML
            return f"invalid XML: {e}"
    return None


class Client:
    def __init__(self, email=None, api_key=None, transport=None, retries=5, sleep=time.sleep):
        self.email = email if email is not None else os.environ.get("NCBI_EMAIL")
        self.api_key = api_key if api_key is not None else os.environ.get("NCBI_API_KEY")
        self.transport = transport or _urlopen
        self.retries, self.sleep = retries, sleep
        self.gap = 0.11 if self.api_key else 0.35

    # ---------- generic ----------
    def _get(self, url, data=None, expect=None, gap=None):
        last = None
        for i in range(self.retries):
            self.sleep(self.gap if gap is None else gap)
            try:
                body = self.transport(url, data) if data else self.transport(url)
                err = _error_in(body, expect)
                if err is None:
                    return body
                last = err
            except Exception as e:  # network / HTTP status
                if getattr(e, "code", None) == 204:      # WoRMS: no content = no match
                    return ""
                last = f"{type(e).__name__}: {e}"
            if i < self.retries - 1:
                self.sleep(2 ** i)
        raise NcbiError(f"{url[:120]} failed {self.retries}x: {last}")

    # ---------- E-utilities ----------
    def eutils(self, endpoint, params, expect=None):
        p = {k: v for k, v in params.items() if v not in (None, "")}
        p["tool"] = "ncbi_genome_inventory"
        if self.email:
            p["email"] = self.email
        if self.api_key:
            p["api_key"] = self.api_key
        query = urllib.parse.urlencode(p)
        # long terms (hundreds of OR clauses) must be POSTed, otherwise HTTP 414
        if len(query) > 1500:
            return self._get(EUTILS + endpoint, data=query.encode(), expect=expect)
        return self._get(EUTILS + endpoint + "?" + query, expect=expect)

    def esearch(self, db, term, retmax=0):
        body = self.eutils("esearch.fcgi", {"db": db, "term": term, "usehistory": "y",
                                            "retmax": retmax, "retmode": "json"}, expect="json")
        r = json.loads(body)["esearchresult"]
        return {"count": int(r.get("count", 0)), "webenv": r.get("webenv", ""),
                "query_key": r.get("querykey", ""), "ids": r.get("idlist", [])}

    def efetch_history(self, db, hist, rettype="", retmode="xml", page=500):
        expect = "runinfo" if rettype == "runinfo" else "xml" if retmode == "xml" else None
        for start in range(0, hist["count"], page):
            yield self.eutils("efetch.fcgi", {"db": db, "WebEnv": hist["webenv"], "query_key": hist["query_key"],
                                              "retstart": start, "retmax": page,
                                              "rettype": rettype, "retmode": retmode}, expect=expect)

    def esummary_history(self, db, hist, page=400):
        for start in range(0, hist["count"], page):
            body = self.eutils("esummary.fcgi", {"db": db, "WebEnv": hist["webenv"], "query_key": hist["query_key"],
                                                 "retstart": start, "retmax": page, "retmode": "json"},
                               expect="json")
            yield json.loads(body)

    def efetch_ids(self, db, ids, rettype="", retmode="xml"):
        return self.eutils("efetch.fcgi", {"db": db, "id": ",".join(ids), "rettype": rettype, "retmode": retmode},
                           expect="xml" if retmode == "xml" else None)

    # ---------- Datasets v2 ----------
    def datasets_pages(self, path, params=None, page_size=1000, max_pages=None):
        """Yield every `reports` item of a paginated Datasets v2 endpoint (at most max_pages pages)."""
        params = dict(params or {})
        params["page_size"] = page_size
        token = None
        pages = 0
        while max_pages is None or pages < max_pages:
            pages += 1
            if token:
                params["page_token"] = token
            url = DATASETS + path + "?" + urllib.parse.urlencode(params, doseq=True)
            j = json.loads(self._get(url, expect="json", gap=0.25))
            yield from j.get("reports", [])
            token = j.get("next_page_token")
            if not token:
                return

    def get_text(self, url):
        """Plain HTTPS GET with retries (NCBI FTP files over HTTPS)."""
        return self._get(url, gap=0.1)

    # ---------- WoRMS ----------
    def worms_match(self, names):
        """AphiaRecordsByMatchNames (<= 50 names per call). Returns one list of candidate records per name."""
        q = urllib.parse.urlencode([("scientificnames[]", n) for n in names] + [("marine_only", "false")])
        body = self._get(WORMS + "AphiaRecordsByMatchNames?" + q, gap=0.2)
        return json.loads(body) if body.strip() else [[] for _ in names]


def resolve_taxid(cl, name):
    """Taxon name -> NCBI taxid; digits returned as-is; missing or ambiguous names raise with candidates."""
    name = name.strip()
    if name.isdigit():
        return name
    r = cl.esearch("taxonomy", f'"{name}"[Scientific Name]', retmax=20)
    if r["count"] == 0:
        raise NcbiError(f"NCBI Taxonomy has no taxon {name!r}; check spelling or pass --taxid")
    if r["count"] == 1 and len(r["ids"]) == 1:
        return r["ids"][0]
    root = xml_root(cl.efetch_ids("taxonomy", r["ids"]))
    lines = [f"  {t.findtext('TaxId')}  {t.findtext('ScientificName')}  ({t.findtext('Lineage', '')[-80:]})"
             for t in root.findall("Taxon")]
    raise NcbiError(f"{name!r} matches {r['count']} NCBI taxa, pass --taxid:\n" + "\n".join(lines))
