#!/usr/bin/env python3
"""Local web interface for genome_inventory.py (standard library only).

    python3 webapp.py                       # http://127.0.0.1:8765, results in ./genome_inventory_results
    python3 webapp.py --data-dir /mnt/p/genome_inventory --port 8800

The page runs the same command-line tool (so CLI and web always give identical results), shows the progress
log, renders the result tables with search/filter, and offers the Excel/TSV files for download.
Listens on 127.0.0.1 only.
"""
import argparse
import json
import re
import subprocess
import sys
import threading
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOL = HERE / "genome_inventory.py"
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ._-]{0,80}$")
TABLES = ["assemblies", "species", "mitogenomes", "qc", "runs", "changes", "goat_pipeline"]
JOBS: "dict[str, subprocess.Popen]" = {}   # taxon -> running job
LOCK = threading.Lock()

PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Genome Inventory</title>
<style>
:root{--bg:#f6f7f9;--card:#fff;--ink:#1d232b;--mute:#5f6b7a;--line:#dde2e8;--acc:#1f6feb;--warn:#b54708;--ok:#067647}
@media (prefers-color-scheme:dark){:root{--bg:#0f1318;--card:#171c22;--ink:#e6eaef;--mute:#9aa6b2;--line:#2a323c;--acc:#58a6ff;--warn:#f79009;--ok:#32d583}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
header{padding:16px 20px;border-bottom:1px solid var(--line);background:var(--card)}h1{margin:0;font-size:18px}
header p{margin:4px 0 0;color:var(--mute)}main{padding:16px 20px;display:grid;grid-template-columns:minmax(0,1fr);gap:16px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px}
.row{display:flex;flex-wrap:wrap;gap:10px;align-items:center}
input[type=text]{padding:7px 9px;border:1px solid var(--line);border-radius:6px;background:var(--bg);color:var(--ink);min-width:220px}
button{padding:7px 14px;border:0;border-radius:6px;background:var(--acc);color:#fff;cursor:pointer;font-weight:600}
button.sec{background:transparent;color:var(--acc);border:1px solid var(--line)}
label{color:var(--mute)}pre{max-height:220px;overflow:auto;background:var(--bg);padding:8px;border-radius:6px;margin:8px 0 0;font-size:12px}
.tabs{display:flex;gap:6px;flex-wrap:wrap;margin:8px 0}.tabs button{background:transparent;color:var(--ink);border:1px solid var(--line)}
.tabs button.on{background:var(--acc);color:#fff;border-color:var(--acc)}
.tbl{overflow:auto;max-height:65vh;border:1px solid var(--line);border-radius:6px}
table{border-collapse:collapse;font-size:12px;white-space:nowrap}th,td{padding:4px 8px;border-bottom:1px solid var(--line);text-align:left}
th{position:sticky;top:0;background:var(--card);cursor:pointer}td.sub{color:var(--warn);font-weight:600}td.same{color:var(--ok)}
a{color:var(--acc)}.mute{color:var(--mute)}
</style></head><body>
<header><h1>NCBI genome inventory</h1><p>Assemblies, raw data (HiFi / Hi-C / RNA / ONT …), same-sample mitogenomes and taxonomy for any taxon.</p></header>
<main>
<section class="card"><div class="row">
<input type="text" id="taxon" placeholder="Taxon name or taxid, e.g. Annelida">
<label><input type="checkbox" id="noworms"> skip WoRMS</label>
<label><input type="checkbox" id="allv"> all assembly versions</label>
<label><input type="checkbox" id="update" checked> update (new assemblies only)</label>
<label><input type="checkbox" id="refresh"> full refresh</label>
<button id="run">Run</button></div>
<pre id="log" hidden></pre></section>
<section class="card"><div class="row"><b>Results</b><span class="mute" id="rmsg"></span></div><div id="results" class="row" style="margin-top:8px"></div></section>
<section class="card" id="view" hidden>
<div class="row"><b id="vtitle"></b><span id="files"></span></div>
<div class="tabs" id="tabs"></div>
<div class="row"><input type="text" id="filter" placeholder="Filter rows (any column)"><span class="mute" id="count"></span>
<label><input type="checkbox" id="keycols" checked> key columns only</label></div>
<div class="tbl" style="margin-top:8px"><table id="tbl"></table></div></section>
</main>
<script>
const KEY={assemblies:["family","genus","species","Assembly Accession","Assembly Paired Assembly Accession","representative_derived","found_in_derived","assembly_content_derived","Assembly Level","assembly_type","Assembly Release Date","Assembly Submitter","Assembly Sequencing Tech","Assembly Stats Total Sequence Length","Assembly Stats Contig N50","Assembly Stats Scaffold N50","data_types_derived","link_evidence_derived","same_individual_types_derived","other_individual_types_derived","HiFi_runs","Hi-C_runs","RNA-seq_runs","ONT_runs","mito_accession","mito_length_bp","mito_source_derived","mito_note_derived"],
species:["class","order","family","species","organism_name","n_assemblies","assemblies","representative_assembly","rep_assembly_level","assembly_linked_types_derived","mito_accession","mito_source_derived","mito_note_derived","mito_from_assembly","goat_sequencing_status"],
mitogenomes:null,qc:null,changes:null,goat_pipeline:null,runs:["species_tax_id","Run","ScientificName","category_derived","category_confidence_derived","category_rule_derived","linked_assemblies_derived","individual_derived","bases_Gb","BioProject","BioSample","Platform","Model","LibraryStrategy","experiment_title","library_protocol"]};
let cur=null,tab="assemblies",data={},sortc=null,sortd=1;
const $=id=>document.getElementById(id);
async function j(u,o){const r=await fetch(u,o);if(!r.ok)throw new Error(await r.text());return r.json()}
async function listResults(){const d=await j("/api/results");$("results").innerHTML="";$("rmsg").textContent=d.length?"":" none yet";
 d.forEach(t=>{const b=document.createElement("button");b.className="sec";b.textContent=t;b.onclick=()=>openRes(t);$("results").appendChild(b)})}
async function openRes(t){cur=t;data={};$("view").hidden=false;$("vtitle").textContent=t;
 const f=await j("/api/files?taxon="+encodeURIComponent(t));
 $("files").innerHTML=f.map(x=>` <a href="/files/${encodeURIComponent(t)}/${encodeURIComponent(x)}">${x}</a>`).join(" ·");
 $("tabs").innerHTML="";["assemblies","species","mitogenomes","qc","runs","changes","goat_pipeline"].forEach(n=>{const b=document.createElement("button");b.textContent=n;b.onclick=()=>show(n);b.id="tab-"+n;$("tabs").appendChild(b)});show("assemblies")}
function parseTSV(s){const L=s.replace(/\r/g,"").split("\n").filter(x=>x.length);const h=L.shift().split("\t");return{h,rows:L.map(l=>l.split("\t"))}}
async function show(n){tab=n;sortc=null;document.querySelectorAll(".tabs button").forEach(b=>b.classList.toggle("on",b.id=="tab-"+n));
 if(!data[n]){const r=await fetch(`/api/table?taxon=${encodeURIComponent(cur)}&name=${n}`);data[n]=r.ok?parseTSV(await r.text()):{h:["(missing)"],rows:[]}}render()}
function render(){const d=data[tab];if(!d)return;let cols=d.h.map((c,i)=>i);
 if($("keycols").checked&&KEY[tab]){const k=KEY[tab].map(c=>d.h.indexOf(c)).filter(i=>i>=0);if(k.length)cols=k}
 const q=$("filter").value.toLowerCase();let rows=d.rows.filter(r=>!q||r.join("\t").toLowerCase().includes(q));
 if(sortc!==null)rows=rows.slice().sort((a,b)=>{const x=a[sortc]||"",y=b[sortc]||"";const nx=parseFloat(x),ny=parseFloat(y);
  return(!isNaN(nx)&&!isNaN(ny)?nx-ny:x.localeCompare(y))*sortd});
 $("count").textContent=`${rows.length} / ${d.rows.length} rows`;
 const esc=s=>String(s).replace(/[&<>]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));
 let h="<thead><tr>"+cols.map(i=>`<th data-i="${i}">${esc(d.h[i])}</th>`).join("")+"</tr></thead><tbody>";
 rows.slice(0,3000).forEach(r=>{h+="<tr>"+cols.map(i=>{const v=r[i]||"";const c=d.h[i]=="mito_source_derived"?(v.startsWith("SUBSTITUTE")?"sub":v?"same":""):"";return`<td class="${c}">${esc(v.length>120?v.slice(0,120)+"…":v)}</td>`}).join("")+"</tr>"});
 $("tbl").innerHTML=h+"</tbody>";$("tbl").querySelectorAll("th").forEach(th=>th.onclick=()=>{const i=+th.dataset.i;sortd=sortc===i?-sortd:1;sortc=i;render()})}
$("filter").oninput=render;$("keycols").onchange=render;
$("run").onclick=async()=>{const t=$("taxon").value.trim();if(!t)return;$("log").hidden=false;$("log").textContent="starting…";
 try{await j("/api/run",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify({taxon:t,no_worms:$("noworms").checked,all_versions:$("allv").checked,update:$("update").checked,refresh:$("refresh").checked})})}
 catch(e){$("log").textContent=e.message;return}poll(t)};
async function poll(t){const s=await j("/api/status?taxon="+encodeURIComponent(t));$("log").textContent=s.log;$("log").scrollTop=1e9;
 if(s.running)setTimeout(()=>poll(t),2000);else{await listResults();if(s.ok)openRes(s.taxon)}}
listResults();
</script></body></html>"""


def safe_name(s):
    s = (s or "").strip()
    if not NAME_RE.match(s):
        raise ValueError("invalid taxon name")
    return s.replace(" ", "_")


class Handler(BaseHTTPRequestHandler):
    data_dir = Path("genome_inventory_results")

    def log_message(self, fmt, *a):            # quiet
        pass

    def send(self, code, body, ctype="application/json; charset=utf-8", extra=None):
        b = body if isinstance(body, bytes) else body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b)))
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(b)

    def result_dir(self, taxon):
        d = (self.data_dir / safe_name(taxon)).resolve()
        if self.data_dir.resolve() not in d.parents:
            raise ValueError("bad path")
        return d

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        q = dict(urllib.parse.parse_qsl(u.query))
        try:
            if u.path == "/":
                return self.send(200, PAGE, "text/html; charset=utf-8")
            if u.path == "/api/results":
                names = sorted(p.name for p in self.data_dir.glob("*") if p.is_dir() and any(p.glob("*.xlsx")))
                return self.send(200, json.dumps(names))
            if u.path == "/api/files":
                d = self.result_dir(q.get("taxon"))
                return self.send(200, json.dumps(sorted(p.name for p in d.iterdir() if p.is_file()
                                                        and p.suffix in (".xlsx", ".tsv", ".log"))))
            if u.path == "/api/table":
                name = q.get("name")
                if name not in TABLES:
                    raise ValueError("unknown table")
                d = self.result_dir(q.get("taxon"))
                hits = sorted(d.glob(f"*_{name}.tsv"))          # <Taxon>_<name>.tsv, folder name may differ
                f = hits[0] if hits else d / "missing"
                if not f.exists():
                    return self.send(404, "missing", "text/plain")
                return self.send(200, f.read_bytes(), "text/tab-separated-values; charset=utf-8")
            if u.path == "/api/status":
                taxon = safe_name(q.get("taxon"))
                d = self.result_dir(taxon)
                log = (d / "run.log").read_text(encoding="utf-8", errors="replace") if (d / "run.log").exists() else ""
                with LOCK:
                    p = JOBS.get(taxon)
                    running = p is not None and p.poll() is None
                    ok = p is not None and p.poll() == 0
                return self.send(200, json.dumps({"taxon": taxon, "running": running, "ok": ok,
                                                  "log": log[-20000:]}))
            if u.path.startswith("/files/"):
                parts = [urllib.parse.unquote(x) for x in u.path.split("/")[2:]]
                if len(parts) != 2 or "/" in parts[1] or parts[1].startswith("."):
                    raise ValueError("bad path")
                f = self.result_dir(parts[0]) / parts[1]
                if not f.is_file():
                    return self.send(404, "missing", "text/plain")
                ctype = ("application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                         if f.suffix == ".xlsx" else "text/plain; charset=utf-8")
                return self.send(200, f.read_bytes(), ctype,
                                 {"Content-Disposition": f'attachment; filename="{f.name}"'})
            return self.send(404, "not found", "text/plain")
        except ValueError as e:
            return self.send(400, str(e), "text/plain")

    def do_POST(self):
        if self.path != "/api/run":
            return self.send(404, "not found", "text/plain")
        try:
            body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
            taxon = safe_name(body.get("taxon"))
        except (ValueError, json.JSONDecodeError) as e:
            return self.send(400, str(e), "text/plain")
        d = self.result_dir(taxon)
        d.mkdir(parents=True, exist_ok=True)
        with LOCK:
            p = JOBS.get(taxon)
            if p is not None and p.poll() is None:
                return self.send(409, "already running", "text/plain")
            cmd = [sys.executable, str(TOOL), body["taxon"].strip(), "-o", str(d)]
            cmd += ["--no-worms"] if body.get("no_worms") else []
            cmd += ["--all-versions"] if body.get("all_versions") else []
            cmd += ["--update"] if body.get("update") else []
            cmd += ["--refresh"] if body.get("refresh") else []
            log = open(d / "run.log", "w", encoding="utf-8")
            JOBS[taxon] = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=str(HERE))
        return self.send(200, json.dumps({"started": taxon}))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default="genome_inventory_results", help="where result folders live")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args(argv)
    Handler.data_dir = Path(a.data_dir)
    Handler.data_dir.mkdir(parents=True, exist_ok=True)
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), Handler)
    print(f"Genome inventory web UI: http://127.0.0.1:{a.port}  (results in {Handler.data_dir.resolve()})")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
