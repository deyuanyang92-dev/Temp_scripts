#!/usr/bin/env python3
"""Run MITOS2 (runmitos) with its Infernal searches in parallel, optionally without its sensitive pass.

MITOS2 first searches all RNA models in Infernal fast mode, then re-runs every model that gave no
hit in sensitive mode (cmsearch --max, E <= 0.1; mitos/scripts/runmitos.py lines ~255-300). Both
passes run one cmsearch per model, one after the other, on one CPU (mitos/mitfi/__init__.py
cmsearch(): "--cpu 1  # enforce serial processing"). For genomes that lack most tRNAs (e.g.
anthozoans with 1-2 tRNAs) the sensitive pass runs ~20 models unfiltered and dominates the run time
(Metridium senile: >10 min).

Here mitfi.cmsearch() is wrapped: a dry pass records the exact cmsearch command lines MITOS2 would
run (its own arguments, nothing rebuilt here), they are run in parallel (--jobs) into the .cmout files
MITOS2 expects, then the original cmsearch() runs and, finding every .cmout present, only runs MiTFi
on them. Results are those of runmitos.

  --no_sensitive   also skip the sensitive pass (mitfi.remove_empty() always False): fastest, but
                   misses weak tRNAs (2026-10-20 benchmark: 28 of 28 tRNAs that only full MITOS2
                   found, e.g. nematode NC_008047.1 trnD/G/N/R/I/S)
  --jobs N         parallel cmsearch processes (default: all CPUs, max 16)

Usage: <MITOS2 env>/bin/python runmitos_fast.py [--no_sensitive] [--jobs N] <runmitos arguments>
"""
import os
import sys
from concurrent.futures import ThreadPoolExecutor

from mitos import mitfi
from mitos.scripts import runmitos

ORIGINAL_CMSEARCH = mitfi.cmsearch


def pop_option(argv, name, has_value):
    """Remove our own option from argv (runmitos must not see it); return its value (True for a flag)."""
    if name not in argv:
        return None
    i = argv.index(name)
    value = argv[i + 1] if has_value else True
    del argv[i:i + (2 if has_value else 1)]
    return value


def parallel_cmsearch(jobs):
    def cmsearch(seqfile, outpath, *args, **kwargs):
        recorded, created = [], []
        real_execute = mitfi.util.execute
        before = set(os.listdir(outpath))

        def record(cmd, *a, **k):
            if cmd and cmd[0] == "cmsearch":
                recorded.append(list(cmd))
            return b"", b""

        mitfi.util.execute = record
        try:
            ORIGINAL_CMSEARCH(seqfile, outpath, *args, **kwargs)     # dry pass: writes empty files
        finally:
            mitfi.util.execute = real_execute
        created = [f for f in os.listdir(outpath) if f not in before]
        for f in created:
            os.remove(os.path.join(outpath, f))

        def run(cmd):
            model, seq = cmd[-2], cmd[-1]
            name = os.path.splitext(os.path.basename(model))[0]
            out, _ = real_execute(cmd, mitfi.InfernalException)
            with open("%s/%s_%s.cmout" % (outpath, os.path.basename(seq), name), "wb") as fh:
                fh.write(out)

        with ThreadPoolExecutor(jobs) as ex:
            list(ex.map(run, recorded))
        return ORIGINAL_CMSEARCH(seqfile, outpath, *args, **kwargs)  # every .cmout present: MiTFi only
    return cmsearch


if __name__ == "__main__":
    argv = sys.argv
    if pop_option(argv, "--no_sensitive", False):
        mitfi.remove_empty = lambda *args, **kwargs: False
        runmitos.mitfi.remove_empty = mitfi.remove_empty
    jobs = int(pop_option(argv, "--jobs", True) or min(os.cpu_count() or 1, 16))
    mitfi.cmsearch = parallel_cmsearch(jobs)
    runmitos.mitfi.cmsearch = mitfi.cmsearch
    sys.exit(runmitos.main())
