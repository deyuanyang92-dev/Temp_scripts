#!/usr/bin/env python3
"""Run MITOS2 (runmitos) without its sensitive-mode fallback for RNA genes.

MITOS2 first searches all RNA models in Infernal fast mode, then re-runs every model that gave no
hit in sensitive mode (cmsearch --max, E <= 0.1; mitos/scripts/runmitos.py lines ~255-300). For
genomes that lack most tRNAs (e.g. anthozoans with 1-2 tRNAs) this fallback runs ~20 models
unfiltered on two copies of the sequence and dominates the run time (Metridium senile: >10 min),
and its relaxed threshold also yields weak hits. Here mitfi.remove_empty() - the test that triggers
the fallback - is replaced to always return False; everything else is MITOS2 unchanged.

Usage: <MITOS2 env>/bin/python runmitos_fast.py <runmitos arguments>
"""
import sys

from mitos import mitfi
from mitos.scripts import runmitos

mitfi.remove_empty = lambda *args, **kwargs: False
runmitos.mitfi.remove_empty = mitfi.remove_empty
sys.exit(runmitos.main())
