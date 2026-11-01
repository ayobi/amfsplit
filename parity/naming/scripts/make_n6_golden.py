#!/usr/bin/env python3
"""Gate 0b step 2 golden: the archived species-clade rows plus the frozen ITS distance rule (N6).

Applies the frozen rule (amf-novelty-benchmark release 0.3.1, benchmarks/m4_external_holdout/its_rule.py,
pool 'all', 99th percentile) to the archived epa_species_clades_* tables fold by fold, exactly as
external_holdout_v2.score_v2 does: training = every culture except the fold's held cultures, and a
firing row takes its placement from the original row's call_genera. Writes
golden/epa_species_clades_n6_{physical_culture,physical_species}.copies.tsv.

    python3 parity/naming/scripts/make_n6_golden.py
"""
import copy
from collections import defaultdict
import hashlib
import shutil
import tempfile
from pathlib import Path
import sys

FIX = Path(__file__).resolve().parent.parent
# bridge.py checks the arm scripts beside itself, as in the frozen M2b layout (code/baseline), so the
# frozen sources are assembled into one temporary folder rather than imported from three.
CODE = Path(tempfile.mkdtemp(prefix='frozen_m2b_'))
for src in [*(FIX / 'frozen_sources').glob('*.py'), *(FIX / 'frozen_sources/m2b').glob('*.py'),
            FIX / 'frozen_sources/m4/its_rule.py']:
    shutil.copy2(src, CODE / src.name)
sys.path.insert(0, str(CODE))
from bridge import Data                      # noqa: E402  frozen M2b loader; checks reference and evidence hashes
from lsu_placement import read_tsv, write_tsv  # noqa: E402  frozen TSV reader and writer
from its_rule import ItsRule                 # noqa: E402  frozen N6

ITS_RULE_SHA = 'de0feb459499755ea955452fa25d46830ec5c3098e72c03f7462c8823ef8dddf'
SOURCES = {'physical_culture': 'a357b563a5def31932085dabe4697f3b49e0f586d56a6c91acefe8bb1cb7cfd4',
           'physical_species': '9b0d4ba422214a763cbf0c5949f8f7b2b163fdfcd8e2fd5ebf81a55db521a53a'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


if sha(FIX / 'frozen_sources/m4/its_rule.py') != ITS_RULE_SHA:
    sys.exit('frozen its_rule.py hash mismatch')
data = Data(FIX / 'inputs/benchmark')
rule = ItsRule(data, 99)
for mode, expected in SOURCES.items():
    source = FIX / 'golden' / f'epa_species_clades_{mode}.copies.tsv'
    if sha(source) != expected:
        sys.exit(f'archived species-clade table changed: {source.name}')
    rows = read_tsv(source)
    check = FIX / 'golden' / f'.roundtrip_{mode}.tsv'          # the frozen writer must reproduce the archive
    write_tsv(check, rows)
    same = sha(check) == expected
    check.unlink()
    if not same:
        sys.exit(f'frozen writer does not reproduce {source.name}')
    folds = defaultdict(list)
    for r in rows:
        folds[r['held_cultures']].append(r)
    fired = 0
    for held, fold in folds.items():
        training = data.cultures - frozenset(h for h in held.split(';') if h)
        before = {r['accession']: copy.deepcopy(r) for r in fold}
        fired += rule.apply(fold, before, training)          # pool 'all', as in score_v2
    out = FIX / 'golden' / f'epa_species_clades_n6_{mode}.copies.tsv'
    write_tsv(out, rows)
    known = [r for r in rows if r['copy_class'] == 'main' and r['status'].startswith('known')]
    print(f'{mode}: {len(folds)} holdout rounds, {len(rows)} rows, {fired} firings; '
          f'main-type known calls left {len(known)}; {out.name} {sha(out)}')
