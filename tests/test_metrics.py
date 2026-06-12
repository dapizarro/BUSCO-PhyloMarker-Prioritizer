import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'compute_alignment_metrics.py'
spec = importlib.util.spec_from_file_location('compute_alignment_metrics', SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_compute_metrics_basic():
    seqs = {'tax1':'ACDE-', 'tax2':'ACDF-', 'tax3':'ACDFG'}
    m = mod.compute_metrics(seqs)
    assert m['n_taxa'] == 3
    assert m['aln_len'] == 5
    assert m['variable_sites'] >= 1
    assert 0 <= m['occupancy'] <= 1
