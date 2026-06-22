"""
Guard that docs/APPARATUS.md stays in sync with records/*.yaml (it is auto-generated
by tools/gen_apparatus.py). If this fails, run `python tools/gen_apparatus.py`.
"""
from tools.gen_apparatus import _OUT, _records, render


def test_apparatus_doc_in_sync():
    expected = render(_records())
    actual = _OUT.read_text(encoding="utf-8")
    assert actual == expected, "docs/APPARATUS.md is stale — run `python tools/gen_apparatus.py`"


def test_every_record_listed():
    names = [r["name"] for r in _records()]
    doc = _OUT.read_text(encoding="utf-8")
    missing = [n for n in names if ("`%s`" % n) not in doc]
    assert not missing, "records missing from APPARATUS.md: %s" % missing
    # the new inferred apparatus quantities are present
    for n in ("raman_mutual_linewidth_25mg", "raman_beam_path_jitter_25mg",
              "mg_sideband_nbar_oc_lf_ma_25mg"):
        assert n in names
