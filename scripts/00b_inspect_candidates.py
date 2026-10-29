"""
Aday GEO serilerinin örnek düzeyi anotasyonunu özetler (doku, hafta, gruplar).

Her seri için GEO'dan tüm GSM kayıtlarının kısa görünümü çekilir; her
`characteristics` anahtarının değer dağılımı ve source_name yazdırılır.

Kullanım: python scripts/00b_inspect_candidates.py GSE46510 GSE108497 ...
Çıktı:    results/dataset_search/characteristics/<GSE>.txt
"""
import collections
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS  # noqa: E402

OUT = RESULTS / "dataset_search" / "characteristics"
OUT.mkdir(parents=True, exist_ok=True)
ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ={}&form=text&view=brief"


def fetch(acc, targ):
    for k in range(4):
        try:
            req = urllib.request.Request(ACC.format(acc, targ), headers={"User-Agent": "pe-specificity-study/0.1"})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception:
            if k == 3:
                raise
            time.sleep(3 * (k + 1))


def summarize(gse):
    series = fetch(gse, "self")
    keep = ("!Series_title", "!Series_pubmed_id", "!Series_overall_design", "!Series_relation",
            "!Series_platform_id", "!Series_supplementary_file")
    head = [l for l in series.splitlines() if l.startswith(keep)]
    gsm = fetch(gse, "gsm")
    samples, cur = [], None
    for line in gsm.splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"_chars": []}
            samples.append(cur)
        elif cur is None or " = " not in line:
            continue
        else:
            k, v = line.split(" = ", 1)
            if k.startswith("!Sample_characteristics"):
                cur["_chars"].append(v)
            elif k in ("!Sample_source_name_ch1", "!Sample_title", "!Sample_platform_id", "!Sample_molecule_ch1"):
                cur[k] = v
    out = [f"######## {gse}  ({len(samples)} örnek)"] + head
    for k in ("!Sample_source_name_ch1", "!Sample_platform_id", "!Sample_molecule_ch1"):
        c = collections.Counter(s.get(k, "") for s in samples)
        out.append(f"{k}: {dict(c.most_common(8))}")
    out.append("!Sample_title örnek: " + " | ".join(s.get("!Sample_title", "") for s in samples[:6]))
    keyed = collections.defaultdict(collections.Counter)
    for s in samples:
        for ch in s["_chars"]:
            key, _, val = ch.partition(":")
            keyed[key.strip() if _ else "(anahtarsız)"][val.strip() if _ else key.strip()] += 1
    for key, c in keyed.items():
        vals = dict(c.most_common(12))
        more = f" … (+{len(c) - 12} değer)" if len(c) > 12 else ""
        out.append(f"  [{key}] n_değer={len(c)}: {vals}{more}")
    text = "\n".join(out)
    (OUT / f"{gse}.txt").write_text(text, encoding="utf-8")
    return text


if __name__ == "__main__":
    for g in sys.argv[1:]:
        try:
            print(summarize(g)[:6000])
        except Exception as e:
            print(f"######## {g}: HATA {e}")
        print()
        time.sleep(0.5)
