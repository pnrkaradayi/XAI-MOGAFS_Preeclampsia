"""GEO series matrix dosyaları için ortak okuma yardımcıları."""
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# Windows konsolu (cp1254) Türkçe dışı karakterlerde çöküyor; çıktıyı UTF-8'e zorla.
for _s in (sys.stdout, sys.stderr):
    if hasattr(_s, "reconfigure"):
        _s.reconfigure(encoding="utf-8", errors="replace")

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"


def _split(line):
    return [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]


def read_sample_annotation(path):
    """Series matrix başlığındaki !Sample_* satırlarını örnek x alan tablosuna çevirir.

    `!Sample_characteristics_ch1` satırları "anahtar: değer" biçimindedir ve aynı
    alan birden çok satıra yayılabilir; anahtar adı sütun adı olarak kullanılır.
    """
    fields = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            cells = _split(line)
            fields[cells[0][len("!Sample_"):]].append(cells[1:])

    n = len(fields["geo_accession"][0])
    out = {}
    for name, rows in fields.items():
        if name.startswith("characteristics"):
            for row in rows:
                for i, cell in enumerate(row):
                    if ":" not in cell:
                        continue
                    key, val = cell.split(":", 1)
                    key = f"ch:{key.strip()}"
                    out.setdefault(key, [None] * n)[i] = val.strip()
        else:
            for j, row in enumerate(rows):
                out[name if j == 0 else f"{name}.{j}"] = row
    return pd.DataFrame(out).set_index("geo_accession")


def read_expression(path):
    """Series matrix veri bloğunu (prob x örnek) okur; 'null' -> NaN."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
        df = pd.read_csv(f, sep="\t", index_col=0, na_values=["null", "NA", ""],
                         comment="!", low_memory=False)
    df.index = df.index.astype(str).str.strip('"')
    return df


def to_num(s):
    return pd.to_numeric(s, errors="coerce")


def slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_").lower()
