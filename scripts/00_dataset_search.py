"""
Adım -1 — Sistematik GEO veri seti taraması (tekrarlanabilir arama stratejisi).

NCBI E-utilities (db=gds) üzerinde önceden tanımlı sorgular çalıştırılır, tüm insan
GSE kayıtlarının özetleri çekilir ve anahtar kelimeyle ön sınıflandırılır:
  * kompartman: maternal kan (tam kan, plazma, serum, PBMC, cfRNA, EV) / doku / kordon kanı
  * veri tipi: ifade (dizi / dizileme), ncRNA, metilasyon, diğer
Ön sınıflandırma yalnızca elemeyi hızlandırır; nihai karar elle inceleme ile verilir.

Çıktı: results/dataset_search/geo_hits.tsv  (+ arama tarihi ve sorgular search_log.txt)
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_utils import RESULTS, pd  # noqa: E402

OUT = RESULTS / "dataset_search"
OUT.mkdir(parents=True, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HUMAN = '"Homo sapiens"[ORGN] AND gse[ETYP]'

QUERIES = {
    "PE_HDP": '(preeclampsia OR pre-eclampsia OR preeclamptic OR "hypertensive disorders of pregnancy" '
              'OR "gestational hypertension" OR HELLP OR eclampsia)',
    "PTB_PPROM": '("preterm birth" OR "preterm delivery" OR "preterm labor" OR "preterm labour" OR PPROM '
                 'OR "premature rupture" OR "spontaneous preterm" OR "preterm parturition")',
    "FGR_SGA": '("fetal growth restriction" OR "intrauterine growth restriction" OR "small for gestational age" '
               'OR "fetal growth" OR FGR OR IUGR)',
    "GDM": '("gestational diabetes" OR "gestational diabetes mellitus" OR GDM)',
    "placental_other": '(stillbirth OR "placental abruption" OR "placental insufficiency" OR "placental dysfunction" '
                       'OR "great obstetrical syndromes" OR "adverse pregnancy outcome" OR "pregnancy complication")',
    "cfRNA_EV": '(pregnancy OR pregnant OR gestation OR gestational) AND ("cell-free RNA" OR cfRNA OR '
                '"circulating RNA" OR "plasma RNA" OR exosome OR exosomes OR exosomal OR "extracellular vesicle" '
                'OR "extracellular vesicles" OR "circulating microRNA")',
    "maternal_blood_longitudinal": '(pregnancy OR pregnant) AND ("whole blood" OR "peripheral blood" OR PBMC '
                                   'OR leukocytes OR plasma OR serum OR "maternal blood") AND (longitudinal '
                                   'OR trimester OR "gestational age" OR "early pregnancy")',
}

BLOOD = r"\b(whole blood|peripheral blood|maternal blood|blood|plasma|serum|pbmcs?|leukocytes?|leucocytes?|" \
        r"monocytes?|neutrophils?|lymphocytes?|t cells?|nk cells?|cell-free|cfrna|cfdna|circulating|" \
        r"exosom\w*|extracellular vesicles?|paxgene|tempus)\b"
TISSUE = r"\b(placenta\w*|decidua\w*|villi|villous|trophoblast\w*|myometri\w*|amnion|chorion\w*|" \
         r"fetal membranes?|endometri\w*|umbilical|huvecs?|cell lines?|biops\w*|adipose|muscle|cervi\w*)\b"
CORD = r"\b(cord blood|umbilical cord blood|neonat\w*|newborns?|infants?|offspring)\b"


def get(url, data=None, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "pe-specificity-study/0.1"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # NCBI ara sıra 429/5xx döner
            if k == tries - 1:
                raise
            time.sleep(2 * (k + 1))


def esearch(term):
    q = urllib.parse.urlencode({"db": "gds", "term": f"({term}) AND {HUMAN}", "retmax": 5000, "retmode": "json"})
    r = get(f"{EU}/esearch.fcgi?{q}")["esearchresult"]
    return r["idlist"], int(r["count"])


def esummary(ids):
    out = []
    for i in range(0, len(ids), 200):
        body = urllib.parse.urlencode({"db": "gds", "id": ",".join(ids[i:i + 200]), "retmode": "json"}).encode()
        res = get(f"{EU}/esummary.fcgi", data=body)["result"]
        out += [res[u] for u in res["uids"]]
        time.sleep(0.4)
    return out


def main():
    log = [f"Arama tarihi: {date.today().isoformat()}", f"Veritabanı: NCBI GEO DataSets (db=gds), filtre: {HUMAN}", ""]
    hits = {}
    for name, term in QUERIES.items():
        ids, n = esearch(term)
        log.append(f"[{name}] {n} kayıt\n    ({term}) AND {HUMAN}")
        for u in ids:
            hits.setdefault(u, set()).add(name)
        time.sleep(0.4)
    log.append(f"\nBirleşim: {len(hits)} benzersiz GSE")
    summ = esummary(list(hits))
    rows = []
    for s in summ:
        text = f"{s.get('title', '')} {s.get('summary', '')}".lower()
        rows.append({
            "gse": s.get("accession"),
            "title": s.get("title", "").strip(),
            "n_samples": int(s.get("n_samples") or 0),
            "type": s.get("gdstype", ""),
            "platforms": ";".join("GPL" + g for g in str(s.get("gpl", "")).split(";") if g),
            "pubmed": ";".join(str(p) for p in s.get("pubmedids", [])),
            "pub_date": s.get("pdat", ""),
            "suppl": s.get("suppfile", ""),
            "queries": ";".join(sorted(hits[s["uid"]])),
            "kw_blood": bool(re.search(BLOOD, text)),
            "kw_tissue": bool(re.search(TISSUE, text)),
            "kw_cord": bool(re.search(CORD, text)),
            "kw_methylation": "methylation" in text or "methylation" in s.get("gdstype", "").lower(),
            "summary": re.sub(r"\s+", " ", s.get("summary", ""))[:1500],
        })
    df = pd.DataFrame(rows).sort_values("n_samples", ascending=False)
    df.to_csv(OUT / "geo_hits.tsv", sep="\t", index=False)
    blood = df[df.kw_blood]
    log.append(f"Kan anahtar kelimesi geçen: {len(blood)}; bunların doku anahtar kelimesi içermeyen: "
               f"{(~blood.kw_tissue).sum()}")
    (OUT / "search_log.txt").write_text("\n".join(log), encoding="utf-8")
    print("\n".join(log))


if __name__ == "__main__":
    main()
