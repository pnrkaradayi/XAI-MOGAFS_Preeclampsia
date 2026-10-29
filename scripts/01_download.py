"""
Adım 0 — Veri indirme.

Tüm dosyalar NCBI GEO'dan, hazır işlenmiş matris olarak indirilir.
Her dosya için boyut ve SHA-256 özeti logs/download_manifest.tsv'ye yazılır.

SON TEST SETİ (GSE192902 Validation-2) data/sealed/ altına indirilir, salt-okunur
yapılır ve SHA-256'sı kaydedilir. Hiçbir analiz betiği bu klasörü okumaz;
yalnızca en son adımda (Adım 7) bir kez açılır. Özet değeri, dosyanın o güne
kadar değişmediğinin kanıtıdır.

Kullanım:
    .venv/Scripts/python.exe scripts/01_download.py            # hepsi
    .venv/Scripts/python.exe scripts/01_download.py GSE85307   # tek seri
"""
import hashlib
import os
import stat
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PLATFORMS = ROOT / "data" / "platforms"
SEALED = ROOT / "data" / "sealed"
MANIFEST = ROOT / "logs" / "download_manifest.tsv"

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def series_url(gse, kind, fname):
    stub = gse[:-3] + "nnn"
    return f"{GEO}/series/{stub}/{gse}/{kind}/{fname}"


def platform_table_url(gpl):
    # Yalnızca platform tablosu (örnek verisi olmadan); GEOquery'nin kullandığı uç nokta.
    return f"{ACC}?acc={gpl}&targ=self&form=text&view=data"


# (grup, hedef klasör, dosya adı, url, rol)
FILES = [
    # --- EĞİTİM: GSE85307 (VDAART, tam kan, aynı kadınlarda mRNA + miRNA)
    ("GSE85307", RAW / "GSE85307", "GSE85307-GPL6244_series_matrix.txt.gz",
     series_url("GSE85307", "matrix", "GSE85307-GPL6244_series_matrix.txt.gz"), "train_mRNA"),
    ("GSE85307", RAW / "GSE85307", "GSE85307-GPL19162_series_matrix.txt.gz",
     series_url("GSE85307", "matrix", "GSE85307-GPL19162_series_matrix.txt.gz"), "train_miRNA"),
    ("GSE85307", RAW / "GSE85307", "GSE85307_miRNA_non-normalized.txt.gz",
     series_url("GSE85307", "suppl", "GSE85307_miRNA_non-normalized.txt.gz"), "train_miRNA_raw"),

    # --- ÖZGÜLLÜK 1: GSE149440 (Tarca 2021; PE / PPROM / sPTD / kontrol)
    ("GSE149440", RAW / "GSE149440", "GSE149440_series_matrix.txt.gz",
     series_url("GSE149440", "matrix", "GSE149440_series_matrix.txt.gz"), "specificity_1"),

    # --- ÖZGÜLLÜK 2: GSE59491 (All Our Babies; sPTL / PPROM / term)
    ("GSE59491", RAW / "GSE59491", "GSE59491_series_matrix.txt.gz",
     series_url("GSE59491", "matrix", "GSE59491_series_matrix.txt.gz"), "specificity_2"),

    # --- BÖLME + SON TEST: GSE192902 (Moufarrej 2022, plazma cfRNA)
    ("GSE192902", RAW / "GSE192902", "GSE192902_series_matrix.txt.gz",
     series_url("GSE192902", "matrix", "GSE192902_series_matrix.txt.gz"), "cfRNA_labels"),
    ("GSE192902", RAW / "GSE192902", "GSE192902_counts_Discovery_postQC.csv.gz",
     series_url("GSE192902", "suppl", "GSE192902_counts_Discovery_postQC.csv.gz"), "cfRNA_discovery"),
    ("GSE192902", RAW / "GSE192902", "GSE192902_counts_Validation1_postQC.csv.gz",
     series_url("GSE192902", "suppl", "GSE192902_counts_Validation1_postQC.csv.gz"), "cfRNA_validation1"),
    ("GSE192902", SEALED, "GSE192902_counts_Validation2_postQC.csv.gz",
     series_url("GSE192902", "suppl", "GSE192902_counts_Validation2_postQC.csv.gz"), "SEALED_final_test"),

    # --- ÖZGÜLLÜK 0: GSE142974 (VDAART, GSE85307 ile aynı kohort/platform/pencere; sPTB vs normal)
    ("GSE142974", RAW / "GSE142974", "GSE142974_series_matrix.txt.gz",
     series_url("GSE142974", "matrix", "GSE142974_series_matrix.txt.gz"), "specificity_0"),

    # --- ÖZGÜLLÜK 3: GSE108497 (Hong 2019; SLE + sağlıklı gebelik; PE / FGR-SGA-fetal ölüm)
    #     Series matrix yalnız anotasyon içerir; ifade verisi ek dosyada (sütunlar çip barkodu).
    ("GSE108497", RAW / "GSE108497", "GSE108497_series_matrix.txt.gz",
     series_url("GSE108497", "matrix", "GSE108497_series_matrix.txt.gz"), "specificity_3_labels"),
    ("GSE108497", RAW / "GSE108497", "GSE108497_normalized_data.txt.gz",
     series_url("GSE108497", "suppl", "GSE108497_normalized_data.txt.gz"), "specificity_3"),
    # Ana dosyada olmayan iki çip (9269325029_J, _K; ikisi de SLE_P_OC) ayrı dosyada.
    ("GSE108497", RAW / "GSE108497", "GSE108497_9269325029_J_and_K_normalized.txt.gz",
     series_url("GSE108497", "suppl", "GSE108497_9269325029_J_and_K_normalized.txt.gz"), "specificity_3_JK"),

    # --- İKİNCİ MÜHÜRLÜ TEST: GSE154377 (Del Vecchio 2021, plazma cfRNA; Moufarrej'in dış doğrulama kohortu)
    #     Etiketler (series matrix) açık; ifade verisi data/sealed/ altında, son testte bir kez değerlendirilir.
    ("GSE154377", RAW / "GSE154377", "GSE154377_series_matrix.txt.gz",
     series_url("GSE154377", "matrix", "GSE154377_series_matrix.txt.gz"), "sealed2_labels"),
    ("GSE154377", SEALED, "GSE154377_RAW.tar",
     series_url("GSE154377", "suppl", "GSE154377_RAW.tar"), "SEALED2_final_test"),

    # --- HAFTA REFERANSI: GSE149436 (KULLANILMIYOR: tamamı GSE149440'ın içinde; kayıt için duruyor)
    ("GSE149436", RAW / "GSE149436", "GSE149436_series_matrix.txt.gz",
     series_url("GSE149436", "matrix", "GSE149436_series_matrix.txt.gz"), "ga_reference"),

    # --- PLATFORM ANOTASYONLARI (prob -> gen sembolü, Adım 3 için)
    ("platforms", PLATFORMS, "GPL6244.annot.gz",
     f"{GEO}/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz", "annot_GSE85307_mRNA"),
    ("platforms", PLATFORMS, "GPL19162_table.txt",
     platform_table_url("GPL19162"), "annot_GSE85307_miRNA"),
    ("platforms", PLATFORMS, "GPL28460_table.txt",
     platform_table_url("GPL28460"), "annot_GSE149440"),
    ("platforms", PLATFORMS, "GPL18964_table.txt",
     platform_table_url("GPL18964"), "annot_GSE59491"),
    ("platforms", PLATFORMS, "GPL10558.annot.gz",
     f"{GEO}/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz", "annot_GSE108497"),

    # --- GEN REFERANSI: Entrez <-> sembol <-> Ensembl (GPL28460/GPL18964 yalnız Entrez,
    #     GSE192902 Ensembl kullanıyor). NCBI günlük günceller; sürüm = manifest tarihi + SHA-256.
    ("reference", ROOT / "data" / "reference", "Homo_sapiens.gene_info.gz",
     "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
     "gene_id_map"),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url, dest, retries=3):
    part = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "pe-specificity-study/0.1"})
            with urllib.request.urlopen(req, timeout=120) as r:
                expected = r.headers.get("Content-Length")
                expected = int(expected) if expected else None
                got = 0
                t0 = time.time()
                with open(part, "wb") as out:
                    while True:
                        buf = r.read(1 << 20)
                        if not buf:
                            break
                        out.write(buf)
                        got += len(buf)
                        if expected:
                            pct = 100 * got / expected
                            mbps = got / 1e6 / max(time.time() - t0, 1e-6)
                            print(f"\r    {got/1e6:7.1f} / {expected/1e6:.1f} MB ({pct:5.1f}%) {mbps:5.1f} MB/s",
                                  end="", flush=True)
                print()
            if expected is not None and got != expected:
                raise IOError(f"eksik indirme: {got} / {expected} bayt")
            if got == 0:
                raise IOError("boş yanıt")
            part.replace(dest)
            return
        except Exception as e:
            print(f"    deneme {attempt}/{retries} başarısız: {e}")
            if part.exists():
                part.unlink()
            if attempt == retries:
                raise
            time.sleep(5 * attempt)


def main():
    only = set(sys.argv[1:])
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    new_manifest = not MANIFEST.exists()
    with open(MANIFEST, "a", encoding="utf-8") as mf:
        if new_manifest:
            mf.write("utc_time\tgroup\trole\tfile\tbytes\tsha256\turl\n")
        for group, folder, fname, url, role in FILES:
            if only and group not in only:
                continue
            folder.mkdir(parents=True, exist_ok=True)
            dest = folder / fname
            if dest.exists():
                print(f"[var]   {dest.relative_to(ROOT)}")
                continue
            print(f"[indir] {dest.relative_to(ROOT)}")
            download(url, dest)
            digest = sha256(dest)
            size = dest.stat().st_size
            if role.startswith("SEALED"):  # SEALED ve SEALED2
                os.chmod(dest, stat.S_IREAD)  # salt-okunur
            ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            mf.write(f"{ts}\t{group}\t{role}\t{dest.relative_to(ROOT).as_posix()}\t{size}\t{digest}\t{url}\n")
            mf.flush()
            print(f"        {size/1e6:.2f} MB  sha256={digest[:16]}…")


if __name__ == "__main__":
    main()
