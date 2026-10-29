# Plazma kolu — kilit kaydı ve son test (Validation-2) ön-kayıt planı

Bu dosya, mühürlü GSE192902 Validation-2 açılmadan ÖNCE yazılmıştır (2026-10-22).
Validation-2 sayım dosyası (`data/sealed/GSE192902_counts_Validation2_postQC.csv.gz`) bu tarihe kadar
hiçbir betik tarafından okunmamıştır; yalnızca SHA-256'sı manifest ile karşılaştırılmıştır.

## Kronoloji

| # | Olay | Dosya / özet |
|---|---|---|
| 1 | Tüm-gen evreninde (7.025 gen) MOGAFS, LASSO5 ve TTEST5 paneli kuruldu ve kilitlendi; Discovery iç içe CV AUC 0,71–0,83, Val-1 AUC 0,85–0,91 | `all_genes/locked_models.json`, SHA-256 `f0546297562c972d33120aafdcaabab03a6edc47ff9231418ec4b3edd69270fb` |
| 2 | Tanılama (`scripts/09b_plasma_diagnostics.py`): seçilen genler çoğunlukla RN7SL (7SL/SRP RNA) psödogenleri. RN7SL kütüphane payı PE'de ~2 kat (Discovery AUC 0,76; Val-1 0,76); hafta ile ilişkisiz; trombosit mRNA'sıyla ilişkisiz (ρ ≈ 0) | `plasma_diagnostics_report.txt` |
| 3 | Yalnız protein kodlayan evrende (6.707 gen) aynı prosedür kilitlendi; ≤5 genlik paneller iç içe CV'de tutmadı (AUC 0,49–0,55) | `protein_coding/locked_models.json`, SHA-256 `b140fe7507fefec65961f780e6e2e1258c0f78ee75dbb1fa35207fef224e9532` |
| 4 | Kullanıcı isteğiyle yöntem varyantları V1–V6 (bileşim özellikleri, kararlılık amacı f5, hafta eşlik değişkeni, tol=0) geliştirme verisinde denendi; önceden yazılmış kurala göre birincil DEĞİŞMEDİ (V0) | `VARIANTS_LOG.md`; kilitler `variants/<kod>/locked_models.json` |
| 5 | Literatür karşılaştırması sonrası (kullanıcı isteği) tekrarlı-örnek eğitimli V7 (≤5) ve V8 (≤18) denendi; kurala göre BİRİNCİL = V7 | `variants/V7/locked_models.json` SHA-256 `6dca6bf65908ec8c…`; `variants/V8/` |

## Validation-2 değerlendirme planı (bir kez açılır)

Uygulama: Val-2'de ≤20 hf, kadın başına en erken örnek (seçim Adım 2'de anotasyondan yapıldı:
26 PE / 61 kontrol). log2(CPM+1); her gen Val-2 içinde (etiketsiz) z-ölçeklenir; kilitli katsayılar
DEĞİŞTİRİLMEDEN uygulanır. Kilitli bir gen yoksa z = 0.

| Rol | Model | Kaynak |
|---|---|---|
| **BİRİNCİL** | MOGAFS V7: tüm genler + bileşim, tekrarlı-örnek eğitimi, ≤5 özellik; kadın düzeyi (kadının en erken ≤20 hf örneği). Val-2 içinde z-ölçekleme TÜM ≤20 hf Val-2 örnekleri üzerinden | `variants/V7/locked_models.json` |
| İkincil (önceki birincil) | MOGAFS V0, tüm-gen evreni (kilit #1) | `all_genes/locked_models.json` |
| İkincil | V7 örnek düzeyi (Moufarrej'in raporlama birimi); V8 (≤18) kadın ve örnek düzeyi | `variants/V7`, `variants/V8` |
| Literatür rakibi | Moufarrej 2022 18 gen, yayımlanmış log OR ağırlıklarıyla (yayımlanmış Val-2 AUC 0.72) | `results/literature/published_signatures.tsv` |
| Birincil hipotez | RN7SL payı (RN7SL* satırlarının ham sayım toplamı / kütüphane toplamı); yön önceden belirli: yüksek = PE | ham sayımlar |
| İkincil | LASSO5 ve TTEST5, tüm-gen evreni (kilit #1) | `all_genes/locked_models.json` |
| İkincil | MOGAFS, LASSO5, TTEST5, protein kodlayan evren (kilit #3) | `protein_coding/locked_models.json` |
| İkincil | Varyant kilitleri V1–V6 (her biri MOGAFS, LASSO5, TTEST5); V4 = protein kodlayan evrende yöntem karşılaştırmasının ana satırı | `variants/V1..V6/locked_models.json` |
| İkincil | Plazma kapsamı ≥ 0.5 olan yayımlanmış imzalar, `scripts/10` kuralıyla (yön ağırlıklı ortalama z): Moufarrej 2022 (işaret ve log OR), Weiner 2022, Weiner 2026, Elovitz 2025 top2, Tarca 2021 LDA4, Myhrer 2024, Lekva 2021 top20 | `results/literature/published_signatures.tsv` (katalog tamamlandı 2026-10-22) |

Raporlanacak: AUC + tabakalı bootstrap %95 GA (2000), hafta dilimine göre AUC (≤12 / 13–20),
hafta-düzeltmeli OR. Birincil sonuç tek sayıdır: MOGAFS V7 panelinin Val-2 kadın düzeyi AUC'si ve %95 GA'sı (güncellendi; gerekçe VARIANTS_LOG.md).

Gerekçe (birincil seçimi): Birincil model, VARIANTS_LOG.md'de varyantlar çalıştırılmadan önce yazılmış seçim kuralıyla
belirlenmiştir (en yüksek Discovery kadın düzeyi iç içe CV AUC'si, Val-1 ≥ 0.70). Kronolojik sıra (V0 → V1–V6 → V7/V8)
ve varyant sayısı makalede raporlanır. RN7SL sinyalinin kökeni (teknik / biyolojik) bu testle kesinleşmez; sınırlılık
olarak tartışılır.

## İkinci mühürlü test: GSE154377 (Del Vecchio 2021) — ön-kayıt (2026-10-22, açılmadan önce)

Kullanıcı kararıyla eklendi. Moufarrej 2022 bu kohortta AUC 0.74 bildirmiştir (kendi örnek tanımıyla).
İfade verisi `data/sealed/GSE154377_RAW.tar` (salt-okunur, SHA-256 manifestte). Yalnız etiketler
(series matrix) ve bir sayım dosyasının BİÇİMİ incelenmiştir (STAR ReadsPerGene, 57.740 Ensembl satırı).

- Örnek seçimi (etiketlerden): kadın başına 1. trimester örneği. Doğum / kordon kanı / gebe olmayan dışlanır.
  Gruplar: Komplikasyonsuz (kontrol, 8), PE/gestasyonel hipertansiyon (8), gestasyonel diyabet (7),
  kronik hipertansiyon (2; yalnız betimleyici).
- Karşılaştırmalar: PE/GH vs kontrol (birincil, bu kohort için); GDM vs kontrol (plasental olmayan komplikasyonda
  "alarm" kontrolü). n çok küçük -> yalnız AUC + %95 GA; hipotez testi değil, tutarlılık kanıtı.
- Ön-işleme: STAR sütunu = tüm örneklerde toplam atanmış okuması en yüksek iplik sütunu (etiketsiz kural);
  Ensembl (eşsiz) -> Entrez; aynı Entrez'de toplam; log2(CPM+1); bileşim payları gen sembolünden; kohort içi z.
  1. trimester -> hafta dilimi ≤12 (eşlik değişkenli modeller için 0).
- Modeller: Val-2 ile AYNI liste (birincil V7).

## Kör ön-işleme ve tek açılış protokolü
`scripts/11_final_test.py`: (1) mühürlü dosyaların SHA-256'sını manifestle doğrular; (2) Val-2 ve Del Vecchio
ifade verisini ETİKETSİZ işler ve özellik matrislerini kaydeder; (3) yalnızca sonra etiketleri birleştirip tüm
kayıtlı modelleri bir kez değerlendirir; (4) açılış zamanını `results/final_test/OPENED.txt`'e yazar. Yeniden
çalıştırma yalnız teknik hata düzeltmesi için ve gerekçesi kaydedilerek yapılabilir (`--rerun-reason`).

## Açılış günlüğü
- 2026-10-22, 1. çalıştırma: mühür ve kilit SHA-256 doğrulamaları geçti; Val-2 etiketsiz ön-işleme tamamlandı;
  Del Vecchio etiketsiz ön-işlemesinde teknik hata (sayım dosyalarında 3 yerine 1 sayım sütunu; KeyError) ile durdu.
  Etiketler birleştirilmeden ve hiçbir performans hesaplanmadan durduğu için (`results/final_test/OPENED.txt`
  oluşmadı) yalnız ayrıştırıcı düzeltildi; başka hiçbir şey değiştirilmedi.
- 2026-10-22, 2. çalıştırma: ayrıştırıcı düzeltmesiyle tamamlandı. Etiketler bir kez birleştirildi (`results/final_test/OPENED.txt`).

## SONUÇ (bir kez açıldı; değiştirilmez)
- **BİRİNCİL:** V7 MOGAFS, Val-2 kadın düzeyi AUC **0.418** [%95 GA 0.29–0.55] (26 PE / 61 kontrol). Birincil hipotez
  DESTEKLENMEDİ.
- Moufarrej 2022 18 gen (aynı kod, yayımlanmış yön/OR): AUC **0.734** [0.61–0.85] (yayımlanmış 0.72 ile uyumlu →
  işleme hattı doğru). DeLong V7 vs Moufarrej: fark −0.31, p = 0.004.
- RN7SL payı (birincil hipotez): Val-2 AUC 0.429 → tekrarlamadı. Del Vecchio: PE/GH 0.78, GDM 0.79 (n=8/8, 7/8) → özgül değil.
- Yorum: Moufarrej 2022 Yöntemler'e göre Discovery ve Validation-1 AYNI Stanford kohortundan dizileme partilerine
  bölünerek oluşturulmuş; Validation-2 farklı kurumlardan. RN7SL sinyali Stanford kohortuna özgü (ön-analitik /
  çalışma düzeyi) bir iz olarak yorumlanır; Val-1 "bağımsız" bir test olmadığı için iç doğrulama iyimser çıktı.
- İkincil (çoklu test; keşifsel): Val-2'de en iyi kendi modellerimiz protein kodlayan evrendeydi (LASSO5 0.62–0.63,
  MOGAFS 0.60); yapısal küçük RNA içeren tüm modeller ~0.42–0.53.
Tam tablo: `results/final_test/final_results.tsv`, `results/final_test/final_report.txt`.
