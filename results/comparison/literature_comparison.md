# Literatür karşılaştırma tablosu — AUC ve doğruluk (2026-10-22, mevcut durum)

Kısaltmalar: DD = dengeli doğruluk = (duyarlılık + özgüllük) / 2 (sınıf oranından bağımsız; karşılaştırma ölçüsü).
Acc = ham doğruluk (sınıf oranına bağlı; Val-2'de herkese "kontrol" demek 0.70 verir). Duy = duyarlılık, Özg = özgüllük.
"Bildirilen" sütunları makalenin kendi değerleridir (kaynaklar: PubMed özet / PMC tam metin; bkz.
`results/literature_specificity/reported_performance.tsv`). "Bizim yeniden testimiz" aynı imzanın bizim bağımsız
kohortlarımızdaki performansıdır (`scripts/10`, `scripts/11`, `scripts/13`). Eşikler yalnız eğitim verisinden (Youden).

## Tablo 1 — Yayımlanmış imzalar: bildirilen ve bağımsız yeniden test

| İmza | Bölme | Bildirilen AUC | Bildirilen Duy / Özg (DD) | Bildirilen doğrulama türü | Bizim yeniden testimiz: AUC [%95 GA] (kohort) | Bizim DD |
|---|---|---|---|---|---|---|
| Moufarrej 2022, 18 gen | plazma | eğitim 0.99; Val-1 0.71; Val-2 0.72; Del Vecchio 0.74 | Val-2: 0.56 / 0.69 (0.63) | aynı çalışma, ayrılmış + dış | **Val-2: 0.73 [0.60–0.84]**; Val-1: 0.66; Del Vecchio: 0.41 (8/8) | Val-2: 0.68 (Duy 0.81 / Özg 0.56) |
| Rasmussen 2022, 7 gen | plazma | 0.82 | Duy 0.75; PPV 0.32 | birini-dışarıda-bırak CV | plazma: test edilemedi (genler GSE192902'de yok); tam kan: 0.49 | — |
| Castillo-Marco 2025, EOPE-T1 (36) | plazma | 0.88 (dış 0.87) | 0.83 / 0.90 (0.87); dış 0.78 / 0.90 (0.84) | %30 ayrılmış; dış veride yeniden ayarlı | plazma: test edilemedi (kapsam < %50); tam kan: 0.40 | — |
| Elovitz 2025 | plazma | 0.88 | — | doğrulama alt grubu | Val-2: 0.35; Val-1: 0.47 | — |
| Weiner 2022, NAMPT + APOA1 | plazma | 0.89 | — | iç | Val-2: 0.41; Val-1: 0.43; Del Vecchio: 0.80 (8/8; GDM'de de 0.59) | — |
| Tarca 2021, 3 gen (11–17 hf) | tam kan | — | 0.83 / 0.94 (0.89) | seçimle aynı veri | bağımsız tam kan: 0.49–0.51 | — |
| Tarca 2021, 4 gen (22–28 hf) | tam kan | — | 0.85 / 0.92 (0.89) | seçimle aynı veri | bağımsız tam kan: 0.51 | — |
| Del Vecchio 2021 | plazma | 0.86 | — | iç (33 kadın) | test edilemedi (kapsam < %50) | — |
| Gunaydin 2025, 12 gen | tam kan | tek gen 0.74–0.87 | — | iç | test edilemedi (kapsam < %50) | — |
| Hong 2019, nötrofil modülü | tam kan | bildirilmemiş | — | — | bağımsız tam kan: 0.48 | — |
| Mirzakhani 2016, 348 gen | tam kan | bildirilmemiş | — | — | bağımsız tam kan: 0.50 | — |
| Hedley 2025, 7 miRNA | serum | 0.88 | — | iç | test edilemedi (eşleşen miRNA verisi yok) | — |
| Timofeeva 2023, miRNA | serum | — | 1.00 / 0.84–1.00 (0.92–1.00) | iç (67 hasta) | test edilemedi | — |

Özet: 26 imza / 228 bağımsız test; FDR < 0.20'de anlamlı hücre: 0. Bağımsız veride tekrarlanan tek imza Moufarrej 2022.

## Tablo 2 — Bu çalışmanın modelleri

| Model | Bölme | Geliştirme (iç içe CV) AUC | İç doğrulama (aynı merkez) AUC | DIŞ test AUC [%95 GA] | Dış DD (Duy / Özg) | Dış Acc |
|---|---|---|---|---|---|---|
| XAI-MOGAFS (özgüllük amaçlı) | tam kan | 0.50 [0.39–0.60] | — | — (tam kanda sinyal yok) | — | — |
| XAI-MOGAFS V0 | plazma | 0.755 | Val-1: 0.852 | Val-2: 0.436 | — | — |
| **XAI-MOGAFS V7 — önceden kayıtlı BİRİNCİL** | plazma | 0.786 | Val-1: 0.821 | **Val-2: 0.418 [0.29–0.55]** | 0.44 (0.19 / 0.69) | 0.54 |
| XAI-MOGAFS-SI (keşifsel, Val-2 sonrası) | plazma | — | — | Val-2 kurumları: 0.78 [0.66–0.86] | 0.72 (0.69 / 0.74) | 0.72 |
| 〃 | plazma | — | — | Stanford: 0.40 [0.26–0.53] | 0.42 (0.38 / 0.46) | 0.44 |
| 〃 | plazma | — | — | Del Vecchio: 0.27 [0.05–0.56] (8/8) | 0.31 | 0.31 |

Notlar:
1. Val-1, Discovery ile aynı Stanford kohortundan (Moufarrej 2022 Yöntemler) — "iç doğrulama" bağımsız değildir.
2. XAI-MOGAFS-SI, Val-2 sonucu görüldükten sonra tasarlanmıştır; üç merkez-dışarıda katı birlikte raporlanmalıdır
   (ortalama ≈ 0.48). Val-2 katındaki 0.78 yeni bağımsız kohortta doğrulanmadan üstünlük kanıtı değildir.
3. Bildirilen değerlerin çoğu iç doğrulamadır; bağımsız yeniden testte 0.82–0.94'ten 0.35–0.51'e iner.
