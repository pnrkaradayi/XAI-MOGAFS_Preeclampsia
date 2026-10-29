# Plazma kolu — yöntem varyantları kaydı

Kullanıcı isteği (2026-10-22): son testten önce yöntemin geliştirilmesi. Bu dosya varyantlar ÇALIŞTIRILMADAN
önce yazılmıştır. Denenen her varyant, sonucu ne olursa olsun buraya eklenir; makalenin ek yöntemlerinde
"denenen varyantlar" olarak raporlanır. Validation-2 bu süreçte AÇILMAZ.

## Varyantlar (hepsi Discovery iç içe CV 5 kat x 3 tekrar + Val-1 iç doğrulama + Disc+Val-1 kilit)

| Kod | Gen evreni | Bileşim özellikleri | Kararlılık amacı (f5) | Hafta eşlik eden değişkeni |
|---|---|---|---|---|
| V0 | tüm genler | – | – | – |  (kilit #1, zaten çalıştı)
| V0pc | protein kodlayan | – | – | – |  (kilit #3, zaten çalıştı)
| V1 | tüm genler | ✓ | – | – |
| V2 | tüm genler | ✓ | ✓ | – |
| V3 | tüm genler | ✓ | ✓ | ✓ |
| V4 | protein kodlayan | ✓ | ✓ | ✓ |

Bileşim özellikleri: ham sayımlardan örnek başına kütüphane payları (log2(% + 0.01)): RN7SL*, RN7SK,
RNA5-8S*/RNA5S*, RNY*, SNOR*/SCARNA*, RNU*, MT-*. Aday özellik olarak genlere eklenir; panel boyutuna (≤5) sayılır.

## ÖN-YAZILI SEÇİM KURALI (sonuçlar görülmeden)
1. Aday: MOGAFS satırı; Val-1 AUC ≥ 0.70 olmalı.
2. Adaylar arasında en yüksek Discovery iç içe CV AUC'si (3 tekrarın ortalaması).
3. 0.01 içinde eşitlikte daha az ek bileşen içeren (daha basit) varyant.
4. Seçilen varyantın kilidi Val-2 için BİRİNCİL olur; V0 (kilit #1) ikincil olarak kalır. Diğer tüm
   varyantların kilitleri de ikincil olarak Val-2'de raporlanır (seçici raporlama yok).

## Sonuçlar

MOGAFS satırları (LASSO5 / TTEST5 aynı özellik uzayında, karşılaştırma için):

| Kod | Discovery iç içe CV AUC (3 tekrar ort.; tekrarlar) | Val-1 AUC | LASSO5 CV | TTEST5 CV | Kilit SHA-256 (ilk 16) |
|---|---|---|---|---|---|
| V0 | 0.755 (0.711, 0.732, 0.820) | 0.852 | 0.710 | 0.827 | f0546297562c972d |
| V0check (V0 yeniden) | 0.755 — paneller V0 ile birebir aynı | 0.852 | 0.710 | 0.827 | (zaman damgası farkı) |
| V0pc | 0.545 (0.573, 0.414, 0.647) | 0.574 | 0.501 | 0.492 | b140fe7507fefec6 |
| V1 | 0.741 (0.707, 0.796, 0.719) | 0.864 | 0.710 | 0.827 | c07d7a18c1e1e823 |
| V2 | 0.741 (0.757, 0.727, 0.741) | 0.790 | 0.710 | 0.827 | 17639114ebc5ffc5 |
| V3 | 0.741 (0.770, 0.721, 0.732) | 0.833 | 0.710 | 0.822 | 71b9aedcab837909 |
| V4 | 0.671 (0.670, 0.630, 0.714) | 0.840 | 0.508 | 0.496 | 47b38b5df7ede641 |

Gözlem: tüm-gen evreninde (tek baskın RN7SL bloğu) TTEST5 > MOGAFS; protein kodlayan evrende (dağınık,
zayıf sinyal) MOGAFS+kararlılık+bileşim (V4) LASSO5/TTEST5'ten ~0.17 yüksek.

## V1–V4 SONUÇLARI GÖRÜLDÜKTEN SONRA EKLENEN varyantlar (açıkça belirtilir)

Gerekçe: Pareto cephesinden seçim kuralı (en iyi iç-CV AUC'nin tol=0.03 içindeki EN KÜÇÜK panel) tüm-gen
evreninde AUC'den feragat ediyor (3 özellik). tol=0 → cephedeki en yüksek iç-CV AUC'li çözüm.

| Kod | Tanım |
|---|---|
| V5 | V2 + tol = 0 (tüm genler, bileşim, kararlılık) |
| V6 | V4 + tol = 0 (protein kodlayan, bileşim, kararlılık, hafta) |

Seçim kuralı DEĞİŞMEDİ (yukarıdaki 1–4); V5/V6 da aynı kurala tabidir.

| Kod | Discovery iç içe CV AUC (tekrarlar) | Val-1 AUC | LASSO5 CV | TTEST5 CV | Kilit SHA-256 (ilk 16) |
|---|---|---|---|---|---|
| V5 | 0.740 (0.724, 0.757, 0.739) | 0.759 | 0.710 | 0.827 | c1812ef68ef05df2 |
| V6 | 0.623 (0.656, 0.543, 0.670) | 0.802 | 0.508 | 0.496 | 619c1310d97e281f |

## KARAR (kural uygulandı, 2026-10-22)
Val-1 ≥ 0.70 koşulunu V0, V1, V2, V3, V4, V5, V6 sağlıyor (V0pc sağlamıyor: 0.574). En yüksek Discovery
CV AUC'si V0 (0.755); sonraki V1/V2/V3 (0.741) farkı 0.014 > 0.01 → **BİRİNCİL = V0 (kilit #1)**.
Varyantların hiçbiri tüm-gen evreninde V0'ı geçmedi; tol=0 (V5, V6) iç-CV'ye aşırı uyum nedeniyle düşürdü.
Yöntem bulgusu: zayıf/dağınık sinyalde (protein kodlayan evren) kararlılık + bileşim eklentili MOGAFS (V4,
0.671) klasik yöntemleri (0.50) belirgin geçiyor; tek baskın blokta (tüm genler) TTEST5 en iyisi (0.827).
Tüm varyant kilitleri Val-2'de ikincil olarak raporlanır.

## Literatür karşılaştırması sonrası EKLENEN varyantlar (kullanıcı isteği, 2026-10-22; çalıştırılmadan önce yazıldı)

Gerekçe: Moufarrej 2022 eğitimde kadın başına birden çok erken örnek kullanmış ve örnek düzeyinde raporlamış.
Aynı veri kullanımını sızıntısız uygulamak: eğitimde her kadının TÜM ≤20 hf örnekleri, tüm CV bölmeleri
KADIN düzeyinde (StratifiedGroupKFold, iç ve dış). Test her zaman eğitimde hiç görülmemiş kadınlarda.

| Kod | Tanım |
|---|---|
| V7 | tüm genler + bileşim, ≤5 özellik (V1'in tekrarlı-örnek eğitimli hâli) |
| V8 | V7 ile aynı, ≤18 özellik (Moufarrej'in 18 genine eşit boyut; İKİNCİL, ≤5 kuralı dışında olduğu açıkça belirtilir) |

Raporlama: (a) kadın düzeyi (kadının en erken örneği; V0 ile aynı) ve (b) örnek düzeyi (tüm ≤20 hf örnekleri;
Moufarrej'in raporlama birimi). Eğitim (apparent) AUC'si de "eğitim" etiketiyle raporlanır.
Seçim kuralı DEĞİŞMEZ (kadın düzeyi Discovery iç içe CV AUC'si; Val-1 ≥ 0.70). V8 birincil olamaz (≤5 kuralı).
Not: GEO yalnız ≤12 / 13–20 hf dilimlerini verdiği için Moufarrej'in "≤16 hf" penceresi birebir kurulamaz.

### V7 / V8 sonuçları

| Kod | Discovery CV AUC — kadın düzeyi (tekrarlar) | örnek düzeyi | Val-1 kadın / örnek | Eğitim (apparent) | LASSO / TTEST CV (kadın) | Kilit SHA-256 (ilk 16) |
|---|---|---|---|---|---|---|
| V7 (≤5) | 0.786 (0.777, 0.774, 0.809) | 0.793 | 0.821 / 0.777 | 0.974 | 0.754 / 0.840 | 6dca6bf65908ec8c |
| V8 (≤18, ikincil) | 0.778 (0.836, 0.729, 0.770) | 0.786 | 0.753 / 0.735 | 0.981 | 0.796 / 0.770 | eda7317beacda60f |

Aynı Val-1 örneklerinde (Discovery'de kurulan modeller): V0 paneli kadın 0.877 / örnek 0.826;
Moufarrej 18 gen (16/18 mevcut, yayımlanmış log OR) kadın 0.648 / örnek 0.605 (yayımlanmış Val-1: 0.71).

## KARAR GÜNCELLEMESİ (aynı kural)
Val-1 ≥ 0.70 olan MOGAFS satırları arasında en yüksek Discovery CV (kadın düzeyi) AUC'si V7 (0.786 > V0 0.755,
fark 0.031 > 0.01) → **BİRİNCİL = V7** (`variants/V7/locked_models.json`). V0 ikincil olarak raporlanır.
V7, V1–V6 sonuçlarından sonra eklenmiştir (yukarıda belirtildi); bu sıralı varyant denemesi makalede açıkça yazılır
ve Discovery CV tahmininde seçim iyimserliği olabileceği sınırlılık olarak belirtilir. Kesin performans Val-2'dir.

## SONRADAN (POST-HOC) — Val-2 açıldıktan sonra (2026-10-22, kullanıcı isteği)
Val-2 birincil sonucu (V7 0.418) DEĞİŞMEZ ve birincil olarak raporlanır. Başarısızlık nedenine (Stanford'a özgü RN7SL
izi) yanıt olarak XAI-MOGAFS-SI (merkez-tutarlı ön-filtre, en kötü merkez AUC amacı, merkez-değişmez SHAP) tasarlandı
ve merkez-dışarıda-bırakma ile değerlendirildi (`scripts/12_site_invariant.py`, `results/site_invariant/`).
Bu analizler keşifseldir: Val-2 ve Del Vecchio sonuçları görüldükten sonra tasarlanmıştır; bağımsız doğrulama için
yeni bir kohort gerekir.

### XAI-MOGAFS-SI merkez-dışarıda-bırakma sonuçları (keşifsel; tüm katlar raporlanır)

| Dışarıda bırakılan | SI-MOGAFS | MOGAFS | META-TTEST5 | LASSO5 | TTEST5 | Moufarrej18 |
|---|---|---|---|---|---|---|
| S2 Val-2 kurumları (26/61) | 0.78 [0.66–0.86] | 0.58 | 0.49 | 0.53 | 0.61 | 0.73 |
| S1 Stanford (29/72) | 0.40 [0.26–0.53] | 0.47 | 0.65 | 0.61 | 0.58 | 0.76 (eğitim kohortu, bağımsız değil) |
| S3 Del Vecchio (8/8) | 0.27 [0.05–0.56] | 0.45 | 0.19 | 0.45 | 0.28 | 0.41 |

Yorum: SI-MOGAFS Val-2 katında iyi, diğer iki katta başarısız; üç kat ortalaması ≈ 0.48. Tutarlı iyileşme YOK.
Muhtemel neden: en-kötü-merkez amacı her katta 8 vakalık S3'e dayanıyor. Val-2 katındaki 0.78 bağımsız yeni veride
doğrulanmadan iddia edilemez (yeni veri talebi: docs/veri_talebi_epostalari.md). Daha fazla varyant DENENMEDİ
(aynı üç merkezde çoklu deneme şans bulgusu üretir).
