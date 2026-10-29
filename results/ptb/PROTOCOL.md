# Kendiliğinden erken doğum (sPTB) paneli — ön-kayıt protokolü (2026-10-23, çalıştırılmadan önce)

Amaç: XAI-MOGAFS'ın klasik özellik seçim yöntemlerine göre bağımsız bir kohorta daha iyi genelleyip
genellemediğini, iki bağımsız tam kan kohortunda temiz bir geliştirme → tek seferlik dış test tasarımıyla sınamak.

- Hedef: kendiliğinden erken doğum bileşiği (sPTD/sPTL + PPROM) vs term kontrol. PE kadınları dışlanır.
- Pencere (eşleştirilmiş): örnekleme 17–23 hf (GSE59491 zaman noktası 1 ile aynı). GSE149440'ta kadının bu
  penceredeki EN ERKEN örneği.
- Geliştirme: GSE149440 (HTA 2.0). Dış test: GSE59491 (HuGene 2.1) — bu kohortta hiçbir sPTB modeli eğitilmedi /
  seçilmedi; sPTB etiketleri yalnız yayımlanmış PE imzalarının değerlendirilmesinde kullanıldı (scripts/10).
- Gen evreni: 5 tam kan kohortunun ortak 17.621 geni; her kohort kendi içinde z-ölçeklenir.
- Yöntemler (aynı ön-işleme, ≤5 gen): XAI-MOGAFS (f1, f3, f4 + kararlılık f5; kararlılık öncelikli seçim, tol 0.03),
  LASSO5, TTEST5; ayrıca tüm genlerle ridge (üst sınır referansı).
- Geliştirme içi tahmin: iç içe CV (dış 5 kat x 3 tekrar).
- Son modeller: GSE149440'ın tamamında eğitilip KİLİTLENİR (SHA-256), sonra GSE59491'e BİR KEZ uygulanır.
- BİRİNCİL sonuç: XAI-MOGAFS'ın GSE59491 AUC'si + %95 GA. İkincil: XAI-MOGAFS vs LASSO5 / TTEST5 DeLong testi;
  alt tipler (sPTL, PPROM) ayrı AUC.
- Tüm sonuçlar, iyi ya da kötü, raporlanır.

## SONUÇ (2026-10-23; protokole uygun, bir kez)
- Geliştirme (GSE149440, iç içe CV): XAI-MOGAFS 0.404, LASSO5 0.394, TTEST5 0.369; ridge (tüm genler) 0.448 → sinyal yok.
- Dış test (GSE59491): XAI-MOGAFS 0.571 [0.48–0.66], LASSO5 0.568, TTEST5 0.596; DeLong farkları anlamsız (p 0.95, 0.59).
- Yorum: erken (17–23 hf) tam kan transkriptomu kendiliğinden erken doğumu da öngörmüyor; PE bulgusuyla tutarlı.
  Birincil hipotez (XAI-MOGAFS daha iyi genelleşir) DESTEKLENMEDİ. Kilit SHA-256 d2ef61bc3aeabac9…
