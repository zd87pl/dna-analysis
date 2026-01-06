# 🧬 Helixight - Otwarte Narzędzie do Analizy Genetycznej

<p align="center">
  <strong>Od pliku BAM do praktycznych wniosków</strong><br>
  Analizuj swoje dane z sekwencjonowania genomu lokalnie, prywatnie i za darmo.
</p>

<p align="center">
  <a href="#-ważne-zastrzeżenie">⚠️ Zastrzeżenie</a> •
  <a href="#-funkcje">Funkcje</a> •
  <a href="#-szybki-start">Szybki Start</a> •
  <a href="#-analizy">Analizy</a> •
  <a href="#-współpraca">Współpraca</a>
</p>

---

## ⚠️ WAŻNE ZASTRZEŻENIE

> **🚨 TO OPROGRAMOWANIE NIE JEST PRZEZNACZONE DO UŻYTKU KLINICZNEGO ANI MEDYCZNEGO 🚨**
>
> Helixight jest udostępniany **wyłącznie do celów edukacyjnych, badawczych i rozrywkowych**.
>
> - ❌ **NIE JEST** urządzeniem medycznym
> - ❌ **NIE JEST** klinicznie zwalidowany
> - ❌ **NIE JEST** zamiennikiem profesjonalnego poradnictwa genetycznego
> - ❌ **NIE JEST** przeznaczony do diagnozowania, leczenia lub zapobiegania jakimkolwiek chorobom
> - ❌ **NIE JEST** zatwierdzony przez FDA, EMA ani żaden organ regulacyjny
>
> **Wyniki produkowane przez to oprogramowanie:**
> - Mogą zawierać błędy, nieścisłości lub fałszywie pozytywne/negatywne wyniki
> - NIE powinny być używane do podejmowania jakichkolwiek decyzji medycznych lub zdrowotnych
> - NIE uwzględniają wszystkich czynników genetycznych, wpływów środowiskowych ani interakcji gen-gen
> - Opierają się na uproszczonych interpretacjach złożonej nauki genetycznej
>
> **Jeśli odkryjesz niepokojące wyniki:**
> - NIE panikuj - to nie są wyniki klasy klinicznej
> - Skonsultuj się z certyfikowanym doradcą genetycznym lub genetykiem medycznym
> - Wykonaj odpowiednie testy kliniczne w akredytowanych laboratoriach
>
> **Używając tego oprogramowania, potwierdzasz że:**
> - Rozumiesz, że jest to tylko dla zabawy i nauki
> - Nie będziesz podejmować żadnych decyzji zdrowotnych na podstawie tych wyników
> - Autorzy nie ponoszą odpowiedzialności za jakiekolwiek działania podjęte na podstawie wyników
>
> **Predyspozycja genetyczna ≠ przeznaczenie.** Większość cech i chorób jest kształtowana przez styl życia, środowisko i złożone interakcje między setkami genów.

---

## 🌟 Funkcje

- **🔒 100% Prywatne** - Wszystkie analizy działają lokalnie na Twoim komputerze. Twoje dane genetyczne nigdy go nie opuszczają.
- **💰 Darmowe i Open Source** - Licencja MIT, bez ukrytych kosztów, bez subskrypcji
- **🖥️ Interaktywne CLI** - Łatwy w użyciu interfejs z menu
- **🌍 Dwujęzyczny** - Wsparcie dla języka polskiego i angielskiego
- **📊 Kompleksowy** - Analiza 500+ wariantów genetycznych w 16 kategoriach
- **🎓 Edukacyjny** - Poznaj swoją genetykę w zabawny, przystępny sposób

---

## 📋 Co Możesz Zbadać?

| Kategoria | Opis | Przykłady |
|-----------|------|-----------|
| 🏃 **Genetyka Sportowa** | Typ włókien mięśniowych, potencjał wytrzymałościowy | ACTN3, ACE, PPARGC1A |
| 🏊 **Predyspozycje Triatlonowe** | VO2max, metabolizm tłuszczów, ryzyko kontuzji | Tendencje Sprint vs Ironman |
| 💊 **Spersonalizowane Wnioski** | Styl treningu, suplementacja, żywienie | Na podstawie MTHFR, COMT, CYP1A2 |
| 🧬 **Kompleksowa Analiza** | 500+ wariantów w 16 domenach | Farmakogenomika, cechy |
| 🔬 **Eksploracja Wariantów** | Znane istotne warianty | Eksploracja badawcza |
| 📊 **Wyniki Poligeniczne** | Połączone efekty wariantów | Edukacyjna eksploracja ryzyka |
| 🎲 **Zabawna Genetyka** | Cechy i wskazówki pochodzenia | Kolor oczu, kofeina, kolendra |

---

## 🚀 Szybki Start

```bash
# Sklonuj repozytorium
git clone https://github.com/helixight/helixight-oss.git
cd helixight-oss

# Zainstaluj zależności (Ubuntu/Debian)
sudo apt install bcftools samtools tabix wget

# Lub użyj instalatora
./install.sh

# Uruchom interaktywne menu
./helixight.sh
```

---

## 📦 Wymagania

| Narzędzie | Cel | Instalacja |
|-----------|-----|------------|
| bcftools | Manipulacja VCF | `apt install bcftools` |
| samtools | Przetwarzanie BAM/SAM | `apt install samtools` |
| tabix | Indeksowanie VCF | `apt install tabix` |
| wget | Pobieranie plików | `apt install wget` |

### Instalacja jedną komendą

**Ubuntu/Debian:**
```bash
sudo apt update && sudo apt install -y bcftools samtools tabix wget gzip curl
```

**macOS (Homebrew):**
```bash
brew install bcftools samtools htslib wget
```

---

## 📖 Użycie

### Tryb Interaktywny (Zalecany)

```bash
./helixight.sh
```

Uruchamia interaktywne menu, gdzie możesz:
1. ✅ Zainstalować zależności
2. ✅ Pobrać genom referencyjny (jeśli masz pliki BAM)
3. ✅ Stworzyć VCF z plików BAM
4. ✅ Uruchomić różne eksploracje genetyczne
5. ✅ Wygenerować raporty podsumowujące

### Bezpośrednie Uruchamianie Skryptów

```bash
# Eksploracja genetyki sportowej
./scripts/athletic_genetics.sh twoje_warianty.vcf.gz

# Analiza predyspozycji triatlonowych
./scripts/triathlon_genetics.sh twoje_warianty.vcf.gz

# Spersonalizowane wnioski
./scripts/personalized_protocol.sh twoje_warianty.vcf.gz

# Kompleksowa eksploracja (500+ wariantów)
./scripts/mega_analysis.sh twoje_warianty.vcf.gz

# Zabawne cechy
./scripts/fun_genetics.sh twoje_warianty.vcf.gz
```

---

## 🧬 Tworzenie VCF z BAM

Jeśli masz plik BAM z sekwencjonowania całego genomu:

### Opcja 1: Interaktywne Menu
```bash
./helixight.sh
# Wybierz opcję 4: Stwórz VCF z pliku BAM
```

### Opcja 2: Ręcznie przez bcftools
```bash
bcftools mpileup -Ou -f GRCh38.fa twoja_probka.bam | \
bcftools call -mv -Oz -o twoje_warianty.vcf.gz

bcftools index -t twoje_warianty.vcf.gz
```

---

## 📊 Dostępne Analizy

### 🏃 Genetyka Sportowa
Eksploruj geny związane z wydajnością sportową:
- **ACTN3** - Tendencje typu włókien mięśniowych
- **ACE** - Markery wydajności sercowo-naczyniowej
- **PPARGC1A** - Biogeneza mitochondriów
- **IL6** - Markery regeneracji i zapalenia
- **COMT** - Odpowiedź na stres (Warrior vs Worrier)

### 🏊 Genetyka Triatlonowa
Kompleksowa eksploracja zorientowana na triathlon:
- Markery potencjału VO2max
- Zdolność oksydacji tłuszczów
- Genetyka progu mleczanowego
- Markery podatności na kontuzje
- Wskaźniki odporności psychicznej
- Punktacja per dyscyplina (Pływanie/Rower/Bieg)

### 💊 Spersonalizowane Wnioski
Edukacyjne rekomendacje oparte na genetyce:
- Sugestie stylu treningowego
- Rozważania suplementacyjne (MTHFR, VDR, itp.)
- Tendencje żywieniowe
- Wnioski dotyczące stylu życia

### 🧬 Mega Analiza
Kompleksowa eksploracja 500+ wariantów:
- Markery rodowe (DNA neandertalskie)
- Eksploracja statusu nosicielstwa
- Chronotyp (sowa vs skowronek)
- Genetyka zmysłów (smak, węch)
- Przewidywanie cech

---

## 📁 Struktura Projektu

```
helixight-oss/
├── helixight.sh          # Główny interaktywny launcher
├── install.sh            # Szybki instalator
├── README.md             # Dokumentacja angielska
├── README_PL.md          # Dokumentacja polska
├── LICENSE               # Licencja MIT
├── scripts/              # Skrypty analityczne
│   ├── athletic_genetics.sh
│   ├── triathlon_genetics.sh
│   ├── personalized_protocol.sh
│   ├── mega_analysis.sh
│   ├── fun_genetics.sh
│   └── ... (20+ skryptów)
└── docs/
    └── variant_database.md
```

---

## 🔬 Nota Naukowa

Warianty analizowane przez ten toolkit opierają się na opublikowanych badaniach i publicznie dostępnych bazach danych, w tym:
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [PharmGKB](https://www.pharmgkb.org/)
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
- [SNPedia](https://www.snpedia.com/)

Jednak nauka genetyczna jest złożona i stale się rozwija. Wiele asocjacji ma małe rozmiary efektów, ograniczoną replikację lub efekty specyficzne dla populacji. To narzędzie zapewnia uproszczony widok edukacyjny i nie powinno być mylone z interpretacją kliniczną.

---

## 🤝 Współpraca

Wkład jest mile widziany! Zapraszamy do tworzenia Pull Requestów.

### Pomysły na Wkład
- [ ] Dodatkowe warianty genetyczne z cytowaniami
- [ ] Więcej skryptów analitycznych
- [ ] Ulepszone wizualizacje
- [ ] Interfejs webowy
- [ ] Kontener Docker
- [ ] Wsparcie dla dodatkowych języków

---

## 📜 Licencja

Licencja MIT - zobacz plik [LICENSE](LICENSE).

**Dodatkowe Warunki:** Używając tego oprogramowania, zgadzasz się z powyższym zastrzeżeniem i potwierdzasz, że nie jest to przeznaczone do użytku klinicznego ani medycznego.

---

## ⚠️ Końcowe Przypomnienie

> **To jest hobbystyczny projekt do eksploracji i edukacji genetycznej.**
>
> 🎓 Używaj go do nauki o genetyce
> 🎮 Używaj go dla zabawy i ciekawości
> 🔬 Używaj go do osobistych badań
>
> 🚫 NIE używaj go do decyzji medycznych
> 🚫 NIE używaj go jako narzędzia diagnostycznego
> 🚫 NIE pomijaj profesjonalnego poradnictwa genetycznego jeśli masz obawy

---

<p align="center">
  Stworzone z 🧬 dla ciekawskich
</p>

<p align="center">
  <sub>Pamiętaj: Twoje geny nie są Twoim przeznaczeniem. Styl życia, środowisko i wybory też mają znaczenie!</sub>
</p>
