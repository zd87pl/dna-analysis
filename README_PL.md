# 🧬 Helixight - Otwarte Narzędzie do Analizy Genetycznej

**Wersja 1.3.0**

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
- **🌍 Dwujęzyczne menu CLI** - menu po polsku i angielsku (`L` w menu lub `HELIXIGHT_LANGUAGE=pl`); raporty z analiz są obecnie mieszanką polskiego i angielskiego, a interfejs webowy jest tylko po angielsku
- **📊 Kompleksowy** - 13 modułów analiz w interfejsie webowym; menu CLI oferuje 12 (dodaje narzędzia STR i EVO2 oparte na BAM, ale nie ma farmakogenomiki, badania nosicielstwa i grupy krwi, które można uruchomić bezpośrednio ze `scripts/`); sama Mega Analiza obejmuje ~80 wariantów w 10 sekcjach
- **🎓 Edukacyjny** - Poznaj swoją genetykę w zabawny, przystępny sposób

---

## 📋 Co Możesz Zbadać?

| Kategoria | Opis | Przykłady |
|-----------|------|-----------|
| 🏃 **Genetyka Sportowa** | Typ włókien mięśniowych, potencjał wytrzymałościowy | ACTN3, ACE, PPARGC1A |
| 🏊 **Predyspozycje Triatlonowe** | VO2max, metabolizm tłuszczów, ryzyko kontuzji | Tendencje Sprint vs Ironman |
| 💊 **Spersonalizowane Wnioski** | Styl treningu, suplementacja, żywienie | Na podstawie MTHFR, COMT, CYP1A2 |
| 💉 **Farmakogenomika** | Metabolizm leków, odpowiedź na leki | CYP2D6, CYP2C19, VKORC1 |
| 🧪 **Badanie nosicielstwa** | Nosicielstwo chorób recesywnych | Mukowiscydoza, anemia sierpowata, Tay-Sachs |
| 🩸 **Grupa krwi** | ABO i czynnik Rh z genotypu | Grupa krwi, zgodność transfuzji |
| 🧬 **Kompleksowa Analiza** | Mega Analiza: ~80 wariantów w 10 sekcjach | Pochodzenie, nosicielstwo, chronotyp, zmysły, długowieczność |
| 🔬 **Eksploracja Wariantów** | Znane istotne warianty | Eksploracja badawcza |
| 📊 **Wyniki Poligeniczne** | Połączone efekty wariantów | Edukacyjna eksploracja ryzyka |
| 🎲 **Zabawna Genetyka** | Cechy i wskazówki pochodzenia | Kolor oczu, kofeina, kolendra |

---

## 🚀 Szybki Start

### Opcja 1: Docker (Zalecana) 🐳

Najprostszy sposób uruchomienia Helixight to Docker - bez instalacji!

```bash
# Sklonuj repozytorium
git clone https://github.com/zd87pl/dna-analysis.git
cd dna-analysis

# Kontener działa jako użytkownik bez uprawnień roota; uruchom go jako TWÓJ
# użytkownik, aby mógł zapisywać do podmontowanych ./data i ./results (bash nie
# eksportuje UID, a nowy terminal lub `sudo docker compose` go gubi). Zapis
# identyfikatorów w .env (ignorowanym przez git) załatwia to na stałe:
echo "UID=$(id -u)" >> .env; echo "GID=$(id -g)" >> .env

# Uruchom interfejs webowy
docker compose up -d

# Otwórz http://localhost:8501 w przeglądarce
```

To wszystko! Interfejs webowy będzie dostępny pod adresem **http://localhost:8501**

Wypróbuj od razu na dołączonej próbce: `data/sample_genome.vcf.gz` pojawia się na liście
"Available VCF files" dzięki montowaniu katalogu `./data` (obraz zawiera też kopię w `/data`
dla zwykłego `docker run`).

#### Wskazówki Docker

```bash
# Podgląd logów
docker compose logs -f

# Zatrzymanie kontenera
docker compose down

# Przebudowa po aktualizacji
docker compose build --no-cache
docker compose up -d
```

- **Dostęp sieciowy**: aplikacja nie ma uwierzytelniania, dlatego `docker-compose.yml`
  wiąże port tylko z `127.0.0.1:8501`. Aby korzystać z niej z innego komputera, ustaw
  `HELIXIGHT_BIND=0.0.0.0` (lub adres konkretnego interfejsu) **i** postaw przed nią
  reverse proxy z uwierzytelnianiem i TLS; nie wystawiaj jej do sieci, której nie ufasz.
- **Właściciel plików**: `docker-compose.yml` uruchamia kontener jako `${UID}:${GID}`, a gdy
  te zmienne nie są ustawione, jako uid 10001. W Linuksie podmontowane `./data` i `./results`
  należą do Ciebie, więc kontener może do nich pisać tylko działając z Twoim uid: wpisz
  `UID=...`/`GID=...` do `.env` obok `docker-compose.yml` (patrz szybki start) albo wykonaj
  `export UID; export GID=$(id -g)` w każdej powłoce przed `docker compose up` (nowy
  terminal lub `sudo docker compose up` gubi eksport). Objaw braku ustawienia: upload kończy
  się `Upload failed: [Errno 13] Permission denied`, a każda analiza (także dołączona
  próbka) `Analysis failed: [Errno 13] Permission denied: /app/results/...`. Osobna
  przyczyna z tym samym objawem: jeśli `./data` lub `./results` utworzyła wcześniejsza
  (rootowa) wersja obrazu, popraw je raz poleceniem
  `sudo chown -R "$(id -u):$(id -g)" data results`.
- **Pamięć podręczna wyników**: każdy wynik analizy (pełny tekst raportu i ustalenia) jest
  zapisywany jako zwykły JSON w `results/.cache/`, dzięki czemu ponowne uruchomienie na tym
  samym pliku jest natychmiastowe. Przycisk **Clear Cache** w panelu bocznym (lub usunięcie
  `results/.cache/`) kasuje go; przełącznik cache w panelu bocznym wyłącza go dla nowych analiz.
- **Rozmiar uploadu**: upload przez przeglądarkę jest ograniczony do 1000 MB
  (`STREAMLIT_SERVER_MAX_UPLOAD_SIZE` w `docker-compose.yml`), ponieważ Streamlit trzyma cały
  plik w pamięci. Większe pliki skopiuj do `./data` i wybierz je z listy "Available VCF files".
- **Pamięć podręczna referencji**: `clinvar_scan.sh` pobiera ClinVar (~70 MB) raz do
  `results/.reference/`; zamontuj tam `./reference` (zakomentowany wolumen w
  `docker-compose.yml`), jeśli chcesz trzymać go poza `results/`. Pobrany plik jest
  weryfikowany sumą `.md5` NCBI; `HELIXIGHT_SKIP_CHECKSUM=1` pomija weryfikację, a
  `HELIXIGHT_CLINVAR_MAX_AGE_DAYS` (domyślnie 90) określa, po ilu dniach kopia jest
  odświeżana. Skrypt wymaga `flock` (util-linux) i `md5sum` (coreutils) oprócz
  `bcftools` i `wget`; uruchomiony z CLI używa katalogu `$HELIXIGHT_REFERENCE_DIR`
  (domyślnie `./clinvar_analysis/`).
- **Układ wyników**: interfejs webowy zapisuje każde uruchomienie w
  `results/<nazwa vcf>_<hash12>/<analysis_id>/`, menu CLI w
  `helixight_results/<nazwa-pliku>/` (raport zbiorczy w `helixight_results/`).
- **Zmienne środowiskowe** ustawiane przez interfejs webowy dla skryptów: `HELIXIGHT_WORK_DIR`
  (`results/.work`, robocze kopie bgzip/indeks), `HELIXIGHT_REFERENCE_DIR`
  (`results/.reference`, cache ClinVar) oraz `HELIXIGHT_ANALYSIS_TIMEOUT` (sekundy na
  analizę, domyślnie 3600, odczytywana ze środowiska kontenera).
- **Pamięć RAM**: plik compose ogranicza kontener do 8 GB (`deploy.resources` dla
  `docker compose`, `mem_limit` dla starego pliku wykonywalnego `docker-compose` v1).

#### Dodawanie własnych danych

Umieść pliki VCF w katalogu `data/` lub wgraj je przez interfejs webowy:

```bash
cp twoje_warianty.vcf.gz ./data/
```

**Obsługiwane formaty:**
- `.vcf` - tekstowy Variant Call Format
- `.vcf.gz` - VCF skompresowany gzip lub bgzip (indeks nie jest wymagany)
- nazwy kontigów z prefiksem `chr` (`chr1`, `chrM`) są akceptowane

Każdy skrypt przed odpytaniem konwertuje wejście do posortowanej, skompresowanej bgzip i
zindeksowanej tabix kopii roboczej z gołymi nazwami kontigów (w `results/.work/` w interfejsie
webowym, a przy bezpośrednim uruchomieniu skryptu w `$HELIXIGHT_WORK_DIR` lub
`$TMPDIR/helixight-<uid>` (domyślnie `/tmp/helixight-<uid>`, tworzony z prawami 700)). Kopia,
jej indeks i log sortowania to pełna kopia Twojego genomu i pozostają tam, dopóki nie usuniesz
tego katalogu. Nie musisz niczego przygotowywać. Jeśli chcesz zrobić to
raz samodzielnie, równoważne polecenie to:

```bash
bcftools sort -Oz -o twoje_warianty.vcf.gz wejscie.vcf && bcftools index -t twoje_warianty.vcf.gz
```

**VCF musi być wywołany względem GRCh38.** Pliki GRCh37/hg19 są odrzucane z czytelnym
komunikatem (kod wyjścia 3) zamiast być cicho źle interpretowane, ponieważ wszystkie
współrzędne w skryptach są GRCh38.

Wyniki analiz zapisywane są w katalogu `results/`.

---

### Opcja 2: Instalacja z linii poleceń

```bash
# Sklonuj repozytorium
git clone https://github.com/zd87pl/dna-analysis.git
cd dna-analysis

# Zainstaluj zależności (Ubuntu/Debian)
sudo apt install bcftools samtools tabix wget bc

# Lub użyj instalatora
./install.sh

# Uruchom interaktywne menu
./helixight.sh
```

---

## 📦 Wymagania

### Dla Dockera (Zalecane)
- Docker Engine 20.10+
- Docker Compose v2.0+ (`docker compose`; stary plik wykonywalny `docker-compose` v1 działa, ale ignoruje `deploy.resources`, dlatego ustawione jest też `mem_limit`)
- minimum 4 GB RAM (8 GB zalecane dla dużych plików VCF)

### Dla instalacji CLI

| Narzędzie | Cel | Instalacja |
|-----------|-----|------------|
| bcftools | Manipulacja VCF | `apt install bcftools` |
| samtools | Przetwarzanie BAM/SAM | `apt install samtools` |
| tabix | Indeksowanie VCF (bgzip/tabix) | `apt install tabix` |
| wget | Pobieranie plików | `apt install wget` |
| bc | Arytmetyka wyników w skryptach | `apt install bc` |
| flock, md5sum | Blokada i suma kontrolna cache ClinVar (`clinvar_scan.sh`) | util-linux / coreutils (preinstalowane w Linuksie) |

### Instalacja jedną komendą

**Ubuntu/Debian:**
```bash
sudo apt update && sudo apt install -y bcftools samtools tabix wget gzip curl bc
```

**macOS (Homebrew):**
```bash
brew install bcftools samtools htslib wget
```

`install.sh` korzysta z tych menedżerów pakietów; na macOS wymaga wcześniej zainstalowanego [Homebrew](https://brew.sh) (w przeciwnym razie kończy się komunikatem). macOS ma wbudowane `bc`.

Uwaga dla macOS: `clinvar_scan.sh` (opcje CLI `20` i `41`) sprawdza pobrany plik ClinVar
przez `md5sum` i serializuje równoległe pobrania przez `flock`; macOS nie ma żadnego z nich.
`brew install coreutils` dostarcza `md5sum` (jako `gmd5sum` lub `md5sum` z gnubin w PATH);
bez niego skrypt używa wbudowanego `md5`. Bez `flock` skan działa, ale dwa jednoczesne
pobrania ClinVar nie są serializowane (wypisywane jest ostrzeżenie). Obraz Docker ma oba.

---

## 📖 Użycie

### Tryb Interaktywny (Zalecany)

```bash
./helixight.sh
```

Uruchamia interaktywne menu z opcjami:
- `1` Instalacja zależności, `2` Sprawdzenie zależności, `3` Pobranie genomu referencyjnego (jeśli masz pliki BAM)
- `4` Utworzenie VCF z plików BAM, `5` Wybór istniejącego pliku VCF
- `10`-`14` Analizy genetyczne (sportowa, triathlon, spersonalizowany protokół, mega analiza, zabawna genetyka)
- `20`-`23` Zdrowie i klinika (skan ClinVar, geny nowotworowe, poligeniczne skale ryzyka, rzadkie warianty)
- `30`-`32` Zaawansowane (analiza STR/powtórzeń i przygotowanie danych EVO2 wymagają BAM; haplogrupa mtDNA)
- `40` Generowanie raportu zbiorczego (tekstowy; konwertowany do PDF tylko gdy zainstalowany jest pandoc z silnikiem PDF, np. LaTeX - gwarantowany eksport PDF jest w interfejsie webowym), `41` Uruchomienie 9 analiz VCF z opcji 10-14 i 20-23, a następnie raportu zbiorczego
- `L` Zmiana języka (EN/PL), `H` Pomoc, `Q` Wyjście

### Bezpośrednie Uruchamianie Skryptów

```bash
# Eksploracja genetyki sportowej
./scripts/athletic_genetics.sh twoje_warianty.vcf.gz

# Analiza predyspozycji triatlonowych
./scripts/triathlon_genetics.sh twoje_warianty.vcf.gz

# Spersonalizowane wnioski
./scripts/personalized_protocol.sh twoje_warianty.vcf.gz

# Kompleksowa eksploracja (~80 wariantów w 10 sekcjach)
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
Kompleksowa eksploracja ~80 wariantów w 10 sekcjach:
- Markery rodowe (DNA neandertalskie)
- Eksploracja statusu nosicielstwa
- Chronotyp (sowa vs skowronek)
- Genetyka zmysłów (smak, węch)
- Przewidywanie cech

---

## 📁 Struktura Projektu

```
dna-analysis/
├── helixight.sh          # Główny interaktywny launcher CLI
├── install.sh            # Szybki instalator dla CLI
├── Dockerfile            # Definicja kontenera Docker
├── docker-compose.yml    # Konfiguracja Docker Compose
├── pytest.ini            # Konfiguracja pytest
├── README.md             # Dokumentacja angielska
├── README_PL.md          # Dokumentacja polska
├── LICENSE               # Licencja MIT
├── .github/
│   ├── workflows/ci.yml  # CI: pytest, shellcheck, uruchomienie skryptów, budowa obrazu Docker
│   └── dependabot.yml    # Kontrola aktualizacji zależności (pip, docker, actions)
├── scripts/              # Skrypty analityczne (bash + bcftools)
│   ├── lib/
│   │   ├── vcf_common.sh # Wspólne funkcje (konwersja wejścia, odczyt genotypów)
│   │   └── GRCh38_COORDINATES.md  # Kanoniczna tabela współrzędnych GRCh38
│   ├── athletic_genetics.sh
│   ├── triathlon_genetics.sh
│   ├── personalized_protocol.sh
│   ├── mega_analysis.sh
│   ├── clinvar_scan.sh
│   ├── cancer_genes_scan.sh
│   └── ... (więcej skryptów)
├── frontend/             # Interfejs webowy (Streamlit)
│   ├── app.py            # Główna aplikacja Streamlit
│   ├── analysis.py       # Wrapper Pythona dla skryptów
│   ├── pdf_report.py     # Generator raportów PDF
│   └── requirements.txt  # Zależności Pythona
├── tests/                # Testy
│   ├── conftest.py       # Środowisko testów (RESULTS_DIR, sys.path)
│   ├── test_analysis.py  # Testy modułu analysis
│   ├── test_pdf_report.py# Testy generowania PDF
│   ├── shell/test_vcf_common.sh  # Testy wspólnej biblioteki shell
│   ├── shell/test_scripts.sh     # Testy regresyjne skryptów (syntetyczne VCF)
│   └── requirements.txt  # Zależności testów
├── data/                 # Tu umieść swoje pliki VCF
│   ├── sample_genome.vcf        # Przykładowe dane (GRCh38, syntetyczne genotypy)
│   ├── sample_genome.vcf.gz     # To samo, skompresowane bgzip
│   └── sample_genome.vcf.gz.tbi # ... i indeks tabix
├── results/              # Wyniki analiz (tu żyją .cache/, .work/, .reference/)
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

## 🛠️ Rozwój / Uruchamianie testów

W CI i w obrazie Docker używany jest Python 3.11.

```bash
# Testy Pythona (wrappery frontendu i raport PDF)
python -m venv .venv && . .venv/bin/activate
pip install -r tests/requirements.txt
python -m pytest

# Testy wspólnej biblioteki shell (wymagają bcftools, bgzip, tabix)
bash tests/shell/test_vcf_common.sh

# Testy regresyjne skryptów (syntetyczne VCF; wymagają bcftools)
bash tests/shell/test_scripts.sh

# Lint skryptów shell
shellcheck -S error scripts/*.sh scripts/lib/*.sh helixight.sh install.sh
```

Skrypty można uruchamiać bez interfejsu webowego; rób to z katalogu roboczego poza
repozytorium (katalogi wynikowe tworzą w bieżącym katalogu):

```bash
mkdir -p /tmp/hx && cd /tmp/hx
HELIXIGHT_WORK_DIR=/tmp/hx/work bash /sciezka/do/dna-analysis/scripts/athletic_genetics.sh /sciezka/do/dna-analysis/data/sample_genome.vcf.gz
```

Każdy skrypt kończy raport linią `HX_DATA_COVERAGE: <znalezione>/<odpytane> markers`;
interfejs webowy odróżnia dzięki niej "brak użytecznych danych w pliku" od prawdziwego wyniku.

CI (`.github/workflows/ci.yml`) uruchamia się przy każdym pushu do `main` i każdym pull
requeście: testy pytest z pokryciem, `shellcheck`, uruchomienie każdego skryptu analiz z
interfejsu webowego na przykładowym VCF (gołe i `chr`-prefiksowane kontigi) oraz budowę
obrazu Docker.

---

## 🤝 Współpraca

Wkład jest mile widziany! Zapraszamy do tworzenia Pull Requestów. Prosimy aktualizować
`README.md` i `README_PL.md` razem.

### Pomysły na Wkład
- [ ] Dodatkowe warianty genetyczne z cytowaniami
- [ ] Więcej skryptów analitycznych
- [ ] Ulepszone wizualizacje
- [x] Interfejs webowy ✅
- [x] Kontener Docker ✅
- [ ] Wsparcie dla dodatkowych języków
- [x] Generowanie raportów PDF ✅
- [x] Metryki jakości VCF ✅
- [x] Pamięć podręczna wyników ✅
- [x] Przykładowe dane testowe ✅
- [x] Interaktywne wykresy wyników ✅

---

## 📝 Historia zmian

### v1.3.0
- **Poprawka**: Normalizacja wejścia - każdy skrypt konwertuje zwykły, gzip, bgzip i `chr`-prefiksowany VCF do posortowanej, zindeksowanej kopii roboczej z gołymi nazwami kontigów (`scripts/lib/vcf_common.sh`) zamiast kończyć się błędem lub cicho zwracać "brak danych"; pliki GRCh37 są odrzucane z czytelnym komunikatem
- **Poprawka**: Poprawione współrzędne GRCh38 i interpretacje (COMT, CYP1A2, APOE, LCT, ALDH2, FTO, TCF7L2, MTHFR, VDR, HFE, CYP2D6, SLCO1B1, TPMT, GBA i inne; kanoniczna tabela w `scripts/lib/GRCh38_COORDINATES.md`); odczyt genotypu sprawdza REF/ALT przed punktacją
- **Poprawka**: Skrypty raportują `HX_DATA_COVERAGE`, dzięki czemu interfejs webowy odróżnia "brak użytecznych danych" od czystego wyniku; naprawione błędy liczników przerywające `pharmacogenomics.sh`, `carrier_screening.sh` i `clinvar_scan.sh` pod `set -e`
- **Poprawka**: Eksport PDF nie zawiesza się już na emoji/Unicode w ustaleniach ani na długich liniach
- **Poprawka**: Walidacja uploadu w interfejsie webowym (oczyszczanie nazwy pliku, sprawdzenie nagłówka VCF, brak podwójnej kopii w pamięci); osobne katalogi robocze per plik, żeby równoległe analizy nie nadpisywały swoich wyników
- **Poprawka**: Przykładowy VCF wygenerowany na nowo - czytelny dla bcftools, unikalne ID/pozycje, współrzędne GRCh38, zamodelowane ACE I/D oraz delecje ABO/CFTR, łagodny rekord BRCA2 i sześć rekordów mtDNA (rCRS), dzięki którym moduły nowotworowy i mtDNA kończą się wynikiem na próbce; dostarczany jako bgzip z indeksem tabix
- **Bezpieczeństwo**: Kontener działa jako użytkownik bez uprawnień roota; port webowy domyślnie związany z `127.0.0.1`; limit uploadu obniżony do 1000 MB; pobranie ClinVar weryfikowane sumą `.md5`
- **Ulepszenie**: Workflow CI (pytest z pokryciem, shellcheck, uruchomienie skryptów na przykładowym VCF, budowa obrazu Docker) oraz Dependabot
- **Ulepszenie**: Uczciwe opisy funkcji - liczby wariantów, numeracja menu i wsparcie językowe w dokumentacji zgodne z kodem
- **Poprawka**: Grupa krwi - rs8176746 (allel B) czytany w orientacji nici plus (G>T), więc genomy B/B i A/B nie są już raportowane jako "A or O"; genomy A/O i A/A w VCF zawierających tylko warianty (brak rekordu rs7853989) dostają wynik zamiast "Unknown"
- **Poprawka**: Niezweryfikowane sondy z zakresu GRCh37 nie generują już wyników: pozycje nazwanych mutacji BRCA1/2 są notatkami informacyjnymi (nie liczą się jako patogenne), Tay-Sachs w `mega_analysis.sh` używa tej samej sondy insercji TATC co `carrier_screening.sh`, COL5A1 rs12722 i NFE2L2 są pomijane, PPARA rs4253778 punktowane w tej samej orientacji (wytrzymałość = G) w obu skryptach fitness, wiersze HBB HbC/IVS-I-1 przesunięte obok rs334; CYP2D6*2 (normalna funkcja) jest informacją, nie czerwonym ustaleniem
- **Poprawka**: Wspólna biblioteka - VCF z niezadeklarowanymi tagami INFO/FORMAT są akceptowane, mieszane rekordy SNP+indel zachowują wywołanie SNP, kopie robocze trafiają do katalogu per użytkownik z prawami 700, brak rozszerzeń tylko dla bash 4 (macOS bash 3.2), a linia `##contig` bez `length=` wywołuje ostrzeżenie "nie można zweryfikować buildu"
- **Poprawka**: Interfejs webowy - pustego VCF (0 rekordów) nie da się analizować, wpisy cache z wersji 1.2.0 są ignorowane (wersja formatu), skrypt przerwany limitem czasu/anulowaniem, którego potomek trzyma potok, jest zabijany po 5 s, `output_files` zawiera tylko pliki bieżącego uruchomienia

### v1.2.0
- **Funkcja**: Analiza farmakogenomiczna (PGx) - genetyka metabolizmu leków (CYP2D6, CYP2C19, CYP2C9, VKORC1 itd.)
- **Funkcja**: Badanie nosicielstwa - status nosiciela chorób recesywnych (mukowiscydoza, anemia sierpowata, Tay-Sachs itd.)
- **Funkcja**: Przewidywanie grupy krwi - ABO i czynnik Rh z wariantów genetycznych
- **Ulepszenie**: Rozszerzony przykładowy VCF do testowania wszystkich analiz
- **Ulepszenie**: 14 modułów analiz

### v1.1.0
- **Funkcja**: Generowanie raportów PDF z profesjonalnym formatowaniem i wykresami
- **Funkcja**: Metryki jakości VCF (liczba wariantów, stosunek Ti/Tv, rozkład genotypów)
- **Funkcja**: Pamięć podręczna wyników, aby nie uruchamiać ponownie analiz na tym samym pliku
- **Funkcja**: Przykładowe dane VCF (`data/sample_genome.vcf`)
- **Funkcja**: Zestaw testów jednostkowych pytest dla modułów analysis i PDF
- **Ulepszenie**: Przełącznik cache w ustawieniach UI
- **Ulepszenie**: Przycisk czyszczenia cache

### v1.0.1
- **Bezpieczeństwo**: Naprawione podatności na wstrzykiwanie poleceń w skryptach shell
- **Bezpieczeństwo**: Walidacja wszystkich ścieżek plików i danych wejściowych użytkownika
- **Bezpieczeństwo**: Pliki blokady zapobiegające wyścigom przy pobieraniu ClinVar
- **Poprawka**: Walidacja uploadu plików `.vcf.gz` w interfejsie webowym
- **Poprawka**: Zagnieżdżenie przycisku pobierania w Streamlit
- **Poprawka**: Potencjalny KeyError przy wyświetlaniu wyników
- **Poprawka**: Nieosiągalny kod przy anulowaniu analizy
- **Poprawka**: Brakująca zmienna środowiskowa `RESULTS_DIR` w Dockerze
- **Ulepszenie**: Obsługa błędów z `set -euo pipefail` we wszystkich skryptach
- **Ulepszenie**: Zastąpienie wartości na sztywno konfigurowalnymi stałymi
- **Ulepszenie**: Zgodność z Pythonem 3.8+ dla adnotacji typów

### v1.0.0
- Pierwsze wydanie z CLI i interfejsem webowym w Dockerze
- 22 skrypty analityczne obejmujące genetykę sportową, zdrowotną i zabawną
- Frontend webowy oparty na Streamlit
- Konteneryzacja Docker z narzędziami bioinformatycznymi

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
