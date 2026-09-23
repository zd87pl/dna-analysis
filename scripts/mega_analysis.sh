#!/bin/bash
#===============================================================================
#
#    ███╗   ███╗███████╗ ██████╗  █████╗      █████╗ ███╗   ██╗ █████╗ ██╗     ██╗███████╗ █████╗ 
#    ████╗ ████║██╔════╝██╔════╝ ██╔══██╗    ██╔══██╗████╗  ██║██╔══██╗██║     ██║╚══███╔╝██╔══██╗
#    ██╔████╔██║█████╗  ██║  ███╗███████║    ███████║██╔██╗ ██║███████║██║     ██║  ███╔╝ ███████║
#    ██║╚██╔╝██║██╔══╝  ██║   ██║██╔══██║    ██╔══██║██║╚██╗██║██╔══██║██║     ██║ ███╔╝  ██╔══██║
#    ██║ ╚═╝ ██║███████╗╚██████╔╝██║  ██║    ██║  ██║██║ ╚████║██║  ██║███████╗██║███████╗██║  ██║
#    ╚═╝     ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝    ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝╚══════╝╚═╝  ╚═╝
#
#    KOMPLEKSOWA ANALIZA GENETYCZNA - WSZYSTKO W JEDNYM!
#    Helixight Genetic Analysis Platform
#
#===============================================================================

# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"
set -eo pipefail

VCF="${1:-}"
hx_require_arg "$VCF" "$(basename "$0") <variants.vcf.gz>"
OUTDIR="mega_analysis_results"
mkdir -p "$OUTDIR"


# Funkcja do sprawdzania SNP (wspolrzedne GRCh38, scripts/lib/GRCh38_COORDINATES.md;
# genotyp pierwszego SNP na pozycji, faza znormalizowana, brak wywolania = pusty)
check_snp() {
    hx_gt "$1"
}

# "REF ALT GT" pierwszego SNP na pozycji (pusty gdy brak rekordu / brak wywolania)
check_snp_full() {
    local out
    out=$(hx_query "$1" '%REF %ALT [%GT]\n' -e 'GT="mis" || TYPE!="snp"')
    printf '%s\n' "${out%%$'\n'*}"
}

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                       ║"
echo "║     ██╗  ██╗███████╗██╗     ██╗██╗  ██╗██╗ ██████╗ ██╗  ██╗████████╗                  ║"
echo "║     ██║  ██║██╔════╝██║     ██║╚██╗██╔╝██║██╔════╝ ██║  ██║╚══██╔══╝                  ║"
echo "║     ███████║█████╗  ██║     ██║ ╚███╔╝ ██║██║  ███╗███████║   ██║                     ║"
echo "║     ██╔══██║██╔══╝  ██║     ██║ ██╔██╗ ██║██║   ██║██╔══██║   ██║                     ║"
echo "║     ██║  ██║███████╗███████╗██║██╔╝ ██╗██║╚██████╔╝██║  ██║   ██║                     ║"
echo "║     ╚═╝  ╚═╝╚══════╝╚══════╝╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝ ╚═╝  ╚═╝   ╚═╝                     ║"
echo "║                                                                                       ║"
echo "║                    M E G A   A N A L I Z A   G E N E T Y C Z N A                      ║"
echo "║                                                                                       ║"
echo "╚═══════════════════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Plik VCF: $VCF"
echo "Wyniki: $OUTDIR/"
echo ""
echo "Rozpoczynam kompleksowa analize... To moze potrwac kilka minut."
echo ""

# Biblioteka normalizuje nazwy kontigow i buduje zindeksowana kopie robocza
hx_prepare_vcf "$VCF"

#===============================================================================
# SEKCJA 1: DNA NEANDERTALCZYKA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🦴 SEKCJA 1: DNA NEANDERTALCZYKA I DENISOWIAŃCZYKA                       ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

NEANDER_SCORE=0

# Znane markery neandertalskie
echo "Sprawdzam markery archaiczne..."
echo ""

# BNC2 - jasna skora (od neandertalczykow)
GT=$(check_snp "9:16817033")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ rs2066807 (BNC2): Wariant neandertalski - jasniejsza skora"
    NEANDER_SCORE=$((NEANDER_SCORE + 1))
fi

# SLC45A2 - pigmentacja
GT=$(check_snp "5:33951693")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ rs16891982 (SLC45A2): Wariant europejski pigmentacji"
    NEANDER_SCORE=$((NEANDER_SCORE + 1))
fi

# OAS1 - odpornosc (od neandertalczykow)
GT=$(check_snp "12:113357193")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ rs10774671 (OAS1): Wariant neandertalski - lepsza odpornosc antywirusowa!"
    NEANDER_SCORE=$((NEANDER_SCORE + 1))
fi

# TLR1/6/10 - odpornosc wrodzona
GT=$(check_snp "4:38775150")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ rs5743618 (TLR1): Wariant neandertalski - odpornosc na bakterie"
    NEANDER_SCORE=$((NEANDER_SCORE + 1))
fi

# STAT2 - odpowiedz na wirusy
GT=$(check_snp "12:56736371")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ rs2066807 (STAT2): Wariant archaiczny - odpowiedz immunologiczna"
    NEANDER_SCORE=$((NEANDER_SCORE + 1))
fi

# SLC16A11 - metabolizm (od neandertalczykow, zwiazany z cukrzyca)
GT=$(check_snp "17:7037074")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ⚠ rs13342232 (SLC16A11): Wariant neandertalski - ryzyko cukrzycy"
    NEANDER_SCORE=$((NEANDER_SCORE + 1))
fi

# EPAS1 - adaptacja wysokogórska (od Denisowian!)
GT=$(check_snp "2:46441523")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ rs142764723 (EPAS1): Wariant DENISOWIANSKI - adaptacja do wysokosci!"
    NEANDER_SCORE=$((NEANDER_SCORE + 1))
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "WYNIK: Znaleziono $NEANDER_SCORE markerow archaicznych"
echo ""
echo "Typowy Europejczyk ma 1-4% DNA neandertalskiego (~50-100 regionow)."
echo "Ta analiza sprawdza tylko wybrane, dobrze zbadane markery."
echo ""

if [ $NEANDER_SCORE -ge 4 ]; then
    echo "→ Masz WIELE wariantow neandertalskich - typowe dla Europejczykow!"
elif [ $NEANDER_SCORE -ge 2 ]; then
    echo "→ Masz kilka wariantow neandertalskich - w normie europejskiej."
else
    echo "→ Malo wykrytych markerow (moze byc kwestia pokrycia VCF)."
fi
echo ""

#===============================================================================
# SEKCJA 2: NOSICIELSTWO CHOROB RECESYWNYCH
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🧬 SEKCJA 2: NOSICIELSTWO CHOROB RECESYWNYCH                             ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Wazne jesli planujesz potomstwo! Nosiciel jest zdrowy, ale moze"
echo "przekazac chorobe dziecku jesli partner tez jest nosicielem."
echo ""

CARRIER_COUNT=0

# Mukowiscydoza (CF) - CFTR
echo "Sprawdzam mukowiscydoze (CFTR)..."
# F508del - najczestsza mutacja (delecja 3 bp ATCT>A, rs113993960: nie SNP,
# wiec nie przez check_snp/hx_gt)
RES=$(hx_gt_allele "7:117559590" ATCT A rs113993960); GT=${RES#*$'\t'}
if [[ "$GT" == *"1"* ]]; then
    echo "  ⚠️ CFTR F508del: NOSICIEL mukowiscydozy!"
    CARRIER_COUNT=$((CARRIER_COUNT + 1))
else
    echo "  ✅ CFTR F508del: Nie wykryto"
fi

# Anemia sierpowata - HBB
echo "Sprawdzam anemie sierpowata (HBB)..."
GT=$(check_snp "11:5227002")
if [[ "$GT" == *"1"* ]]; then
    echo "  ⚠️ HBB rs334: NOSICIEL anemii sierpowatej!"
    CARRIER_COUNT=$((CARRIER_COUNT + 1))
else
    echo "  ✅ HBB rs334: Nie wykryto"
fi

# Talasemia beta
GT=$(check_snp "11:5226774")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ⚠️ HBB: Wariant talasemii beta wykryty"
    CARRIER_COUNT=$((CARRIER_COUNT + 1))
fi

# Fenyloketonuria (PKU) - PAH
echo "Sprawdzam fenyloketonurie (PAH)..."
# R408W, GRCh38 12:102840632 (GRCh37 was 12:103234258)  # TODO(build): unverified GRCh38 (table row is 'unverified')
RES=$(hx_gt_allele "12:102840632" G A rs5030858); GT=${RES#*$'\t'}
if [[ "$GT" == *"1"* ]]; then
    echo "  ⚠️ PAH rs5030858: NOSICIEL fenyloketonurii!"
    CARRIER_COUNT=$((CARRIER_COUNT + 1))
else
    echo "  ✅ PAH: Nie wykryto"
fi

# Choroba Tay-Sachsa - HEXA
echo "Sprawdzam chorobe Tay-Sachsa (HEXA)..."
# c.1274_1277dupTATC (rs76173977): 4-bp insertion, not a SNP. The former probe
# 15:72638892 was the GRCh37 locus (~292 kb from GRCh38 HEXA); this is the
# position carrier_screening.sh uses (insertion anchored at 72346579/80).
# TODO(build): unverified GRCh38 position
REC=$(hx_gt_indel "15:72346579-72346580" 4)
GT=""
if [[ -n "$REC" ]]; then
    read -r _ HEXA_REF HEXA_ALT HEXA_GT <<<"$REC"
    # only an insertion of exactly 4 bp (relative to REF) counts
    if [[ ${#HEXA_ALT} -eq $(( ${#HEXA_REF} + 4 )) ]]; then
        GT=$HEXA_GT
    fi
fi
if [[ "$GT" == *"1"* ]]; then
    echo "  ⚠️ HEXA: Wariant Tay-Sachsa wykryty!"
    CARRIER_COUNT=$((CARRIER_COUNT + 1))
else
    echo "  ✅ HEXA: Nie wykryto"
fi

# Hemochromatoza - HFE
echo "Sprawdzam hemochromatoze (HFE)..."
GT=$(check_snp "6:26092913")  # C282Y
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
    echo "  🔴 HFE C282Y: HOMOZYGOTA - ryzyko hemochromatozy!"
    CARRIER_COUNT=$((CARRIER_COUNT + 2))
elif [[ "$GT" == *"1"* ]]; then
    echo "  🟡 HFE C282Y: NOSICIEL hemochromatozy"
    CARRIER_COUNT=$((CARRIER_COUNT + 1))
else
    echo "  ✅ HFE C282Y: Nie wykryto"
fi

GT=$(check_snp "6:26090951")  # H63D
if [[ "$GT" == *"1"* ]]; then
    echo "  🟡 HFE H63D: Wariant - lekko podwyzszone zelazo"
fi

# SMA - SMN1
echo "Sprawdzam SMA (SMN1)..."
# SMA wymaga analizy CNV, wiec tylko przyblizenie
echo "  ℹ️ SMA wymaga analizy liczby kopii (CNV) - niedostepne w podstawowym VCF"

# Niedobor alfa-1 antytrypsyny - SERPINA1
echo "Sprawdzam niedobor alfa-1 antytrypsyny..."
GT=$(check_snp "14:94378610")  # Z allele  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ⚠️ SERPINA1 Z: Wariant niedoboru A1AT"
    CARRIER_COUNT=$((CARRIER_COUNT + 1))
fi
GT=$(check_snp "14:94376747")  # S allele  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  🟡 SERPINA1 S: Lagodny wariant A1AT"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if [ $CARRIER_COUNT -eq 0 ]; then
    echo "✅ Nie wykryto nosicielstwa badanych chorob recesywnych"
else
    echo "⚠️ Wykryto $CARRIER_COUNT wariantow nosicielstwa"
    echo "→ Jesli planujesz potomstwo, rozważ badanie partnera"
fi
echo ""

#===============================================================================
# SEKCJA 3: ZDROWIE PSYCHICZNE
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🧠 SEKCJA 3: GENETYKA ZDROWIA PSYCHICZNEGO                               ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "UWAGA: Zdrowie psychiczne to zlozony fenotyp - geny to tylko czesc obrazu!"
echo "Srodowisko, doswiadczenia i styl zycia maja ogromne znaczenie."
echo ""

# ADHD
echo "📊 ADHD / Koncentracja:"
echo "─────────────────────────────────────────"

# DRD4 - receptor dopaminy D4
GT=$(check_snp "11:637269")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  DRD4 7R: Wariant zwiazany z ADHD i poszukiwaniem nowosci"
fi

# DAT1/SLC6A3 - transporter dopaminy
GT=$(check_snp "5:1392791")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  SLC6A3: Wariant transportera dopaminy"
fi

# SNAP25
GT=$(check_snp "20:10268206")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  SNAP25 rs363050: Wariant zwiazany z funkcjami poznawczymi"
fi

echo ""

# DEPRESJA
echo "📊 Depresja / Nastroj:"
echo "─────────────────────────────────────────"

# 5-HTTLPR / SLC6A4
GT=$(check_snp "17:30194319")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  SLC6A4 (transporter serotoniny): Wariant S - wrazliwosc na stres"
    echo "    → Wieksze znaczenie wsparcia spolecznego i terapii"
fi

# BDNF
GT=$(check_snp "11:27658369")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  BDNF Val66Val: Optymalna plastycznosc mozgu"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  BDNF Val66Met: Heterozygota - cwiczenia szczegolnie wazne!"
            ;;
        "1/1"|"1|1")
            echo "  BDNF Met66Met: Nizsza plastycznosc - cwiczenia KRYTYCZNE"
            echo "    → Regularna aktywnosc fizyczna podnosi BDNF!"
            ;;
    esac
fi

# COMT
GT=$(check_snp "22:19963748")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  COMT Val/Val: 'Warrior' - odporny na stres, mniej ruminacji"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  COMT Val/Met: Zrownowazony profil stresowy"
            ;;
        "1/1"|"1|1")
            echo "  COMT Met/Met: 'Worrier' - wieksza wrazliwosc na stres"
            echo "    → Mindfulness i techniki relaksacyjne szczegolnie pomocne"
            ;;
    esac
fi

# FKBP5 - odpowiedz na stres / PTSD
GT=$(check_snp "6:35656217")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  FKBP5 rs1360780: Wariant wrazliwosci na traume"
    echo "    → Wczesna interwencja po traumie szczegolnie wazna"
fi

echo ""

# UZALEZNIENIA
echo "📊 Predyspozycje do uzaleznien:"
echo "─────────────────────────────────────────"

# DRD2 Taq1A
GT=$(check_snp "11:113400106")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  DRD2 Taq1A (rs1800497): Mniej receptorow dopaminy"
    echo "    → Wieksza potrzeba stymulacji, ostroznosc z uzywkami"
fi

# OPRM1 - receptor opioidowy
GT=$(check_snp "6:154039662")
if [[ "$GT" == *"1"* ]]; then
    echo "  OPRM1 rs1799971: Zmieniona odpowiedz na opioidy"
fi

# ADH1B - metabolizm alkoholu. rs1229984 T>C: REF T = His48 = ADH1B*2 (szybki),
# ALT C = Arg48 = ADH1B*1 (czesty u Europejczykow); DOS = liczba alleli C.
# TODO(build): unverified GRCh38 orientation (table row is 'unverified')
RES=$(hx_gt_allele 4:99318162 T C rs1229984)
DOS=${RES%%$'\t'*}
case "$DOS" in
    0)
        echo "  ADH1B rs1229984: T/T (*2/*2, His48) → Szybki metabolizm alkoholu"
        echo "    → Mniejsze ryzyko alkoholizmu (nieprzyjemne objawy)"
        ;;
    1)
        echo "  ADH1B rs1229984: T/C (*1/*2) → Przyspieszony metabolizm alkoholu (jeden allel His48)"
        ;;
    2)
        echo "  ADH1B rs1229984: C/C (*1/*1, Arg48) → Standardowy (wolniejszy) metabolizm alkoholu, typowy dla Europejczykow"
        echo "    → Potencjalnie wieksze ryzyko uzaleznienia"
        ;;
esac

# CHRNA5 - nikotyna
GT=$(check_snp "15:78882925")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  CHRNA5 rs16969968: Wieksza wrazliwosc na nikotyne"
    echo "    → Wieksze ryzyko uzaleznienia od papierosow"
fi

# GABRA2 - alkohol
GT=$(check_snp "4:46252826")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  GABRA2 rs279858: Wariant zwiazany z alkoholizmem"
fi

echo ""

# SCHIZOFRENIA / CHOROBA DWUBIEGUNOWA (tylko informacyjnie)
echo "📊 Inne markery (informacyjne):"
echo "─────────────────────────────────────────"

# CACNA1C - choroba dwubiegunowa
GT=$(check_snp "12:2345295")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  CACNA1C rs1006737: Wariant zwiazany z zaburzeniami nastroju"
fi

# ZNF804A - schizofrenia
GT=$(check_snp "2:185778428")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ZNF804A rs1344706: Wariant ryzyka psychozy (maly efekt)"
fi

echo ""

#===============================================================================
# SEKCJA 4: CHRONOTYP - SOWA CZY SKOWRONEK
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🌙 SEKCJA 4: CHRONOTYP - ZEGAR BIOLOGICZNY                               ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

MORNING_SCORE=0
EVENING_SCORE=0

# PER2
GT=$(check_snp "2:239186871")  # TODO(build): unverified GRCh38
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
    echo "  PER2 rs12914385: T/T → Silna tendencja do wczesnego wstawania"
    MORNING_SCORE=$((MORNING_SCORE + 2))
elif [[ "$GT" == *"1"* ]]; then
    MORNING_SCORE=$((MORNING_SCORE + 1))
fi

# CLOCK
GT=$(check_snp "4:56294068")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  CLOCK rs1801260: Wariant → Tendencja do poznego chodzenia spac"
    EVENING_SCORE=$((EVENING_SCORE + 2))
fi

# PER3
GT=$(check_snp "1:7870203")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  PER3 rs2304672: Wariant → Krotszy sen, wieczorna aktywnosc"
    EVENING_SCORE=$((EVENING_SCORE + 1))
fi

# ARNTL/BMAL1
GT=$(check_snp "11:13387651")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ARNTL rs11121022: Wariant regulacji rytmu dobowego"
fi

# ADA - sen glęboki
GT=$(check_snp "20:44619933")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ADA rs73598374: Wariant → Wiecej snu glębokiego!"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
if [ $MORNING_SCORE -gt $EVENING_SCORE ]; then
    echo "WYNIK: 🐦 SKOWRONEK (ranny ptaszek)"
    echo "→ Najlepsza produktywnosc: 8:00-12:00"
    echo "→ Unikaj waznych zadan po 20:00"
elif [ $EVENING_SCORE -gt $MORNING_SCORE ]; then
    echo "WYNIK: 🦉 SOWA (nocny marek)"
    echo "→ Najlepsza produktywnosc: 16:00-22:00"
    echo "→ Nie planuj waznych spotkan przed 10:00"
else
    echo "WYNIK: 🕐 TYP POSREDNI"
    echo "→ Elastyczny rytm - dostosuj do stylu zycia"
fi
echo ""

#===============================================================================
# SEKCJA 5: NUTRIGENOMIKA - SPERSONALIZOWANA DIETA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🥗 SEKCJA 5: NUTRIGENOMIKA - DIETA DOPASOWANA DO GENOW                   ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "📊 Metabolizm makroskladnikow:"
echo "─────────────────────────────────────────"

# FTO - otylosc
GT=$(check_snp "16:53786615")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  FTO rs9939609: T/T → ✅ Niskie ryzyko otylosci genetycznej"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  FTO rs9939609: T/A → 🟡 Srednie ryzyko (+1.5 kg sredniej masy)"
            ;;
        "1/1"|"1|1")
            echo "  FTO rs9939609: A/A → 🔴 Podwyzszone ryzyko (+3 kg)"
            echo "    → CWICZENIA niweluja ten efekt calkowicie!"
            ;;
    esac
fi

# TCF7L2 - cukrzyca / weglowodany
GT=$(check_snp "10:112998590")
if [[ "$GT" == *"1"* ]]; then
    echo "  TCF7L2 rs7903146: Wariant → Gorsza tolerancja weglowodanow"
    echo "    → Dieta niskoglikemiczna, mniej cukrow prostych"
fi

# APOA2 - tluszcze nasycone
GT=$(check_snp "1:161222292")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  APOA2 rs5082: Wariant → Wrazliwosc na tluszcze nasycone"
    echo "    → Ogranicz: maslo, tlusty ser, czerwone mieso (<22g/dzien)"
fi

# PPARG - tluszcze
GT=$(check_snp "3:12351626")
if [[ "$GT" == *"1"* ]]; then
    echo "  PPARG rs1801282: Wariant Pro12Ala → Lepsza insulinowrazliwosc"
fi

echo ""
echo "📊 Witaminy i mineraly:"
echo "─────────────────────────────────────────"

# MTHFR - folian
GT=$(check_snp "1:11796321")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  MTHFR C677T: C/C → Normalny metabolizm folianu"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  MTHFR C677T: C/T → 🟡 Obnizony metabolizm folianu (65%)"
            echo "    → Kwas foliowy pozostaje skuteczny; L-metylofolian jest opcja, nie koniecznoscia"
            ;;
        "1/1"|"1|1")
            echo "  MTHFR C677T: T/T → 🔴 Znacznie obnizony metabolizm (30%)"
            echo "    → Kwas foliowy pozostaje skuteczny; ACMG nie zaleca zmian postepowania na podstawie MTHFR C677T"
            echo "    → Osoby w ciazy / planujace ciaze NIE powinny odstawiac kwasu foliowego; L-metylofolian jest opcja, nie koniecznoscia"
            echo "    → Wazne: B12 jako metylokobalamina"
            ;;
    esac
fi

# BCMO1 - witamina A z beta-karotenu
GT=$(check_snp "16:81264597")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  BCMO1 rs12934922: Wariant → Slaba konwersja beta-karotenu"
    echo "    → Jedz retinol (watrobka, jaja) zamiast marchewki"
fi

# VDR - witamina D
GT=$(check_snp "12:47846052")
if [[ "$GT" == *"1"* ]]; then
    echo "  VDR rs1544410: Wariant → Wieksza potrzeba witaminy D"
    echo "    → Suplementacja 2000-4000 IU dziennie"
fi

# GC - bialko wiazace wit. D
GT=$(check_snp "4:72618334")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  GC rs2282679: Wariant → Nizszy poziom wit. D w surowicy"
fi

# FADS1 - omega-3
GT=$(check_snp "11:61567029")  # TODO(build): unverified GRCh38
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
    echo "  FADS1 rs174546: T/T → Slaba konwersja ALA do EPA/DHA"
    echo "    → Jedz ryby lub suplementuj omega-3 (EPA/DHA bezposrednio)"
fi

# FUT2 - witamina B12
GT=$(check_snp "19:49206674")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  FUT2 rs602662: Non-secretor → Nizsze wchłanianie B12"
    echo "    → Rozważ suplementacje B12"
fi

# SLC23A1 - witamina C
GT=$(check_snp "5:139387768")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  SLC23A1 rs33972313: Wariant → Nizszy poziom wit. C"
fi

echo ""
echo "📊 Kofeina i uzywki:"
echo "─────────────────────────────────────────"

# CYP1A2 - kofeina
GT=$(check_snp "15:74749576")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  CYP1A2 rs762551: A/A → ☕ SZYBKI metabolizer kofeiny"
            echo "    → Kawa moze byc korzystna dla serca!"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  CYP1A2 rs762551: A/C → Sredni metabolizer"
            echo "    → Max 2-3 kawy dziennie"
            ;;
        "1/1"|"1|1")
            echo "  CYP1A2 rs762551: C/C → 🐢 WOLNY metabolizer kofeiny"
            echo "    → Ogranicz kawe! Ryzyko sercowe przy >2 filizankach"
            echo "    → Unikaj kofeiny po 14:00"
            ;;
    esac
fi

# ADORA2A - wrazliwosc na kofeine
GT=$(check_snp "22:24441816")  # TODO(build): unverified GRCh38
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
    echo "  ADORA2A rs5751876: T/T → Wrazliwosc na kofeine (lęk, bezsennosc)"
fi

# ADH1B i ALDH2 - alkohol (rs1229984: T = His48 = *2 szybki, C = Arg48 = *1)
# TODO(build): unverified GRCh38 orientation (table row is 'unverified')
RES=$(hx_gt_allele 4:99318162 T C rs1229984)
DOS=${RES%%$'\t'*}
ADH1B_FAST=0
case "$DOS" in
    0)
        ADH1B_FAST=1
        echo "  ADH1B rs1229984: T/T (*2/*2, His48) → Szybki metabolizm alkoholu"
        ;;
    1)
        # shellcheck disable=SC2034  # informational flag, not used in scoring yet
        ADH1B_FAST=1
        echo "  ADH1B rs1229984: T/C (*1/*2) → Przyspieszony metabolizm alkoholu (jeden allel His48)"
        ;;
    2)
        echo "  ADH1B rs1229984: C/C (*1/*1, Arg48) → Standardowy (wolniejszy) metabolizm alkoholu"
        ;;
esac

GT=$(check_snp "12:111803962")
if [[ "$GT" == *"1"* ]]; then
    echo "  ALDH2 rs671: ⚠️ Asian flush - nietolerancja alkoholu!"
    echo "    → Alkohol zwieksza ryzyko raka przelku - unikaj!"
fi

echo ""

#===============================================================================
# SEKCJA 6: COVID-19 I ODPORNOSC
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🦠 SEKCJA 6: COVID-19 I ODPORNOSC                                        ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# ACE2 - receptor wirusa
echo "📊 Podatnosc na COVID-19:"
echo "─────────────────────────────────────────"

# 3p21.31 locus - najsilniejszy GWAS hit dla ciezkiego COVID
GT=$(check_snp "3:45823240")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  rs11385942 (3p21.31): Wariant → Podwyzszone ryzyko ciezkiego przebiegu!"
    echo "    → Ten region pochodzi od neandertalczykow"
fi

# ABO - grupa krwi
GT=$(check_snp "9:136137106")  # TODO(build): unverified GRCh38
echo "  Grupa krwi ABO: (wymaga pelnej analizy)"
echo "    → Grupa 0: nieco nizsze ryzyko"
echo "    → Grupa A: nieco wyzsze ryzyko"

# IFNAR2 - receptor interferonu
GT=$(check_snp "21:34603203")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  IFNAR2 rs2236757: Wariant → Slabsza odpowiedz interferonowa"
fi

# OAS1 - odpowiedz antywirusowa (od neandertalczykow!)
GT=$(check_snp "12:113357193")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  OAS1 rs10774671: ✅ Wariant OCHRONNY (neandertalski!)"
    echo "    → Lepsza odpowiedz antywirusowa"
fi

# TYK2 - ciężkosc przebiegu
GT=$(check_snp "19:10463118")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  TYK2 rs74956615: Wariant ryzyka ciezkiego przebiegu"
fi

echo ""
echo "📊 Ogolna odpornosc:"
echo "─────────────────────────────────────────"

# IL6 - zapalenie
GT=$(check_snp "7:22727026")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  IL6 rs1800795: G/G → Nizsza odpowiedz zapalna"
            ;;
        "1/1"|"1|1")
            echo "  IL6 rs1800795: C/C → Silniejsza odpowiedz zapalna"
            echo "    → Moze byc dobra na infekcje, ale więkze ryzyko burzy cytokinowej"
            ;;
    esac
fi

# TLR4 - odpornosc wrodzona
GT=$(check_snp "9:120475302")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  TLR4 rs4986790: Wariant → Zmieniona odpowiedz na bakterie gram-ujemne"
fi

echo ""

#===============================================================================
# SEKCJA 7: ZMYSLY - SMAK, WECH, WZROK, SLUCH
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  👃 SEKCJA 7: ZMYSLY - SMAK, WECH, WZROK, SLUCH                           ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "📊 Smak:"
echo "─────────────────────────────────────────"

# TAS2R38 - gorycz (PROP/PTC). rs713598 A49P, C>G: gen na nici minus, wiec
# plus-strand C = Ala49 (AVI, non-taster), G = Pro49 (PAV, taster);
# DOS = liczba alleli G = liczba alleli PAV.
# TODO(build): unverified GRCh38; orientation assumes table row C>G with REF C = Ala49 (AVI), TAS2R38 minus strand
RES=$(hx_gt_allele 7:141972804 C G rs713598)
DOS=${RES%%$'\t'*}
if [ -n "$DOS" ]; then
    case $DOS in
        2)
            echo "  TAS2R38 rs713598: PAV/PAV → SUPERSMAKOWANIE!"
            echo "    → Bardzo wrazliwy na gorycz (brokuly, piwo, grejpfrut)"
            ;;
        1)
            echo "  TAS2R38 rs713598: PAV/AVI → Srednia wrazliwosc na gorycz"
            ;;
        0)
            echo "  TAS2R38 rs713598: AVI/AVI → Nie czujesz goryczki PROP/PTC"
            echo "    → Brokuly i piwo smakuja lagodniej"
            ;;
    esac
fi

# OR6A2 - kolendra jak mydlo?
GT=$(check_snp "11:57820524")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  OR6A2 rs72921001: Wariant → KOLENDRA SMAKUJE JAK MYDLO! 🧼"
else
    echo "  OR6A2 rs72921001: Kolendra smakuje normalnie"
fi

# TAS1R2 - slodki smak
GT=$(check_snp "1:18924070")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  TAS1R2 rs35874116: Wariant → Zmniejszona wrazliwosc na slodycz"
    echo "    → Mozesz potrzebowac wiecej cukru aby cos bylo 'slodkie'"
fi

echo ""
echo "📊 Wech:"
echo "─────────────────────────────────────────"

# OR7D4 - androstenon (zapach dzika/pizmowy)
GT=$(check_snp "19:9350580")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  OR7D4 rs61729907: Wariant → Nie czujesz zapachu androstenonu"
    echo "    → Wieprzowina moze pachnieć inaczej niz dla innych"
fi

# OR2M7 - aspiragus (mocz po szparagach)
GT=$(check_snp "1:248441820")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  OR2M7: Wariant → Nie czujesz zapachu moczu po szparagach"
fi

# Anosmia na izonitryle (zapach rybi)
GT=$(check_snp "11:57761519")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  TAAR5 region: Wariant receptora wechowego"
fi

echo ""
echo "📊 Wzrok:"
echo "─────────────────────────────────────────"

# GJD2 - krotkowzrocznosc
GT=$(check_snp "15:35005886")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  GJD2 rs634990: Wariant → Predyspozycja do krotkowzrocznosci"
fi

# RASGRF1 - krotkowzrocznosc
GT=$(check_snp "15:79378016")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  RASGRF1 rs8027411: Wariant → Wyzsże ryzyko krotkowzrocznosci"
    echo "    → Czas na swiežym powietrzu w dziecinstwie ochronny!"
fi

# OPN1LW/OPN1MW - widzenie barw
echo "  Widzenie barw: (geny na chr X - specjalna analiza potrzebna)"

echo ""
echo "📊 Sluch:"
echo "─────────────────────────────────────────"

# GRM7 - presbyacusis (utrata sluchu z wiekiem)
GT=$(check_snp "3:7644503")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  GRM7 rs11928865: Wariant → Wieksze ryzyko utraty sluchu z wiekiem"
    echo "    → Chron sluch przed halas em!"
fi

# DFNA5/GSDME
GT=$(check_snp "7:24771936")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  GSDME rs2802253: Wariant zwiazany ze sluchem"
fi

echo ""

#===============================================================================
# SEKCJA 8: SKORA, WLOSY, STARZENIE
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  💇 SEKCJA 8: SKORA, WLOSY, STARZENIE                                     ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "📊 Lysienie (mezczyzni):"
echo "─────────────────────────────────────────"

BALD_SCORE=0

# AR - receptor androgenowy (glowny gen!)
GT=$(check_snp "X:67723340")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  AR rs6152: Wariant → PODWYZSZONE ryzyko lysienia"
    BALD_SCORE=$((BALD_SCORE + 3))
fi

# 20p11 - drugi najsilniejszy locus
GT=$(check_snp "20:21995360")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  rs1160312 (20p11): Wariant lysienia"
    BALD_SCORE=$((BALD_SCORE + 2))
fi

# 1p36
GT=$(check_snp "1:9961401")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  rs12565727 (1p36): Wariant lysienia"
    BALD_SCORE=$((BALD_SCORE + 1))
fi

# PAX1
GT=$(check_snp "20:21696679")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  rs2180439 (PAX1): Wariant lysienia"
    BALD_SCORE=$((BALD_SCORE + 1))
fi

if [ $BALD_SCORE -eq 0 ]; then
    echo "  → Niskie genetyczne ryzyko lysienia ✅"
elif [ $BALD_SCORE -le 3 ]; then
    echo "  → Umiarkowane ryzyko lysienia 🟡"
else
    echo "  → Wysokie genetyczne ryzyko lysienia 🔴"
    echo "    → Finasteryd/minoxidil moga byc skuteczne"
fi

echo ""
echo "📊 Siwienie:"
echo "─────────────────────────────────────────"

# IRF4 - siwienie
GT=$(check_snp "6:396321")
if [[ "$GT" == *"1"* ]]; then
    echo "  IRF4 rs12203592: Wariant → Wczesniejsze siwienie"
fi

echo ""
echo "📊 Starzenie skory:"
echo "─────────────────────────────────────────"

# MMP1 - kolagen
GT=$(check_snp "11:102660769")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  MMP1 rs1799750: 2G allel → Szybsza degradacja kolagenu"
    echo "    → Retinol, wit. C, SPF 30+ codziennie!"
fi

# COL1A1 - jakosc kolagenu
GT=$(check_snp "17:50201632")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  COL1A1 rs1800012: Wariant → Slabszy kolagen"
fi

# MC1R - UV i pieg i
GT=$(check_snp "16:89919436")
if [[ "$GT" == *"1"* ]]; then
    echo "  MC1R rs1805007: Wariant rudy → Wyzsza wrazliwosc na UV!"
    echo "    → SPF 50+, unikaj slonca 11:00-15:00"
fi

echo ""
echo "📊 Woskowina uszna i zapach ciala:"
echo "─────────────────────────────────────────"

# ABCC11
GT=$(check_snp "16:48224287")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ABCC11 rs17822931: C/C → Mokra woskowina, normalny zapach ciala"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  ABCC11 rs17822931: C/T → Mokra woskowina (allel mokry jest dominujacy)"
            ;;
        "1/1"|"1|1")
            echo "  ABCC11 rs17822931: T/T → SUCHA woskowina, MNIEJSZY zapach ciala"
            echo "    → Dezodorant moze byc zbedny! (typowe w Azji)"
            ;;
    esac
fi

echo ""

#===============================================================================
# SEKCJA 9: DLUGOWIECZNOSC
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🎂 SEKCJA 9: DLUGOWIECZNOSC                                              ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

LONGEVITY_SCORE=0

# FOXO3 - glowny gen dlugowiecznosci!
GT=$(check_snp "6:108989256")  # TODO(build): unverified GRCh38
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  🌟 FOXO3 rs2802292: G/G → KORZYSTNY wariant dlugowiecznosci!"
            echo "    → Zwiazany z dozyciem 100 lat"
            echo "    → Lepsza odpowiedz na stres oksydacyjny"
            LONGEVITY_SCORE=$((LONGEVITY_SCORE + 3))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  FOXO3 rs2802292: G/T → Jeden korzystny allel"
            LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
            ;;
        "1/1"|"1|1")
            echo "  FOXO3 rs2802292: T/T → Standardowy wariant"
            ;;
    esac
fi

# APOE (wczesniej analizowane, ale kluczowe)
APOE1=$(check_snp "19:44908684")
APOE2=$(check_snp "19:44908822")

if [[ "$APOE1" != *"1"* ]] && [[ "$APOE2" != *"1"* ]]; then
    echo "  APOE e3/e3: ✅ OPTYMALNY dla dlugowiecznosci"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 2))
elif [[ "$APOE2" == *"1"* ]] && [[ "$APOE1" != *"1"* ]]; then
    echo "  APOE e2 obecny: ✅ KORZYSTNY (ale uwaga na lipidy)"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 2))
elif [[ "$APOE1" == *"1"* ]]; then
    echo "  APOE e4 obecny: ⚠️ Zwiazany z krotszym zyciem"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE - 2))
fi

# CETP - HDL
GT=$(check_snp "16:57015091")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  CETP rs5882: Wariant → Wyzszy HDL (dobry cholesterol)"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
fi

# TERT - telomeraza
GT=$(check_snp "5:1287194")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  TERT rs2736100: Wariant → Potencjalnie dluzsze telomery"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
fi

# IL6 - niskie zapalenie korzystne
GT=$(check_snp "7:22727026")
if [[ "$GT" == "0/0" ]] || [[ "$GT" == "0|0" ]]; then
    echo "  IL6 G/G: Niski poziom zapalenia → Korzystne dla dlugowiecznosci"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
fi

# SIRT3 - sirtuiny
GT=$(check_snp "11:236236")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  SIRT3: Wariant → Lepszy metabolizm mitochondrialny"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "WYNIK DLUGOWIECZNOSCI: $LONGEVITY_SCORE / 9 punktow"
echo ""
if [ $LONGEVITY_SCORE -ge 6 ]; then
    echo "🏆 BARDZO KORZYSTNY profil genetyczny dla dlugiego zycia!"
elif [ $LONGEVITY_SCORE -ge 3 ]; then
    echo "✅ DOBRY profil - zdrowy styl zycia moze wiele dodac"
elif [ $LONGEVITY_SCORE -ge 0 ]; then
    echo "🟡 STANDARDOWY profil - styl zycia ma kluczowe znaczenie"
else
    echo "⚠️ Niektore warianty ryzyka - styl zycia jeszcze WAZNIEJSZY"
fi
echo ""

#===============================================================================
# SEKCJA 10: SPERSONALIZOWANA SUPLEMENTACJA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  💊 SEKCJA 10: SPERSONALIZOWANY PLAN SUPLEMENTACJI                        ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Na podstawie Twoich wariantow genetycznych:"
echo ""

# MTHFR
GT=$(check_snp "1:11796321")
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ L-METYLOFOLIAN (opcja, nie koniecznosc - kwas foliowy pozostaje skuteczny)"
    echo "    Dawka: 400-800 mcg/dzien"
    echo "    Powod: MTHFR - obnizony metabolizm folianu"
    echo ""
fi

# VDR
GT=$(check_snp "12:47846052")
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ WITAMINA D3"
    echo "    Dawka: 2000-4000 IU/dzien (zalezenie od poziomu)"
    echo "    Powod: VDR - wieksze zapotrzebowanie na wit. D"
    echo ""
fi

# FADS1
GT=$(check_snp "11:61567029")  # TODO(build): unverified GRCh38
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
    echo "  ✓ OMEGA-3 (EPA + DHA)"
    echo "    Dawka: 1-2g EPA+DHA/dzien"
    echo "    Powod: FADS1 - slaba konwersja z ALA"
    echo ""
fi

# BCMO1
GT=$(check_snp "16:81264597")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ WITAMINA A (retinol, nie beta-karoten)"
    echo "    Dawka: 2500-5000 IU/dzien lub z diety (watrobka)"
    echo "    Powod: BCMO1 - slaba konwersja beta-karotenu"
    echo ""
fi

# FUT2
GT=$(check_snp "19:49206674")  # TODO(build): unverified GRCh38
if [[ "$GT" == *"1"* ]]; then
    echo "  ✓ WITAMINA B12 (metylokobalamina)"
    echo "    Dawka: 500-1000 mcg/dzien"
    echo "    Powod: FUT2 - nizsze wchlanianie B12"
    echo ""
fi

# SOD2 + ogolne antyoksydanty
GT=$(check_snp "6:159692840")
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "0/0" ]]; then
    echo "  ✓ Rozważ: NAC lub GLUTATION"
    echo "    Dawka: 600-1200 mg NAC/dzien"
    echo "    Powod: SOD2 homozygota - wsparcie antyoksydacyjne"
    echo ""
fi

# COMT Met/Met
GT=$(check_snp "22:19963748")
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
    echo "  ✓ MAGNEZ (cytrynian lub glicynian)"
    echo "    Dawka: 200-400 mg/dzien"
    echo "    Powod: COMT Met/Met - wsparcie dla układu nerwowego"
    echo ""
fi

# CYP1A2 wolny
GT=$(check_snp "15:74749576")
if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
    echo "  ⚠️ UWAGA: OGRANICZ KOFEINE"
    echo "    Max: 1-2 kawy/dzien, nie po 14:00"
    echo "    Powod: CYP1A2 wolny metabolizer"
    echo ""
fi

echo ""
echo "UWAGA: To sugestie na podstawie genow. Skonsultuj z lekarzem!"
echo ""

#===============================================================================
# SEKCJA 11: PODSUMOWANIE MEGA ANALIZY
#===============================================================================

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                       ║"
echo "║                     P O D S U M O W A N I E   M E G A   A N A L I Z Y                 ║"
echo "║                                                                                       ║"
echo "╚═══════════════════════════════════════════════════════════════════════════════════════╝"
echo ""

# Zapisz raport do pliku
REPORT="$OUTDIR/MEGA_RAPORT_$(date +%Y%m%d_%H%M%S).txt"

cat << EOF > "$REPORT"
================================================================================
                    MEGA ANALIZA GENETYCZNA - RAPORT
                    Wygenerowano: $(date)
                    VCF: $VCF
================================================================================

PODSUMOWANIE KLUCZOWYCH WYNIKOW:
--------------------------------

1. DNA NEANDERTALCZYKA: $NEANDER_SCORE markerow archaicznych

2. NOSICIELSTWO CHOROB: $CARRIER_COUNT wariantow wykrytych

3. CHRONOTYP: $([ $MORNING_SCORE -gt $EVENING_SCORE ] && echo "SKOWRONEK" || ([ $EVENING_SCORE -gt $MORNING_SCORE ] && echo "SOWA" || echo "POSREDNI"))

4. DLUGOWIECZNOSC: $LONGEVITY_SCORE / 9 punktow

5. LYSIENIE (genetyczne): $BALD_SCORE punktow ryzyka

================================================================================

Ten raport ma charakter informacyjny i edukacyjny.
Nie zastepuje konsultacji z lekarzem genetykiem.

Wygenerowano przez: Helixight Genetic Analysis Platform
================================================================================
EOF

echo "TWOJ GENETYCZNY PROFIL W PIGULCE:"
echo "─────────────────────────────────────────────────────────────────────────────"
echo ""
echo "  🦴 DNA Neandertalczyka:   $NEANDER_SCORE markerow"
echo "  🧬 Nosicielstwo chorob:   $CARRIER_COUNT wariantow"
echo "  🌙 Chronotyp:             $([ $MORNING_SCORE -gt $EVENING_SCORE ] && echo "🐦 SKOWRONEK" || ([ $EVENING_SCORE -gt $MORNING_SCORE ] && echo "🦉 SOWA" || echo "🕐 POSREDNI"))"
echo "  🎂 Dlugowiecznosc:        $LONGEVITY_SCORE / 9 pkt"
echo "  💇 Lysienie (genetycznie): $BALD_SCORE pkt ryzyka"
echo ""
echo "─────────────────────────────────────────────────────────────────────────────"
echo ""
echo "Pelny raport zapisany: $REPORT"
echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  Dziekujemy za uzycie Helixight Mega Analizy!                             ║"
echo "║  Pamietaj: Geny to nie przeznaczenie - styl zycia ma ogromne znaczenie!  ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
hx_report_coverage
