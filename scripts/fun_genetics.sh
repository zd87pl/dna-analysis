#!/bin/bash
#===============================================================================
# ZAAWANSOWANA ANALIZA GENETYCZNA - CIEKAWE ODKRYCIA
# Rzeczy ktorych jeszcze nie sprawdzilismy!
#===============================================================================

# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"
set -eo pipefail

VCF="${1:-}"
hx_require_arg "$VCF" "$(basename "$0") <variants.vcf.gz>"

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║           ZAAWANSOWANA ANALIZA GENETYCZNA                                 ║"
echo "║           Ciekawe odkrycia z Twojego DNA                                  ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Wspolrzedne GRCh38 (scripts/lib/GRCh38_COORDINATES.md); biblioteka normalizuje
# nazwy kontigow i buduje zindeksowana kopie robocza. hx_gt zwraca genotyp
# pierwszego SNP (faza znormalizowana, brak wywolania = pusty).
hx_prepare_vcf "$VCF"

#===============================================================================
# 1. DNA NEANDERTALCZYKA
#===============================================================================

echo "🦴 DNA NEANDERTALCZYKA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Markery neandertalskie (im wiecej wariantow, tym wiecej DNA neandertalczyka):"
echo ""

# Znane warianty neandertalskie
NEANDER_COUNT=0

# rs2066807 - BNC2 (jasna skora, od neandertalczykow)
result=$(hx_gt 9:16817033)  # TODO(build): unverified GRCh38
if [ "$result" == "0/1" ] || [ "$result" == "1/1" ]; then
    echo "  ✓ rs2066807 (BNC2): Wariant neandertalski - jasna skora"
    NEANDER_COUNT=$((NEANDER_COUNT + 1))
fi

# rs10490924 - zwiazany z neandertalczykami
result=$(hx_gt 10:124214448)  # TODO(build): unverified GRCh38
if [ "$result" == "0/1" ] || [ "$result" == "1/1" ]; then
    echo "  ✓ rs10490924: Wariant neandertalski (zwiazany z AMD)"
    NEANDER_COUNT=$((NEANDER_COUNT + 1))
fi

# rs3917862 - SELP (krzepniecie krwi, od neandertalczykow)
result=$(hx_gt 1:169573470)  # TODO(build): unverified GRCh38
if [ "$result" == "0/1" ] || [ "$result" == "1/1" ]; then
    echo "  ✓ rs3917862 (SELP): Wariant neandertalski - krzepniecie krwi"
    NEANDER_COUNT=$((NEANDER_COUNT + 1))
fi

# rs1042602 - TYR (pigmentacja)
result=$(hx_gt 11:88911696)  # TODO(build): unverified GRCh38
if [ "$result" == "0/1" ] || [ "$result" == "1/1" ]; then
    echo "  ✓ rs1042602 (TYR): Wariant pigmentacji (czesc neandertalska)"
    NEANDER_COUNT=$((NEANDER_COUNT + 1))
fi

# rs4988235 - LCT (tolerancja laktozy - NIE neandertalska, ale ciekawa)
# Neandertalczycy NIE mieli tolerancji laktozy

echo ""
echo "Znalezione markery neandertalskie: $NEANDER_COUNT"
echo "(Typowy Europejczyk: 1-4% DNA neandertalskiego, ~10-20 znanych markerow)"
echo ""

#===============================================================================
# 2. CHRONOTYP - SOWA CZY SKOWRONEK?
#===============================================================================

echo "🌙 CHRONOTYP - SOWA CZY SKOWRONEK?"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

MORNING_SCORE=0
EVENING_SCORE=0

# rs12914385 - PER2 (zegar biologiczny)
result=$(hx_gt 2:239186871)  # TODO(build): unverified GRCh38
if [ "$result" == "1/1" ]; then
    echo "  rs12914385 (PER2): T/T → Tendencja do wczesnego wstawania"
    MORNING_SCORE=$((MORNING_SCORE + 2))
elif [ "$result" == "0/1" ]; then
    MORNING_SCORE=$((MORNING_SCORE + 1))
fi

# rs1801260 - CLOCK (glowny gen zegara)
result=$(hx_gt 4:56294068)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs1801260 (CLOCK): Wariant → Tendencja do poznego chodzenia spac"
    EVENING_SCORE=$((EVENING_SCORE + 2))
fi

# rs2304672 - PER3 (dlugosc snu)
result=$(hx_gt 1:7870203)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs2304672 (PER3): Wariant → Krotszy sen, wieksza aktywnosc wieczorem"
    EVENING_SCORE=$((EVENING_SCORE + 1))
fi

# rs11121022 - ARNTL (regulacja rytmu)
result=$(hx_gt 11:13387651)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs11121022 (ARNTL): Wariant → Wplywa na rytm dobowy"
fi

echo ""
if [ $MORNING_SCORE -gt $EVENING_SCORE ]; then
    echo "Wynik: 🐦 Prawdopodobnie SKOWRONEK (ranny ptaszek)"
elif [ $EVENING_SCORE -gt $MORNING_SCORE ]; then
    echo "Wynik: 🦉 Prawdopodobnie SOWA (nocny marek)"
else
    echo "Wynik: 🕐 Typ posredni"
fi
echo ""

#===============================================================================
# 3. EMPATIA I ZACHOWANIA SPOLECZNE
#===============================================================================

echo "🧠 EMPATIA I ZACHOWANIA SPOLECZNE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# rs53576 - OXTR (receptor oksytocyny - "hormonu milosci")
GT=$(hx_gt 3:8804371)  # TODO(build): unverified GRCh38

if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  rs53576 (OXTR): G/G → Wysoka empatia, latwiejsze nawiazywanie wiezji"
            echo "    • Lepsze rozpoznawanie emocji u innych"
            echo "    • Wieksze wsparcie spoleczne"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  rs53576 (OXTR): G/A → Sredni poziom empatii"
            ;;
        "1/1"|"1|1")
            echo "  rs53576 (OXTR): A/A → Nizsza empatia, wieksza niezaleznosc"
            echo "    • Mniejsza potrzeba wsparcia spolecznego"
            echo "    • Lepsza odpornosc na stres spoleczny"
            ;;
    esac
fi
echo ""

# rs1800497 - DRD2 (receptor dopaminy - Taq1A)
result=$(hx_gt 11:113400106)  # TODO(build): unverified GRCh38
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs1800497 (DRD2 Taq1A): Normalny - standardowa odpowiedz na nagrody"
            ;;
        "0/1"|"1/1")
            echo "  rs1800497 (DRD2 Taq1A): Wariant A1 → Mniej receptorow dopaminy"
            echo "    • Wieksza potrzeba stymulacji"
            echo "    • Potencjalnie wieksze ryzyko uzaleznien"
            ;;
    esac
fi
echo ""

#===============================================================================
# 4. METABOLIZM I DIETA
#===============================================================================

echo "🥗 NUTRIGENOMIKA - DIETA DOPASOWANA DO GENOW"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# FTO - otylosc
result=$(hx_gt 16:53786615)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs9939609 (FTO): T/T → Niskie ryzyko otylosci"
            ;;
        "0/1"|"0|1")
            echo "  rs9939609 (FTO): T/A → Srednie ryzyko otylosci (+1.5 kg sredniej masy)"
            ;;
        "1/1"|"1|1")
            echo "  rs9939609 (FTO): A/A → Podwyzszone ryzyko otylosci (+3 kg sredniej masy)"
            echo "    → Rekomendacja: Wieksza aktywnosc fizyczna niweluje ten efekt!"
            ;;
    esac
fi

# APOA2 - nasycone tluszcze
result=$(hx_gt 1:161222292)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs5082 (APOA2): Wariant → Wrazliwosc na tluszcze nasycone"
    echo "    → Rekomendacja: Ogranicz tluszcze nasycone (<22g/dzien)"
fi

# TCF7L2 - cukrzyca
result=$(hx_gt 10:112998590)
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs7903146 (TCF7L2): Wariant → Podwyzszone ryzyko cukrzycy typu 2"
    echo "    → Rekomendacja: Dieta niskoglikemiczna, regularne posilki"
fi

# FADS1 - omega-3
result=$(hx_gt 11:61567029)  # TODO(build): unverified GRCh38
if [ -n "$result" ]; then
    case $result in
        "1/1"|"1|1")
            echo "  rs174546 (FADS1): Wariant → Slaba konwersja ALA do EPA/DHA"
            echo "    → Rekomendacja: Jedz ryby lub suplementuj omega-3 (EPA/DHA)"
            ;;
    esac
fi

# BCMO1 - witamina A z beta-karotenu
result=$(hx_gt 16:81264597)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs12934922 (BCMO1): Wariant → Slaba konwersja beta-karotenu do wit. A"
    echo "    → Rekomendacja: Jedz zrodla retinolu (watrobka, jaja) zamiast marchewki"
fi

echo ""

#===============================================================================
# 5. STARZENIE SIE I DLUGOWIECZNOSC
#===============================================================================

echo "🎂 DLUGOWIECZNOSC I STARZENIE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

LONGEVITY_SCORE=0

# FOXO3 - glowny gen dlugowiecznosci
result=$(hx_gt 6:108989256)  # TODO(build): unverified GRCh38
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs2802292 (FOXO3): G/G → KORZYSTNY wariant dlugowiecznosci!"
            echo "    • Zwiazany z dozycie 100 lat"
            echo "    • Lepsza odpowiedz na stres oksydacyjny"
            LONGEVITY_SCORE=$((LONGEVITY_SCORE + 2))
            ;;
        "0/1"|"0|1")
            echo "  rs2802292 (FOXO3): G/T → Jeden korzystny allel"
            LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
            ;;
        "1/1"|"1|1")
            echo "  rs2802292 (FOXO3): T/T → Standardowy wariant"
            ;;
    esac
fi

# CETP - cholesterol HDL
result=$(hx_gt 16:57015091)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs5882 (CETP): Wariant → Wyzszy HDL (dobry cholesterol)"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
fi

# APOE (juz wczesniej sprawdzane, ale podsumujmy)
APOE1=$(hx_gt 19:44908684)
APOE2=$(hx_gt 19:44908822)
if [[ "$APOE1" == "0/0"* ]] && [[ "$APOE2" == "0/0"* ]]; then
    echo "  APOE: e3/e3 → OPTYMALNY dla dlugowiecznosci"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 2))
elif [[ "$APOE2" == *"1"* ]]; then
    echo "  APOE: Obecny e2 → KORZYSTNY dla dlugowiecznosci"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
fi

# TERT - telomeraza
result=$(hx_gt 5:1287194)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs2736100 (TERT): Wariant → Potencjalnie dluzsze telomery"
    LONGEVITY_SCORE=$((LONGEVITY_SCORE + 1))
fi

echo ""
echo "Wynik dlugowiecznosci: $LONGEVITY_SCORE / 7"
if [ $LONGEVITY_SCORE -ge 5 ]; then
    echo "→ Bardzo korzystny profil genetyczny dla dlugiego zycia!"
elif [ $LONGEVITY_SCORE -ge 3 ]; then
    echo "→ Dobry profil genetyczny"
else
    echo "→ Standardowy profil (styl zycia ma wieksże znaczenie!)"
fi
echo ""

#===============================================================================
# 6. WRAZLIWOSC NA BOL
#===============================================================================

echo "💉 WRAZLIWOSC NA BOL"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# COMT (juz mamy z poprzedniej analizy)
result=$(hx_gt 22:19963748)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs4680 (COMT): Val/Val → Nizsza wrazliwosc na bol"
            echo "    • Mniejsze zapotrzebowanie na opioidy po operacjach"
            ;;
        "0/1"|"0|1")
            echo "  rs4680 (COMT): Val/Met → Srednia wrazliwosc na bol"
            ;;
        "1/1"|"1|1")
            echo "  rs4680 (COMT): Met/Met → Wyzsza wrazliwosc na bol"
            echo "    • Wieksze zapotrzebowanie na leki przeciwbolowe"
            ;;
    esac
fi

# OPRM1 - receptor opioidowy
result=$(hx_gt 6:154039662)
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs1799971 (OPRM1): Wariant G → Wieksza potrzeba opioidow"
    echo "    • Mniejsza odpowiedz na morfine"
fi

# SCN9A - kanal sodowy
result=$(hx_gt 2:167138656)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs6746030 (SCN9A): Wariant → Zmieniona percepcja bolu"
fi

echo ""

#===============================================================================
# 7. ALERGIE I UKLAD ODPORNOSCIOWY
#===============================================================================

echo "🤧 ALERGIE I ODPORNOSC"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# IL13 - alergie
result=$(hx_gt 5:131995964)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs20541 (IL13): Wariant → Podwyzszone ryzyko alergii i astmy"
fi

# FLG - filagryna (atopowe zapalenie skory)
result=$(hx_gt 1:152285861)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs61816761 (FLG): Wariant → Ryzyko AZS (atopowe zapalenie skory)"
fi

# HLA-DQ (juz sprawdzane - celiakia)
# Podsumowanie

echo ""

#===============================================================================
# 8. SKORA I STARZENIE SKORY
#===============================================================================

echo "✨ SKORA I STARZENIE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# MMP1 - kolagen (zmarszczki)
result=$(hx_gt 11:102660769)  # TODO(build): unverified GRCh38
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs1799750 (MMP1): 1G/1G → Wolniejsze starzenie skory"
            ;;
        "0/1"|"1/1")
            echo "  rs1799750 (MMP1): 2G allel → Szybsza degradacja kolagenu"
            echo "    → Rekomendacja: Krem z retinolem, wit. C, ochrona UV"
            ;;
    esac
fi

# STXBP5L - cellulit
result=$(hx_gt 3:120819158)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs7636981 (STXBP5L): Wariant → Wieksze ryzyko cellulitu (kobiety)"
fi

# IRF4 - piegi/slonce
result=$(hx_gt 6:396321)
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs12203592 (IRF4): Wariant → Wieksze ryzyko uszkodzen slonecznych"
    echo "    → Rekomendacja: SPF 30+ codziennie!"
fi

echo ""

#===============================================================================
# 9. WZROK
#===============================================================================

echo "👓 WZROK"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# GJD2 - krotkowzrocznosc
result=$(hx_gt 15:35005886)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs634990 (GJD2): Wariant → Predyspozycja do krotkowzrocznosci"
fi

# RASGRF1 - krotkowzrocznosc
result=$(hx_gt 15:79378016)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs8027411 (RASGRF1): Wariant → Zwiekszone ryzyko krotkowzrocznosci"
    echo "    → Rekomendacja: Czas na swiezym powietrzu zmniejsza ryzyko!"
fi

echo ""

#===============================================================================
# 10. REAKCJE NA SUBSTANCJE
#===============================================================================

echo "🧪 REAKCJE NA ROZNE SUBSTANCJE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# OR6A2 - kolendra (smakuje jak mydlo?)
result=$(hx_gt 11:57820524)  # TODO(build): unverified GRCh38
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs72921001 (OR6A2): A/A → Kolendra smakuje normalnie"
            ;;
        "0/1"|"1/1")
            echo "  rs72921001 (OR6A2): Wariant C → Kolendra moze smakowac jak MYDLO!"
            ;;
    esac
fi

# CYP1A2 - kofeina (juz sprawdzane)
# ADH1B - alkohol (juz sprawdzane)

# ABCC11 - woskowina uszna / zapach ciala
result=$(hx_gt 16:48224287)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs17822931 (ABCC11): C/C → Mokra woskowina, normalny zapach ciala"
            ;;
        "0/1")
            echo "  rs17822931 (ABCC11): C/T → Mokra woskowina (allel mokry jest dominujacy)"
            ;;
        "1/1"|"1|1")
            echo "  rs17822931 (ABCC11): T/T → Sucha woskowina, mniejszy zapach ciala"
            echo "    • Ten wariant sprawia, ze dezodorant jest mniej potrzebny!"
            ;;
    esac
fi

# ACTN3 - sprint vs wytrzymalosc (juz sprawdzane)

echo ""

#===============================================================================
# 11. UNIKALNE ZDOLNOSCI
#===============================================================================

echo "🦸 UNIKALNE ZDOLNOSCI GENETYCZNE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# BDNF - plastycznosc mozgu, uczenie sie
result=$(hx_gt 11:27658369)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs6265 (BDNF): Val/Val → Optymalna plasycznosc mozgu"
            echo "    • Lepsze uczenie sie nowych umiejetnosci motorycznych"
            ;;
        "0/1"|"1/1")
            echo "  rs6265 (BDNF): Met carrier → Inna plastycznosc mozgu"
            echo "    • Lepsza pamiec epizodyczna"
            echo "    • Cwiczenia fizyczne szczegolnie korzystne!"
            ;;
    esac
fi

# KIBRA - pamiec
result=$(hx_gt 5:167791876)  # TODO(build): unverified GRCh38
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  rs17070145 (KIBRA): C/C → Standardowa pamiec"
            ;;
        "0/1"|"1/1")
            echo "  rs17070145 (KIBRA): T allel → LEPSZA pamiec epizodyczna!"
            echo "    • Lepsze zapamietywanie wydarzen i faktow"
            ;;
    esac
fi

# SNAP25 - IQ / zdolnosci poznawcze
result=$(hx_gt 20:10268206)  # TODO(build): unverified GRCh38
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  rs363050 (SNAP25): Wariant → Zwiazany z wyzszym IQ"
fi

echo ""

#===============================================================================
# PODSUMOWANIE
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                         ANALIZA ZAKONCZONA                                ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Te wyniki to tylko czesc Twojego genetycznego potencjalu!"
echo "Pamietaj: geny to nie przeznaczenie - styl zycia ma ogromny wplyw."
echo ""
hx_report_coverage
