#!/bin/bash
#===============================================================================
# QUICK GENOME INSIGHTS
# Natychmiastowe interesujące odkrycia z VCF
#
# Współrzędne: GRCh38, nazwy kontigów bez prefiksu "chr" (biblioteka
# scripts/lib/vcf_common.sh normalizuje plik wejściowy przed zapytaniami).
# Kanoniczna tabela: scripts/lib/GRCh38_COORDINATES.md
#===============================================================================

# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"
set -eo pipefail

VCF="${1:-}"
hx_require_arg "$VCF" "$(basename "$0") <variants.vcf.gz>"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║              QUICK GENOME INSIGHTS                                ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

hx_prepare_vcf "$VCF"

# gt_at REGION REF ALT [RSID]: znormalizowany genotyp rekordu, którego REF/ALT
# zgadzają się z tabelą kanoniczną (pusty, gdy brak rekordu lub inne allele).
# Genotyp wyprowadzany z dawki oczekiwanego allelu ALT (0/0, 0/1, 1/1), tak by
# rekordy wieloalleliczne (np. G>A,C z GT 0/2) nie były mylone z nosicielstwem.
gt_at() {
    local out d
    out=$(hx_gt_allele "$@")
    [[ -n $out ]] || { printf '\n'; return 0; }
    d=${out%%$'\t'*}
    case $d in
        0) printf '0/0\n' ;;
        1) printf '0/1\n' ;;
        2) printf '1/1\n' ;;
        *) printf '%s\n' "${out#*$'\t'}" ;;
    esac
}

# gt_letters GT REF ALT: "0/1" + C/T -> "C/T" (litery z nici plus, jak w VCF)
gt_letters() {
    local a=${1%%/*} b=${1#*/}
    [[ $a == 0 ]] && a=$2 || a=$3
    [[ $b == 0 ]] && b=$2 || b=$3
    printf '%s/%s' "$a" "$b"
}

#===============================================================================
# PODSTAWOWE STATYSTYKI
#===============================================================================

echo "📊 PODSTAWOWE STATYSTYKI"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

TOTAL=$(bcftools view -H "$HX_VCF" | wc -l | tr -d ' ')
SNPS=$(bcftools view -v snps -H "$HX_VCF" | wc -l | tr -d ' ')
INDELS=$(bcftools view -v indels -H "$HX_VCF" | wc -l | tr -d ' ')

printf "Wszystkie warianty: %'d\n" "$TOTAL"
printf "SNPs:               %'d\n" "$SNPS"
printf "Indele:             %'d\n" "$INDELS"
echo ""

#===============================================================================
# PŁEĆ GENETYCZNA
#===============================================================================

echo "👤 PŁEĆ GENETYCZNA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Błąd bcftools przerywa skrypt (pipefail); 0 oznacza tu prawdziwy brak rekordów.
X_COUNT=$(hx_query X '%POS\n' | wc -l | tr -d ' ')
Y_COUNT=$(hx_query Y '%POS\n' | wc -l | tr -d ' ')

echo "Warianty chr X: $X_COUNT"
echo "Warianty chr Y: $Y_COUNT"

if [ "$X_COUNT" -eq 0 ] && [ "$Y_COUNT" -eq 0 ]; then
    echo "Wniosek: brak danych o chr X/Y - nie można ustalić płci"
    SEX=""
elif [ "$Y_COUNT" -gt "500" ]; then
    echo "Wniosek: ♂ MĘŻCZYZNA (XY)"
    SEX="M"
else
    echo "Wniosek: ♀ KOBIETA (XX)"
    SEX="F"
fi
echo ""

#===============================================================================
# KOLOR OCZU
#===============================================================================

echo "👁️ KOLOR OCZU"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs12913832 - główny marker koloru oczu (HERC2, 15:28120472 A>G)
GT=$(gt_at 15:28120472 A G rs12913832)

if [ -n "$GT" ]; then
    case $GT in
        "0/0")
            echo "rs12913832: A/A → 🟤 Prawdopodobnie BRĄZOWE oczy"
            ;;
        "0/1")
            echo "rs12913832: A/G → 🟢 Prawdopodobnie ZIELONE lub PIWNE oczy"
            ;;
        "1/1")
            echo "rs12913832: G/G → 🔵 Prawdopodobnie NIEBIESKIE oczy"
            ;;
        *)
            echo "rs12913832: $GT (nietypowy genotyp)"
            ;;
    esac
else
    echo "Brak danych dla rs12913832"
fi
echo ""

#===============================================================================
# KOLOR WŁOSÓW
#===============================================================================

echo "💇 KOLOR WŁOSÓW / RYZYKO RUDOŚCI"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# MC1R - rude włosy: liczba alleli T rs1805007 (R151C, 16:89919436 GRCh38,
# jedyny marker MC1R w tabeli kanonicznej; stare pozycje 16:89985844-89986144
# były GRCh37).
MC1R_COUNT=0
GT=$(gt_at 16:89919436 C T rs1805007)
if [ -n "$GT" ]; then
    d=$(hx_alt_dosage "$GT"); [ "$d" != NA ] && MC1R_COUNT=$((MC1R_COUNT + d))
fi

if [ "$MC1R_COUNT" -gt "0" ]; then
    echo "Warianty MC1R: $MC1R_COUNT → 🔴 Możliwe RUDE włosy lub piegi"
else
    echo "Warianty MC1R: 0 → Brak predyspozycji do rudości"
fi

# KITLG - blond
BLOND=$(hx_gt 12:89328335)  # TODO(build): unverified GRCh38
if [ "$BLOND" == "0/1" ] || [ "$BLOND" == "1/1" ]; then
    echo "rs12821256: wariant → 👱 Skłonność do JASNYCH włosów"
fi
echo ""

#===============================================================================
# TOLERANCJA LAKTOZY
#===============================================================================

echo "🥛 TOLERANCJA LAKTOZY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs4988235 (LCT/MCM6 -13910C>T): 2:135851076 G>A na nici plus; A = tolerancja
GT=$(gt_at 2:135851076 G A rs4988235)

if [ -n "$GT" ]; then
    case $GT in
        "0/0")
            echo "rs4988235: C/C → ❌ NIETOLERANCJA laktozy (typ afrykański/azjatycki)"
            ;;
        "0/1")
            echo "rs4988235: C/T → ⚠️ CZĘŚCIOWA tolerancja"
            ;;
        "1/1")
            echo "rs4988235: T/T → ✅ PEŁNA tolerancja laktozy (typ europejski)"
            ;;
    esac
else
    echo "Brak danych"
fi
echo ""

#===============================================================================
# METABOLIZM KOFEINY
#===============================================================================

echo "☕ METABOLIZM KOFEINY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs762551 (CYP1A2): 15:74749576 A>C; A = *1F szybki metabolizer
CAF_REF=A; CAF_ALT=C
GT=$(gt_at 15:74749576 "$CAF_REF" "$CAF_ALT" rs762551)

if [ -n "$GT" ]; then
    LETTERS=$(gt_letters "$GT" "$CAF_REF" "$CAF_ALT")
    case $GT in
        "0/0")
            echo "rs762551: $LETTERS → ⚡ SZYBKI metabolizer - kawa bezpieczna"
            ;;
        "0/1")
            echo "rs762551: $LETTERS → ⚠️ ŚREDNI metabolizer - ogranicz do 2-3 kaw"
            ;;
        "1/1")
            echo "rs762551: $LETTERS → 🐢 WOLNY metabolizer - ogranicz kawę!"
            ;;
    esac
else
    echo "Brak danych"
fi
echo ""

#===============================================================================
# METABOLIZM ALKOHOLU
#===============================================================================

echo "🍷 METABOLIZM ALKOHOLU"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Asian flush - ALDH2 rs671: 12:111803962 G>A
GT=$(gt_at 12:111803962 G A rs671)

if [ -n "$GT" ]; then
    if [ "$GT" != "0/0" ]; then
        echo "rs671 ALDH2: wariant → 🔴 Asian flush! Zwiększone ryzyko raka przy alkoholu"
    else
        echo "rs671 ALDH2: G/G → ✅ Normalny metabolizm aldehydu octowego"
    fi
else
    echo "rs671: Brak danych (prawdopodobnie OK jeśli europejskie pochodzenie)"
fi
echo ""

#===============================================================================
# APOE - ALZHEIMER / SERCE
#===============================================================================

echo "🧠 APOE (Alzheimer / Serce)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs429358 19:44908684 T>C (C = ε4), rs7412 19:44908822 C>T (T = ε2)
APOE1=$(gt_at 19:44908684 T C rs429358)
APOE2=$(gt_at 19:44908822 C T rs7412)

if [ -n "$APOE1" ]; then
    echo "rs429358 (ε4): $APOE1 ($(gt_letters "$APOE1" T C))"
else
    echo "rs429358 (ε4): brak danych"
fi
if [ -n "$APOE2" ]; then
    echo "rs7412 (ε2): $APOE2 ($(gt_letters "$APOE2" C T))"
else
    echo "rs7412 (ε2): brak danych"
fi

# Interpretacja (uproszczona)
if [[ "$APOE1" == "0/0" ]] && [[ "$APOE2" == "0/0" ]]; then
    echo "→ Prawdopodobnie ε3/ε3 ✅ OPTYMALNY"
elif [[ "$APOE1" == *"1"* ]]; then
    echo "→ ⚠️ Obecny allel ε4 - podwyższone ryzyko Alzheimera"
fi
echo ""

#===============================================================================
# MTHFR
#===============================================================================

echo "🧬 MTHFR (Metabolizm folianów)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs1801133 (C677T, nić minus): 1:11796321 G>A; A = 677T
GT=$(gt_at 1:11796321 G A rs1801133)

if [ -n "$GT" ]; then
    case $GT in
        "0/0")
            echo "rs1801133 C677T: C/C → ✅ Normalny enzym"
            ;;
        "0/1")
            echo "rs1801133 C677T: C/T → ⚠️ 35% redukcja - rozważ metylofolian"
            ;;
        "1/1")
            echo "rs1801133 C677T: T/T → 🔴 70% redukcja - POTRZEBNY metylofolian!"
            ;;
    esac
else
    echo "Brak danych"
fi
echo ""

#===============================================================================
# ŁYSIENIE (dla mężczyzn)
#===============================================================================

if [ "$SEX" == "M" ]; then
    echo "👨‍🦲 ŁYSIENIE ANDROGENOWE"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    BALD=$(hx_gt 20:21985252)  # TODO(build): unverified GRCh38

    if [ -n "$BALD" ]; then
        if [ "$BALD" != "0/0" ]; then
            echo "rs2180439: wariant → 📈 Zwiększone ryzyko łysienia"
        else
            echo "rs2180439: brak wariantu → Niższe ryzyko łysienia"
        fi
    fi
    echo ""
fi

#===============================================================================
# SMAK GORZKI (brokuły, kawa)
#===============================================================================

echo "🥦 PERCEPCJA SMAKU GORZKIEGO"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs713598 A49P (7:141972804 GRCh38, C>G): TAS2R38 leży na nici minus, więc
# plus-strand C = Ala49 (AVI, non-taster), G = Pro49 (PAV, taster); genotyp
# liczony jako dawka G = liczba alleli PAV.
# TODO(build): unverified GRCh38 orientation (table row is 'unverified')
BITTER=$(gt_at 7:141972804 C G rs713598)

if [ -n "$BITTER" ]; then
    case $BITTER in
        "1/1")
            echo "TAS2R38: PAV/PAV → 😖 Super-taster! Silnie czujesz gorzki smak"
            ;;
        "0/1")
            echo "TAS2R38: PAV/AVI → 😐 Średnio czujesz gorzki smak"
            ;;
        "0/0")
            echo "TAS2R38: AVI/AVI → 😋 Non-taster - nie czujesz gorzkiego!"
            ;;
    esac
fi
echo ""

#===============================================================================
# WOSKOWINA USZNA
#===============================================================================

echo "👂 TYP WOSKOWINY USZNEJ"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs17822931 (ABCC11): 16:48224287 C>T; C/C i C/T = mokra (dominująca), T/T = sucha
EW_REF=C; EW_ALT=T
GT=$(gt_at 16:48224287 "$EW_REF" "$EW_ALT" rs17822931)

if [ -n "$GT" ]; then
    LETTERS=$(gt_letters "$GT" "$EW_REF" "$EW_ALT")
    case $GT in
        "0/0")
            echo "rs17822931: $LETTERS → 💧 Mokra woskowina (typ europejski)"
            ;;
        "0/1")
            echo "rs17822931: $LETTERS → 💧 Mokra woskowina (allel mokry jest dominujący)"
            ;;
        "1/1")
            echo "rs17822931: $LETTERS → 🔸 Sucha woskowina (typ azjatycki)"
            ;;
    esac
fi
echo ""

#===============================================================================
# KICHANIE NA SŁOŃCE
#===============================================================================

echo "☀️ ODRUCH KICHANIA NA SŁOŃCE (photic sneeze)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Pozycja GRCh38 rs10427255 (~2:145.37 Mb, okolice ZEB2) nie została
# zweryfikowana; poprzednio odczytywano tu genotyp LCT rs4988235 (2:135851076).
echo "rs10427255: pozycja niezweryfikowana - pominięto"
echo ""

#===============================================================================
# SIŁA MIĘŚNI
#===============================================================================

echo "💪 GENETYKA SPORTOWA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# rs1815739 (ACTN3 R577X): 11:66560624 C>T
GT=$(gt_at 11:66560624 C T rs1815739)

if [ -n "$GT" ]; then
    case $GT in
        "0/0")
            echo "ACTN3 rs1815739: C/C (R/R) → 🏋️ ELITARNA genetyka siłowa!"
            ;;
        "0/1")
            echo "ACTN3 rs1815739: C/T (R/X) → 🏃 Dobra mieszanka siła/wytrzymałość"
            ;;
        "1/1")
            echo "ACTN3 rs1815739: T/T (X/X) → 🚴 Genetyka wytrzymałościowa"
            ;;
    esac
fi
echo ""

#===============================================================================
# PODSUMOWANIE
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║                    ANALIZA ZAKOŃCZONA                             ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "Więcej analiz dostępnych przez: ./helixight_local_toolkit.sh $VCF"
echo ""
hx_report_coverage
