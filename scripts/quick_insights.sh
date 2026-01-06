#!/bin/bash
#===============================================================================
# QUICK GENOME INSIGHTS
# Natychmiastowe interesujące odkrycia z VCF
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║              QUICK GENOME INSIGHTS                                ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

if [ ! -f "$VCF" ]; then
    echo "Błąd: Nie znaleziono pliku $VCF"
    exit 1
fi

#===============================================================================
# PODSTAWOWE STATYSTYKI
#===============================================================================

echo "📊 PODSTAWOWE STATYSTYKI"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

TOTAL=$(bcftools view -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
SNPS=$(bcftools view -v snps -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
INDELS=$(bcftools view -v indels -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')

printf "Wszystkie warianty: %'d\n" $TOTAL
printf "SNPs:               %'d\n" $SNPS
printf "Indele:             %'d\n" $INDELS
echo ""

#===============================================================================
# PŁEĆ GENETYCZNA
#===============================================================================

echo "👤 PŁEĆ GENETYCZNA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

X_COUNT=$(bcftools view -r X -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
Y_COUNT=$(bcftools view -r Y -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')

echo "Warianty chr X: $X_COUNT"
echo "Warianty chr Y: $Y_COUNT"

if [ "$Y_COUNT" -gt "500" ]; then
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

# rs12913832 - główny marker koloru oczu
EYE_MAIN=$(bcftools query -r 15:28365618 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$EYE_MAIN" ]; then
    GT=$(echo $EYE_MAIN | awk '{print $3}')
    case $GT in
        "0/0"|"0|0")
            echo "rs12913832: A/A → 🔵 Prawdopodobnie NIEBIESKIE oczy"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "rs12913832: A/G → 🟢 Prawdopodobnie ZIELONE lub PIWNE oczy"
            ;;
        "1/1"|"1|1")
            echo "rs12913832: G/G → 🟤 Prawdopodobnie BRĄZOWE oczy"
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

# MC1R - rude włosy
MC1R_COUNT=0
for pos in 89985844 89985940 89986091 89986117 89986130 89986144; do
    result=$(bcftools view -r "16:$pos" -H "$VCF" 2>/dev/null | wc -l)
    MC1R_COUNT=$((MC1R_COUNT + result))
done

if [ "$MC1R_COUNT" -gt "0" ]; then
    echo "Warianty MC1R: $MC1R_COUNT → 🔴 Możliwe RUDE włosy lub piegi"
else
    echo "Warianty MC1R: 0 → Brak predyspozycji do rudości"
fi

# KITLG - blond
BLOND=$(bcftools query -r 12:89328335 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ "$BLOND" == "0/1" ] || [ "$BLOND" == "1/1" ]; then
    echo "rs12821256: wariant → 👱 Skłonność do JASNYCH włosów"
fi
echo ""

#===============================================================================
# TOLERANCJA LAKTOZY
#===============================================================================

echo "🥛 TOLERANCJA LAKTOZY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

LCT=$(bcftools query -r 2:136608646 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$LCT" ]; then
    GT=$(echo $LCT | awk '{print $3}')
    case $GT in
        "0/0"|"0|0")
            echo "rs4988235: C/C → ❌ NIETOLERANCJA laktozy (typ afrykański/azjatycki)"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "rs4988235: C/T → ⚠️ CZĘŚCIOWA tolerancja"
            ;;
        "1/1"|"1|1")
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

CAFFEINE=$(bcftools query -r 15:75041917 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$CAFFEINE" ]; then
    GT=$(echo $CAFFEINE | awk '{print $3}')
    case $GT in
        "0/0"|"0|0")
            echo "rs762551: A/A → ⚡ SZYBKI metabolizer - kawa bezpieczna"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "rs762551: A/C → ⚠️ ŚREDNI metabolizer - ogranicz do 2-3 kaw"
            ;;
        "1/1"|"1|1")
            echo "rs762551: C/C → 🐢 WOLNY metabolizer - ogranicz kawę!"
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

# Asian flush - ALDH2
ALDH2=$(bcftools query -r 12:112241766 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$ALDH2" ]; then
    GT=$(echo $ALDH2 | awk '{print $3}')
    if [ "$GT" != "0/0" ] && [ "$GT" != "0|0" ]; then
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

APOE1=$(bcftools query -r 19:45411941 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)  # rs429358
APOE2=$(bcftools query -r 19:45412079 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)  # rs7412

echo "rs429358 (ε4): $APOE1"
echo "rs7412 (ε2): $APOE2"

# Interpretacja (uproszczona)
if [[ "$APOE1" == "0/0"* ]] && [[ "$APOE2" == "0/0"* ]]; then
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

MTHFR=$(bcftools query -r 1:11856378 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$MTHFR" ]; then
    GT=$(echo $MTHFR | awk '{print $3}')
    case $GT in
        "0/0"|"0|0")
            echo "rs1801133 C677T: C/C → ✅ Normalny enzym"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "rs1801133 C677T: C/T → ⚠️ 35% redukcja - rozważ metylofolian"
            ;;
        "1/1"|"1|1")
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
    
    BALD=$(bcftools query -r 20:21985252 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
    
    if [ -n "$BALD" ]; then
        if [ "$BALD" != "0/0" ] && [ "$BALD" != "0|0" ]; then
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

BITTER=$(bcftools query -r 7:141672604 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$BITTER" ]; then
    case $BITTER in
        "0/0"|"0|0")
            echo "TAS2R38: PAV/PAV → 😖 Super-taster! Silnie czujesz gorzki smak"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "TAS2R38: PAV/AVI → 😐 Średnio czujesz gorzki smak"
            ;;
        "1/1"|"1|1")
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

EARWAX=$(bcftools query -r 16:48258198 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$EARWAX" ]; then
    GT=$(echo $EARWAX | awk '{print $3}')
    case $GT in
        "0/0"|"0|0")
            echo "rs17822931: T/T → 💧 Mokra woskowina (typ europejski)"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "rs17822931: T/C → Mieszany typ"
            ;;
        "1/1"|"1|1")
            echo "rs17822931: C/C → 🔸 Sucha woskowina (typ azjatycki)"
            ;;
    esac
fi
echo ""

#===============================================================================
# KICHANIE NA SŁOŃCE
#===============================================================================

echo "☀️ ODRUCH KICHANIA NA SŁOŃCE (photic sneeze)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

SNEEZE=$(bcftools query -r 2:135851076 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$SNEEZE" ]; then
    if [ "$SNEEZE" != "0/0" ] && [ "$SNEEZE" != "0|0" ]; then
        echo "rs10427255: wariant → 🤧 TAK - kichasz patrząc na słońce!"
    else
        echo "rs10427255: brak → Nie kichasz na słońce"
    fi
fi
echo ""

#===============================================================================
# SIŁA MIĘŚNI
#===============================================================================

echo "💪 GENETYKA SPORTOWA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

ACTN3=$(bcftools query -r 11:66560624 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)

if [ -n "$ACTN3" ]; then
    GT=$(echo $ACTN3 | awk '{print $3}')
    case $GT in
        "0/0"|"0|0")
            echo "ACTN3 rs1815739: C/C (R/R) → 🏋️ ELITARNA genetyka siłowa!"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "ACTN3 rs1815739: C/T (R/X) → 🏃 Dobra mieszanka siła/wytrzymałość"
            ;;
        "1/1"|"1|1")
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
