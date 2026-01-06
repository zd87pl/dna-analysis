#!/bin/bash
#===============================================================================
# POLYGENIC RISK SCORES (PRS) - FIXED VERSION
# Wielogenowe wyniki ryzyka dla zlozonych chorob
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║           POLYGENIC RISK SCORES (PRS)                                     ║"
echo "║           Wielogenowe ryzyko zlozonych chorob                             ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "PRS laczy wiele wariantow genetycznych, z ktorych kazdy ma maly efekt,"
echo "w jeden wynik ryzyka. To bardziej zaawansowane niz pojedyncze SNP!"
echo ""

#===============================================================================
# PRS: CHOROBA WIENCOWA (CAD)
#===============================================================================

echo "❤️ CHOROBA WIENCOWA (CAD) - Polygenic Risk Score"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

CAD_SCORE=0
CAD_MAX=20

echo "Sprawdzam markery ryzyka CAD..."
echo ""

# Format: pozycja|rsid|OR
CAD_SNPS="
9:22125503|rs1333049|1.30
1:109818530|rs17465637|1.15
10:44765623|rs1746048|1.12
6:160961137|rs6922269|1.13
2:203454130|rs2306374|1.11
21:35593836|rs9982601|1.15
6:12903957|rs12526453|1.10
1:56962821|rs11206510|1.13
19:11202306|rs1122608|1.15
6:161010118|rs3798220|1.42
"

while IFS='|' read -r pos rsid odds_ratio; do
    [ -z "$pos" ] && continue
    
    result=$(bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
    
    if [ -n "$result" ]; then
        case $result in
            "0/1"|"0|1"|"1/0"|"1|0")
                CAD_SCORE=$((CAD_SCORE + 1))
                echo "  $rsid: Heterozygota (+1 pkt, OR=$odds_ratio)"
                ;;
            "1/1"|"1|1")
                CAD_SCORE=$((CAD_SCORE + 2))
                echo "  $rsid: Homozygota ryzyka (+2 pkt, OR=$odds_ratio)"
                ;;
        esac
    fi
done <<< "$CAD_SNPS"

echo ""
CAD_PERCENTILE=$((CAD_SCORE * 100 / CAD_MAX))
echo "Wynik CAD PRS: $CAD_SCORE / $CAD_MAX punktow (~$CAD_PERCENTILE percentyl)"

if [ $CAD_SCORE -le 4 ]; then
    echo "Interpretacja: ✅ NISKIE ryzyko genetyczne CAD"
elif [ $CAD_SCORE -le 10 ]; then
    echo "Interpretacja: 🟡 SREDNIE ryzyko genetyczne CAD"
else
    echo "Interpretacja: 🔴 PODWYZSZONE ryzyko genetyczne CAD"
    echo "→ Rekomendacja: Regularne badania lipidow, zdrowa dieta, aktywnosc fizyczna"
fi
echo ""

#===============================================================================
# PRS: CUKRZYCA TYPU 2 (T2D)
#===============================================================================

echo "🩸 CUKRZYCA TYPU 2 (T2D) - Polygenic Risk Score"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

T2D_SCORE=0
T2D_MAX=16

echo "Sprawdzam markery ryzyka T2D..."
echo ""

T2D_SNPS="
10:114758349|rs7903146|1.37
3:185511687|rs4402960|1.14
6:20679709|rs1801282|1.14
9:22134094|rs10811661|1.20
10:12307894|rs12779790|1.10
11:92708710|rs5215|1.11
16:53820527|rs9939609|1.15
8:118184783|rs13266634|1.12
"

while IFS='|' read -r pos rsid odds_ratio; do
    [ -z "$pos" ] && continue
    
    result=$(bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
    
    if [ -n "$result" ]; then
        case $result in
            "0/1"|"0|1"|"1/0"|"1|0")
                T2D_SCORE=$((T2D_SCORE + 1))
                echo "  $rsid: Heterozygota (+1 pkt, OR=$odds_ratio)"
                ;;
            "1/1"|"1|1")
                T2D_SCORE=$((T2D_SCORE + 2))
                echo "  $rsid: Homozygota ryzyka (+2 pkt, OR=$odds_ratio)"
                ;;
        esac
    fi
done <<< "$T2D_SNPS"

echo ""
T2D_PERCENTILE=$((T2D_SCORE * 100 / T2D_MAX))
echo "Wynik T2D PRS: $T2D_SCORE / $T2D_MAX punktow (~$T2D_PERCENTILE percentyl)"

if [ $T2D_SCORE -le 4 ]; then
    echo "Interpretacja: ✅ NISKIE ryzyko genetyczne T2D"
elif [ $T2D_SCORE -le 8 ]; then
    echo "Interpretacja: 🟡 SREDNIE ryzyko genetyczne T2D"
else
    echo "Interpretacja: 🔴 PODWYZSZONE ryzyko genetyczne T2D"
    echo "→ Rekomendacja: Dieta niskoglikemiczna, kontrola masy ciala, regularna glukoza"
fi
echo ""

#===============================================================================
# PRS: RAK PROSTATY (tylko mezczyzni)
#===============================================================================

echo "🔵 RAK PROSTATY - Polygenic Risk Score"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

PROSTATE_SCORE=0
PROSTATE_MAX=14

echo "Sprawdzam markery ryzyka raka prostaty..."
echo ""

PROSTATE_SNPS="
8:128355618|rs16901979|1.79
8:128093297|rs6983267|1.26
17:69108753|rs1859962|1.20
10:51549496|rs10993994|1.25
7:97614628|rs6465657|1.12
2:43407453|rs2660753|1.15
22:43500212|rs9623117|1.18
"

while IFS='|' read -r pos rsid odds_ratio; do
    [ -z "$pos" ] && continue
    
    result=$(bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
    
    if [ -n "$result" ]; then
        case $result in
            "0/1"|"0|1"|"1/0"|"1|0")
                PROSTATE_SCORE=$((PROSTATE_SCORE + 1))
                echo "  $rsid: Heterozygota (+1 pkt, OR=$odds_ratio)"
                ;;
            "1/1"|"1|1")
                PROSTATE_SCORE=$((PROSTATE_SCORE + 2))
                echo "  $rsid: Homozygota ryzyka (+2 pkt, OR=$odds_ratio)"
                ;;
        esac
    fi
done <<< "$PROSTATE_SNPS"

echo ""
PROSTATE_PERCENTILE=$((PROSTATE_SCORE * 100 / PROSTATE_MAX))
echo "Wynik Prostate PRS: $PROSTATE_SCORE / $PROSTATE_MAX punktow"

if [ $PROSTATE_SCORE -le 3 ]; then
    echo "Interpretacja: ✅ NISKIE ryzyko genetyczne"
elif [ $PROSTATE_SCORE -le 7 ]; then
    echo "Interpretacja: 🟡 SREDNIE ryzyko genetyczne"
else
    echo "Interpretacja: 🔴 PODWYZSZONE ryzyko genetyczne"
    echo "→ Rekomendacja: PSA i badania od 45 r.z. (zamiast 50)"
fi
echo ""

#===============================================================================
# PRS: CHOROBA ALZHEIMERA
#===============================================================================

echo "🧠 CHOROBA ALZHEIMERA - Polygenic Risk Score"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

ALZ_SCORE=0
ALZ_MAX=14

# APOE jest dominujacy - warto go osobno
APOE1=$(bcftools query -r 19:45411941 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
APOE2=$(bcftools query -r 19:45412079 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)

if [[ "$APOE1" == *"1"* ]]; then
    # APOE e4 obecny
    if [[ "$APOE1" == "1/1" ]] || [[ "$APOE1" == "1|1" ]]; then
        ALZ_SCORE=$((ALZ_SCORE + 6))
        echo "  APOE e4/e4: WYSOKIE ryzyko (+6 pkt, OR~12)"
    else
        ALZ_SCORE=$((ALZ_SCORE + 3))
        echo "  APOE e4 heterozygota: Podwyzszone ryzyko (+3 pkt, OR~3)"
    fi
elif [[ "$APOE2" == *"1"* ]]; then
    echo "  APOE e2 obecny: OCHRONNY wariant (-1 pkt)"
    ALZ_SCORE=$((ALZ_SCORE - 1))
else
    echo "  APOE e3/e3: Standardowe ryzyko (0 pkt)"
fi

# Inne SNP
ALZ_SNPS="
8:27219987|rs744373|1.15
2:127892810|rs6733839|1.20
11:121435587|rs11218343|0.90
6:47452463|rs9271192|1.15
"

while IFS='|' read -r pos rsid odds_ratio; do
    [ -z "$pos" ] && continue
    
    result=$(bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
    
    if [ -n "$result" ] && [ "$result" != "0/0" ] && [ "$result" != "0|0" ]; then
        # Sprawdz czy ochronny (OR < 1)
        is_protective=$(echo "$odds_ratio < 1" | bc -l 2>/dev/null || echo "0")
        if [ "$is_protective" == "1" ]; then
            echo "  $rsid: Wariant OCHRONNY (-1 pkt)"
            ALZ_SCORE=$((ALZ_SCORE - 1))
        else
            echo "  $rsid: Wariant ryzyka (+1 pkt)"
            ALZ_SCORE=$((ALZ_SCORE + 1))
        fi
    fi
done <<< "$ALZ_SNPS"

echo ""
echo "Wynik Alzheimer PRS: $ALZ_SCORE / $ALZ_MAX punktow"

if [ $ALZ_SCORE -le 0 ]; then
    echo "Interpretacja: ✅ NISKIE ryzyko genetyczne (korzystne warianty)"
elif [ $ALZ_SCORE -le 3 ]; then
    echo "Interpretacja: 🟢 STANDARDOWE ryzyko genetyczne"
elif [ $ALZ_SCORE -le 6 ]; then
    echo "Interpretacja: 🟡 UMIARKOWANIE PODWYZSZONE ryzyko"
else
    echo "Interpretacja: 🔴 ZNACZNIE PODWYZSZONE ryzyko genetyczne"
    echo "→ Rekomendacja: Aktywnosc fizyczna, dieta srodziemnomorska, aktywnosc umyslowa"
fi
echo ""

#===============================================================================
# PODSUMOWANIE PRS
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    PODSUMOWANIE POLYGENIC RISK SCORES                     ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Tabela
printf "%-25s %-15s %-20s\n" "Choroba" "Wynik" "Interpretacja"
echo "─────────────────────────────────────────────────────────────────────────────"

printf "%-25s %-15s " "Choroba wiencowa (CAD)" "$CAD_SCORE/$CAD_MAX"
if [ $CAD_SCORE -le 4 ]; then echo "✅ Niskie"; elif [ $CAD_SCORE -le 10 ]; then echo "🟡 Srednie"; else echo "🔴 Wysokie"; fi

printf "%-25s %-15s " "Cukrzyca typu 2" "$T2D_SCORE/$T2D_MAX"
if [ $T2D_SCORE -le 4 ]; then echo "✅ Niskie"; elif [ $T2D_SCORE -le 8 ]; then echo "🟡 Srednie"; else echo "🔴 Wysokie"; fi

printf "%-25s %-15s " "Rak prostaty" "$PROSTATE_SCORE/$PROSTATE_MAX"
if [ $PROSTATE_SCORE -le 3 ]; then echo "✅ Niskie"; elif [ $PROSTATE_SCORE -le 7 ]; then echo "🟡 Srednie"; else echo "🔴 Wysokie"; fi

printf "%-25s %-15s " "Choroba Alzheimera" "$ALZ_SCORE/$ALZ_MAX"
if [ $ALZ_SCORE -le 0 ]; then echo "✅ Niskie"; elif [ $ALZ_SCORE -le 3 ]; then echo "🟢 Standardowe"; elif [ $ALZ_SCORE -le 6 ]; then echo "🟡 Podwyzszone"; else echo "🔴 Wysokie"; fi

echo ""
echo "UWAGA: To uproszczone PRS (10-20 SNP). Kliniczne PRS uzywaja tysiecy SNP"
echo "i sa bardziej dokladne. Wyniki maja charakter informacyjny."
echo ""
