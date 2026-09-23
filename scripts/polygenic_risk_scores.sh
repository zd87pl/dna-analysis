#!/bin/bash
#===============================================================================
# POLYGENIC RISK SCORES (PRS)
# Wielogenowe wyniki ryzyka dla zlozonych chorob
#
# Wspolrzedne: GRCh38, nazwy kontigow bez prefiksu "chr" (biblioteka
# scripts/lib/vcf_common.sh normalizuje plik wejsciowy przed zapytaniami).
# Kanoniczna tabela: scripts/lib/GRCh38_COORDINATES.md
#===============================================================================

# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"
set -eo pipefail

VCF="${1:-}"
hx_require_arg "$VCF" "$(basename "$0") <variants.vcf.gz>"

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║           POLYGENIC RISK SCORES (PRS)                                     ║"
echo "║           Wielogenowe ryzyko zlozonych chorob                             ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "PRS laczy wiele wariantow genetycznych, z ktorych kazdy ma maly efekt,"
echo "w jeden wynik ryzyka. To bardziej zaawansowane niz pojedyncze SNP!"
echo ""

hx_prepare_vcf "$VCF"

#===============================================================================
# POMOCNICZE: liczba kopii allelu ryzyka
#===============================================================================

# prs_dosage POS RSID RISK_ALLELE [REF ALT]
# Wypisuje liczbe kopii (0/1/2) allelu ryzyka w rekordzie na pozycji POS,
# "NA" gdy nie ma rekordu (lub rekord ma inne REF/ALT niz tabela - wtedy
# biblioteka wypisuje ostrzezenie na stderr) albo "MISMATCH" gdy allel ryzyka
# nie jest ani REF, ani ALT rekordu (inna nic / build / wspolrzedna).
# Wiersze tabeli bez REF/ALT czytaja allele z rekordu w pliku VCF.
prs_dosage() {
    local pos=$1 rsid=$2 risk ref=${4:-} alt=${5:-}
    risk=$(printf '%s' "$3" | tr '[:lower:]' '[:upper:]')
    local start rec p r a out dosage
    if [[ -z $ref || -z $alt ]]; then
        start=${pos#*:}
        ref=""
        alt=""
        rec=$(hx_query "$pos" '%POS\t%REF\t%ALT\n')
        while IFS=$'\t' read -r p r a; do
            [[ $p == "$start" ]] || continue
            ref=$r
            alt=$a
            break
        done <<<"$rec"
        if [[ -z $ref ]]; then
            # brak rekordu: placeholder, zeby hx_gt_allele policzylo zapytanie
            ref=N
            alt=N
        elif [[ ",$(printf '%s' "$alt" | tr '[:lower:]' '[:upper:]')," == *",$risk,"* ]]; then
            alt=$risk
        else
            alt=${alt%%,*}
        fi
    fi
    out=$(hx_gt_allele "$pos" "$ref" "$alt" "$rsid")
    if [[ -z $out ]]; then
        printf 'NA\n'
        return 0
    fi
    dosage=${out%%$'\t'*}
    if [[ $(printf '%s' "$alt" | tr '[:lower:]' '[:upper:]') == "$risk" ]]; then
        printf '%s\n' "$dosage"
    elif [[ $(printf '%s' "$ref" | tr '[:lower:]' '[:upper:]') == "$risk" ]]; then
        printf '%s\n' "$((2 - dosage))"
    else
        printf 'MISMATCH\n'
    fi
}

#===============================================================================
# PRS: CHOROBA WIENCOWA (CAD)
#===============================================================================

echo "❤️ CHOROBA WIENCOWA (CAD) - Polygenic Risk Score"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

CAD_SCORE=0
CAD_MAX=20

# Kluczowe SNP dla CAD (uproszczone - prawdziwe PRS ma tysiace SNP)
# Format: rsid|allel_ryzyka|OR|REF|ALT (REF/ALT puste = czytane z rekordu VCF)
declare -A CAD_SNPS=(
    ["9:22125504"]="rs1333049|C|1.3|G|C"     # 9p21 - najsilniejszy!  # TODO(build): unverified GRCh38
    ["1:109818530"]="rs17465637|C|1.15"      # MIA3   # TODO(build): unverified GRCh38
    ["10:44765623"]="rs1746048|C|1.12"       # CXCL12 # TODO(build): unverified GRCh38
    ["6:160961137"]="rs6922269|A|1.13"       # MTHFD1L # TODO(build): unverified GRCh38
    ["2:203454130"]="rs2306374|C|1.11"       # WDR12  # TODO(build): unverified GRCh38
    ["21:35593836"]="rs9982601|T|1.15"       # MRPS6  # TODO(build): unverified GRCh38
    ["6:12903957"]="rs12526453|C|1.10"       # PHACTR1 # TODO(build): unverified GRCh38
    ["1:56962821"]="rs11206510|T|1.13"       # PCSK9  # TODO(build): unverified GRCh38
    ["19:11202306"]="rs1122608|G|1.15"       # LDLR   # TODO(build): unverified GRCh38
    ["6:161010118"]="rs3798220|C|1.42"       # LPA    # TODO(build): unverified GRCh38
)

echo "Sprawdzam markery ryzyka CAD..."
echo ""

for pos in "${!CAD_SNPS[@]}"; do
    IFS='|' read -r rsid risk_allele odds_ratio ref alt <<< "${CAD_SNPS[$pos]}"

    dosage=$(prs_dosage "$pos" "$rsid" "$risk_allele" "$ref" "$alt")

    # Prosta punktacja: 0, 1 lub 2 kopie allelu ryzyka
    case $dosage in
        NA)
            echo "  $rsid: nie genotypowano (0 pkt)"
            ;;
        MISMATCH)
            echo "  $rsid: allel ryzyka $risk_allele nie pasuje do REF/ALT - pominieto"
            ;;
        0)
            echo "  $rsid: brak allelu ryzyka (0 pkt, OR=$odds_ratio)"
            ;;
        1)
            CAD_SCORE=$((CAD_SCORE + 1))
            echo "  $rsid: Heterozygota (+1 pkt, OR=$odds_ratio)"
            ;;
        2)
            CAD_SCORE=$((CAD_SCORE + 2))
            echo "  $rsid: Homozygota ryzyka (+2 pkt, OR=$odds_ratio)"
            ;;
    esac
done

echo ""
CAD_PCT=$((CAD_SCORE * 100 / CAD_MAX))
echo "Wynik CAD PRS: $CAD_SCORE / $CAD_MAX punktow (${CAD_PCT}% maksimum; nie jest to percentyl populacyjny)"

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

# Format: rsid|allel_ryzyka|OR|REF|ALT (REF/ALT puste = czytane z rekordu VCF)
declare -A T2D_SNPS=(
    ["10:112998590"]="rs7903146|T|1.37|C|T"  # TCF7L2 - najsilniejszy!
    ["3:185511687"]="rs4402960|T|1.14"       # IGF2BP2 # TODO(build): unverified GRCh38
    ["3:12351626"]="rs1801282|C|1.14|C|G"    # PPARG (allel ryzyka C = REF)
    ["9:22134095"]="rs10811661|T|1.20|T|C"   # CDKN2A/B (allel ryzyka T = REF) # TODO(build): unverified GRCh38
    ["10:12307894"]="rs12779790|G|1.10"      # CDC123 # TODO(build): unverified GRCh38
    ["11:92708710"]="rs5215|C|1.11"          # KCNJ11 # TODO(build): unverified GRCh38
    ["16:53786615"]="rs9939609|A|1.15|T|A"   # FTO
    ["8:117172544"]="rs13266634|C|1.12|C|T"  # SLC30A8 (allel ryzyka C = REF) # TODO(build): unverified GRCh38
)

echo "Sprawdzam markery ryzyka T2D..."
echo ""

for pos in "${!T2D_SNPS[@]}"; do
    IFS='|' read -r rsid risk_allele odds_ratio ref alt <<< "${T2D_SNPS[$pos]}"

    dosage=$(prs_dosage "$pos" "$rsid" "$risk_allele" "$ref" "$alt")

    case $dosage in
        NA)
            echo "  $rsid: nie genotypowano (0 pkt)"
            ;;
        MISMATCH)
            echo "  $rsid: allel ryzyka $risk_allele nie pasuje do REF/ALT - pominieto"
            ;;
        0)
            echo "  $rsid: brak allelu ryzyka (0 pkt)"
            ;;
        1)
            T2D_SCORE=$((T2D_SCORE + 1))
            echo "  $rsid: Heterozygota (+1 pkt)"
            ;;
        2)
            T2D_SCORE=$((T2D_SCORE + 2))
            echo "  $rsid: Homozygota ryzyka (+2 pkt)"
            ;;
    esac
done

echo ""
T2D_PCT=$((T2D_SCORE * 100 / T2D_MAX))
echo "Wynik T2D PRS: $T2D_SCORE / $T2D_MAX punktow (${T2D_PCT}% maksimum; nie jest to percentyl populacyjny)"

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

# Format: rsid|allel_ryzyka|OR|REF|ALT (REF/ALT puste = czytane z rekordu VCF)
declare -A PROSTATE_SNPS=(
    ["8:128355618"]="rs16901979|A|1.79"   # 8q24 region # TODO(build): unverified GRCh38
    ["8:128093297"]="rs6983267|G|1.26"    # 8q24        # TODO(build): unverified GRCh38
    ["17:69108753"]="rs1859962|G|1.20"    # HNF1B       # TODO(build): unverified GRCh38
    ["10:51549496"]="rs10993994|T|1.25"   # MSMB        # TODO(build): unverified GRCh38
    ["7:97614628"]="rs6465657|C|1.12"     # LMTK2       # TODO(build): unverified GRCh38
    ["2:43407453"]="rs2660753|T|1.15"     # intergenic  # TODO(build): unverified GRCh38
    ["22:43500212"]="rs9623117|C|1.18"    # BIK         # TODO(build): unverified GRCh38
)

echo "Sprawdzam markery ryzyka raka prostaty..."
echo ""

for pos in "${!PROSTATE_SNPS[@]}"; do
    IFS='|' read -r rsid risk_allele odds_ratio ref alt <<< "${PROSTATE_SNPS[$pos]}"

    dosage=$(prs_dosage "$pos" "$rsid" "$risk_allele" "$ref" "$alt")

    case $dosage in
        NA)
            echo "  $rsid: nie genotypowano (0 pkt)"
            ;;
        MISMATCH)
            echo "  $rsid: allel ryzyka $risk_allele nie pasuje do REF/ALT - pominieto"
            ;;
        0)
            echo "  $rsid: brak allelu ryzyka (0 pkt)"
            ;;
        1)
            PROSTATE_SCORE=$((PROSTATE_SCORE + 1))
            echo "  $rsid: Heterozygota (+1 pkt)"
            ;;
        2)
            PROSTATE_SCORE=$((PROSTATE_SCORE + 2))
            echo "  $rsid: Homozygota ryzyka (+2 pkt)"
            ;;
    esac
done

echo ""
PROSTATE_PCT=$((PROSTATE_SCORE * 100 / PROSTATE_MAX))
echo "Wynik Prostate PRS: $PROSTATE_SCORE / $PROSTATE_MAX punktow (${PROSTATE_PCT}% maksimum; nie jest to percentyl populacyjny)"

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
# rs429358 19:44908684 T>C (C = allel e4), rs7412 19:44908822 C>T (T = allel e2)
APOE_E4=$(prs_dosage 19:44908684 rs429358 C T C)
APOE_E2=$(prs_dosage 19:44908822 rs7412 T C T)

case $APOE_E4 in
    2)
        ALZ_SCORE=$((ALZ_SCORE + 6))
        echo "  APOE e4/e4: WYSOKIE ryzyko (+6 pkt, OR~12)"
        ;;
    1)
        ALZ_SCORE=$((ALZ_SCORE + 3))
        echo "  APOE e4 heterozygota: Podwyzszone ryzyko (+3 pkt, OR~3)"
        ;;
    *)
        if [[ $APOE_E2 == 1 || $APOE_E2 == 2 ]]; then
            echo "  APOE e2 obecny: OCHRONNY wariant (-1 pkt)"
            ALZ_SCORE=$((ALZ_SCORE - 1))
        elif [[ $APOE_E4 == 0 || $APOE_E2 == 0 ]]; then
            echo "  APOE e3/e3: Standardowe ryzyko (0 pkt)"
        else
            echo "  APOE (rs429358/rs7412): nie genotypowano (0 pkt)"
        fi
        ;;
esac

# Inne SNP
# Format: rsid|allel_efektu|OR|REF|ALT (OR < 1 = wariant ochronny)
declare -A ALZ_SNPS=(
    ["8:27219987"]="rs744373|T|1.15"      # BIN1  # TODO(build): unverified GRCh38
    ["2:127892810"]="rs6733839|T|1.20"    # BIN1  # TODO(build): unverified GRCh38
    ["11:121435587"]="rs11218343|C|0.90"  # SORL1 (ochronny) # TODO(build): unverified GRCh38
    ["6:47452463"]="rs9271192|C|1.15"     # HLA-DRB5 # TODO(build): unverified GRCh38
)

for pos in "${!ALZ_SNPS[@]}"; do
    IFS='|' read -r rsid risk_allele odds_ratio ref alt <<< "${ALZ_SNPS[$pos]}"

    dosage=$(prs_dosage "$pos" "$rsid" "$risk_allele" "$ref" "$alt")

    case $dosage in
        NA)
            echo "  $rsid: nie genotypowano (0 pkt)"
            ;;
        MISMATCH)
            echo "  $rsid: allel ryzyka $risk_allele nie pasuje do REF/ALT - pominieto"
            ;;
        0)
            echo "  $rsid: brak allelu efektu (0 pkt)"
            ;;
        *)
            case $odds_ratio in
                0.*)
                    echo "  $rsid: Wariant OCHRONNY (-$dosage)"
                    ALZ_SCORE=$((ALZ_SCORE - dosage))
                    ;;
                *)
                    echo "  $rsid: Wariant ryzyka (+$dosage)"
                    ALZ_SCORE=$((ALZ_SCORE + dosage))
                    ;;
            esac
            ;;
    esac
done

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
[ $CAD_SCORE -le 4 ] && echo "✅ Niskie" || ([ $CAD_SCORE -le 10 ] && echo "🟡 Srednie" || echo "🔴 Wysokie")

printf "%-25s %-15s " "Cukrzyca typu 2" "$T2D_SCORE/$T2D_MAX"
[ $T2D_SCORE -le 4 ] && echo "✅ Niskie" || ([ $T2D_SCORE -le 8 ] && echo "🟡 Srednie" || echo "🔴 Wysokie")

printf "%-25s %-15s " "Rak prostaty" "$PROSTATE_SCORE/$PROSTATE_MAX"
[ $PROSTATE_SCORE -le 3 ] && echo "✅ Niskie" || ([ $PROSTATE_SCORE -le 7 ] && echo "🟡 Srednie" || echo "🔴 Wysokie")

printf "%-25s %-15s " "Choroba Alzheimera" "$ALZ_SCORE/$ALZ_MAX"
[ $ALZ_SCORE -le 0 ] && echo "✅ Niskie" || ([ $ALZ_SCORE -le 3 ] && echo "🟢 Standardowe" || ([ $ALZ_SCORE -le 6 ] && echo "🟡 Podwyzszone" || echo "🔴 Wysokie"))

echo ""
echo "UWAGA: To uproszczone PRS (10-20 SNP). Kliniczne PRS uzywaja tysiecy SNP"
echo "i sa bardziej dokladne. Wyniki maja charakter informacyjny."
echo ""
hx_report_coverage
