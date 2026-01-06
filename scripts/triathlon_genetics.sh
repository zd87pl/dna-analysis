#!/bin/bash
#===============================================================================
#
#    ████████╗██████╗ ██╗ █████╗ ████████╗██╗  ██╗██╗      ██████╗ ███╗   ██╗
#    ╚══██╔══╝██╔══██╗██║██╔══██╗╚══██╔══╝██║  ██║██║     ██╔═══██╗████╗  ██║
#       ██║   ██████╔╝██║███████║   ██║   ███████║██║     ██║   ██║██╔██╗ ██║
#       ██║   ██╔══██╗██║██╔══██║   ██║   ██╔══██║██║     ██║   ██║██║╚██╗██║
#       ██║   ██║  ██║██║██║  ██║   ██║   ██║  ██║███████╗╚██████╔╝██║ ╚████║
#       ╚═╝   ╚═╝  ╚═╝╚═╝╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝
#
#     ██████╗ ███████╗███╗   ██╗███████╗████████╗██╗ ██████╗███████╗
#    ██╔════╝ ██╔════╝████╗  ██║██╔════╝╚══██╔══╝██║██╔════╝██╔════╝
#    ██║  ███╗█████╗  ██╔██╗ ██║█████╗     ██║   ██║██║     ███████╗
#    ██║   ██║██╔══╝  ██║╚██╗██║██╔══╝     ██║   ██║██║     ╚════██║
#    ╚██████╔╝███████╗██║ ╚████║███████╗   ██║   ██║╚██████╗███████║
#     ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝   ╚═╝   ╚═╝ ╚═════╝╚══════╝
#
#    KOMPLEKSOWA ANALIZA PREDYSPOZYCJI DO TRIATHLONU
#    🏊 Pływanie • 🚴 Kolarstwo • 🏃 Bieganie
#
#    Helixight Genetic Analysis Platform
#
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"
OUTDIR="triathlon_analysis"
mkdir -p "$OUTDIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                       ║"
echo "║     🏊 TRIATHLON GENETIC ANALYSIS 🚴 🏃                                               ║"
echo "║                                                                                       ║"
echo "║     Kompleksowa analiza predyspozycji genetycznych do triathlonu                      ║"
echo "║     Sprint • Olympic • Half-Ironman • Ironman                                         ║"
echo "║                                                                                       ║"
echo "╚═══════════════════════════════════════════════════════════════════════════════════════╝"
echo ""

#===============================================================================
# FUNKCJE
#===============================================================================

check_snp() {
    local pos="$1"
    bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1
}

# Wyniki dla każdej kategorii
VO2MAX_SCORE=0
VO2MAX_MAX=14

FAT_METABOLISM_SCORE=0
FAT_METABOLISM_MAX=10

LACTATE_SCORE=0
LACTATE_MAX=8

ENDURANCE_FIBER_SCORE=0
ENDURANCE_FIBER_MAX=10

INJURY_RESISTANCE_SCORE=0
INJURY_RESISTANCE_MAX=8

RECOVERY_SCORE=0
RECOVERY_MAX=10

MENTAL_SCORE=0
MENTAL_MAX=10

THERMOREGULATION_SCORE=0
THERMO_MAX=6

# Wyniki per dyscyplina
SWIM_SCORE=0
BIKE_SCORE=0
RUN_SCORE=0

# Szczegółowe wyniki
declare -a FINDINGS
declare -a STRENGTHS
declare -a WEAKNESSES
declare -a RECOMMENDATIONS

#===============================================================================
# SEKCJA 1: WYDAJNOŚĆ TLENOWA (VO2max)
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🫁 WYDAJNOŚĆ TLENOWA (VO2max potential)                                  ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "VO2max to kluczowy parametr w triathlonie - określa pułap tlenowy."
echo ""

# PPARGC1A (PGC-1α) - główny regulator biogenezy mitochondriów
echo "Analizuję PPARGC1A (biogeneza mitochondriów)..."
GT=$(check_snp "4:23814039")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ PPARGC1A rs8192678: G/G - Optymalna biogeneza mitochondriów"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 3))
            SWIM_SCORE=$((SWIM_SCORE + 1))
            BIKE_SCORE=$((BIKE_SCORE + 1))
            RUN_SCORE=$((RUN_SCORE + 1))
            FINDINGS+=("PPARGC1A G/G: Świetna zdolność do tworzenia nowych mitochondriów")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  🟡 PPARGC1A rs8192678: G/A - Standardowa biogeneza"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 2))
            ;;
        "1/1"|"1|1")
            echo "  🔴 PPARGC1A rs8192678: A/A - Słabsza adaptacja mitochondrialna"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 1))
            WEAKNESSES+=("Wolniejsza adaptacja mitochondrialna - wymaga więcej czasu na budowanie bazy tlenowej")
            RECOMMENDATIONS+=("Dłuższe bloki treningowe bazowe (12-16 tygodni) przed intensywnością")
            ;;
    esac
fi

# VEGFA - angiogeneza (tworzenie nowych naczyń)
echo "Analizuję VEGFA (angiogeneza)..."
GT=$(check_snp "6:43770613")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ VEGFA rs2010963: G/G - Wysoka produkcja VEGF"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 2))
            FINDINGS+=("VEGFA G/G: Lepsza kapilaryzacja mięśni = efektywniejszy transport tlenu")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  🟡 VEGFA rs2010963: G/C - Standardowa angiogeneza"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 1))
            ;;
        "1/1"|"1|1")
            echo "  🔴 VEGFA rs2010963: C/C - Niższa produkcja VEGF"
            WEAKNESSES+=("Słabsza angiogeneza - wolniejsze budowanie kapilaryzacji")
            ;;
    esac
fi

# HIF1A - adaptacja do hipoksji
echo "Analizuję HIF1A (adaptacja do hipoksji)..."
GT=$(check_snp "14:62207556")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  🟡 HIF1A rs11549465: C/C - Standardowa odpowiedź na hipoksję"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 1))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0"|"1/1"|"1|1")
            echo "  ✅ HIF1A rs11549465: T carrier - Lepsza adaptacja do hipoksji!"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 2))
            STRENGTHS+=("HIF1A: Lepsza adaptacja do niedoboru tlenu - przewaga w wysokich intensywnościach")
            RECOMMENDATIONS+=("Rozważ trening wysokościowy lub maski hipoksyczne")
            ;;
    esac
fi

# EPAS1 (HIF-2α) - adaptacja wysokościowa (wariant Denisowiański!)
echo "Analizuję EPAS1 (adaptacja wysokościowa)..."
GT=$(check_snp "2:46441523")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  🌟 EPAS1: Wariant adaptacji wysokościowej (Denisowiański!)"
    VO2MAX_SCORE=$((VO2MAX_SCORE + 2))
    STRENGTHS+=("EPAS1: Rzadki wariant adaptacji wysokościowej - potencjalna przewaga w endurance!")
fi

# NRF1 - funkcja mitochondriów
echo "Analizuję NRF1 (funkcja mitochondriów)..."
GT=$(check_snp "7:129613162")
if [ -n "$GT" ]; then
    if [[ "$GT" != "0/0" ]] && [[ "$GT" != "0|0" ]]; then
        echo "  ✅ NRF1 rs2402970: Wariant korzystny dla wydolności"
        VO2MAX_SCORE=$((VO2MAX_SCORE + 2))
    else
        echo "  🟡 NRF1 rs2402970: Standardowy"
        VO2MAX_SCORE=$((VO2MAX_SCORE + 1))
    fi
fi

# ADRB2 - receptory beta-adrenergiczne (bronchodilatacja, lipoliza)
echo "Analizuję ADRB2 (receptory adrenergiczne)..."
GT=$(check_snp "5:148826877")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ ADRB2 Gly16: Lepsza bronchodilatacja podczas wysiłku"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 2))
            FINDINGS+=("ADRB2 Gly16: Lepsze rozszerzanie oskrzeli = łatwiejszy oddech przy intensywnym wysiłku")
            ;;
        *)
            echo "  🟡 ADRB2 Arg16: Standardowa odpowiedź"
            VO2MAX_SCORE=$((VO2MAX_SCORE + 1))
            ;;
    esac
fi

echo ""
echo "   VO2max Score: $VO2MAX_SCORE / $VO2MAX_MAX"
printf "   ["
for ((i=0; i<VO2MAX_SCORE; i++)); do printf "█"; done
for ((i=VO2MAX_SCORE; i<VO2MAX_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 2: METABOLIZM TŁUSZCZÓW (kluczowe dla długich dystansów!)
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🔥 METABOLIZM TŁUSZCZÓW (Fat Oxidation)                                  ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "W Ironmanie spalasz ~8000-10000 kcal. Efektywne spalanie tłuszczów = klucz!"
echo ""

# PPARA - główny regulator metabolizmu tłuszczów
echo "Analizuję PPARA (metabolizm kwasów tłuszczowych)..."
GT=$(check_snp "22:46150500")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  🔴 PPARA rs4253778: G/G - Niższa oksydacja tłuszczów"
            FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 1))
            WEAKNESSES+=("PPARA G/G: Słabsze spalanie tłuszczów - większe uzależnienie od glikogenu")
            RECOMMENDATIONS+=("Trening nisko-glikogenowy (train low) 1-2x/tydzień buduje fat adaptation")
            RECOMMENDATIONS+=("Dieta ketogeniczna/LCHF może pomóc kompensować PPARA")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  🟡 PPARA rs4253778: G/C - Średnia oksydacja tłuszczów"
            FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 2))
            ;;
        "1/1"|"1|1")
            echo "  ✅ PPARA rs4253778: C/C - Wysoka oksydacja tłuszczów!"
            FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 3))
            STRENGTHS+=("PPARA C/C: Świetne spalanie tłuszczów - ogromna przewaga na długich dystansach!")
            BIKE_SCORE=$((BIKE_SCORE + 2))
            RUN_SCORE=$((RUN_SCORE + 2))
            ;;
    esac
fi

# PPARD - utylizacja tłuszczów w mięśniach
echo "Analizuję PPARD (utylizacja tłuszczów)..."
GT=$(check_snp "6:35381192")
if [ -n "$GT" ]; then
    if [[ "$GT" != "0/0" ]] && [[ "$GT" != "0|0" ]]; then
        echo "  ✅ PPARD rs2016520: Wariant korzystny - lepsza utylizacja tłuszczów"
        FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 2))
        FINDINGS+=("PPARD: Zwiększona ekspresja genów metabolizmu tłuszczów w mięśniach")
    else
        echo "  🟡 PPARD rs2016520: Standardowy"
        FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 1))
    fi
fi

# FABP2 - wchłanianie kwasów tłuszczowych
echo "Analizuję FABP2 (transport kwasów tłuszczowych)..."
GT=$(check_snp "4:119502441")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ FABP2 Ala54: Normalne wchłanianie tłuszczów"
            FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0"|"1/1"|"1|1")
            echo "  🟡 FABP2 Thr54: Zwiększone wchłanianie - uważaj na tłuszcze nasycone"
            FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 1))
            RECOMMENDATIONS+=("FABP2 Thr54: Priorytet tłuszcze nienasycone (oliwa, orzechy, awokado)")
            ;;
    esac
fi

# ADIPOQ - adiponektyna (wrażliwość na insulinę, spalanie tłuszczu)
echo "Analizuję ADIPOQ (adiponektyna)..."
GT=$(check_snp "3:186570892")
if [ -n "$GT" ]; then
    if [[ "$GT" == "0/0" ]] || [[ "$GT" == "0|0" ]]; then
        echo "  ✅ ADIPOQ rs1501299: G/G - Wyższa adiponektyna"
        FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 2))
        FINDINGS+=("ADIPOQ: Wysoka adiponektyna = lepsza wrażliwość na insulinę i spalanie tłuszczu")
    else
        echo "  🟡 ADIPOQ rs1501299: T carrier - Niższa adiponektyna"
        FAT_METABOLISM_SCORE=$((FAT_METABOLISM_SCORE + 1))
    fi
fi

echo ""
echo "   Fat Metabolism Score: $FAT_METABOLISM_SCORE / $FAT_METABOLISM_MAX"
printf "   ["
for ((i=0; i<FAT_METABOLISM_SCORE; i++)); do printf "█"; done
for ((i=FAT_METABOLISM_SCORE; i<FAT_METABOLISM_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 3: TOLERANCJA MLECZANU
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  ⚡ TOLERANCJA MLECZANU (Lactate Threshold)                               ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# MCT1 (SLC16A1) - transporter mleczanu
echo "Analizuję MCT1/SLC16A1 (transport mleczanu)..."
GT=$(check_snp "1:113454754")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ MCT1 rs1049434: A/A - Efektywny transport mleczanu"
            LACTATE_SCORE=$((LACTATE_SCORE + 3))
            STRENGTHS+=("MCT1 A/A: Szybkie usuwanie mleczanu z mięśni - możesz dłużej utrzymać wysoką intensywność")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  🟡 MCT1 rs1049434: A/T - Średni transport mleczanu"
            LACTATE_SCORE=$((LACTATE_SCORE + 2))
            ;;
        "1/1"|"1|1")
            echo "  🔴 MCT1 rs1049434: T/T - Wolniejszy transport mleczanu"
            LACTATE_SCORE=$((LACTATE_SCORE + 1))
            WEAKNESSES+=("MCT1 T/T: Wolniejsze usuwanie mleczanu - szybsze zakwaszenie")
            RECOMMENDATIONS+=("Więcej treningu progowego (tempo runs, sweet spot) dla kompensacji MCT1")
            ;;
    esac
fi

# AMPD1 - deaminaza AMP (metabolizm energetyczny)
echo "Analizuję AMPD1 (metabolizm energetyczny)..."
GT=$(check_snp "1:115227107")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ AMPD1 rs17602729: C/C - Normalna aktywność AMPD"
            LACTATE_SCORE=$((LACTATE_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  🟡 AMPD1 rs17602729: C/T - Obniżona aktywność (nosiciel)"
            LACTATE_SCORE=$((LACTATE_SCORE + 1))
            ;;
        "1/1"|"1|1")
            echo "  🔴 AMPD1 rs17602729: T/T - Niedobór AMPD (rzadki)"
            WEAKNESSES+=("AMPD1 T/T: Niedobór może powodować szybsze zmęczenie mięśni")
            ;;
    esac
fi

# CKM - kinaza kreatynowa
echo "Analizuję CKM (kinaza kreatynowa)..."
GT=$(check_snp "19:45812091")
if [ -n "$GT" ]; then
    if [[ "$GT" != "0/0" ]] && [[ "$GT" != "0|0" ]]; then
        echo "  ✅ CKM rs8111989: Wariant korzystny dla wydolności"
        LACTATE_SCORE=$((LACTATE_SCORE + 2))
    else
        echo "  🟡 CKM rs8111989: Standardowy"
        LACTATE_SCORE=$((LACTATE_SCORE + 1))
    fi
fi

echo ""
echo "   Lactate Score: $LACTATE_SCORE / $LACTATE_MAX"
printf "   ["
for ((i=0; i<LACTATE_SCORE; i++)); do printf "█"; done
for ((i=LACTATE_SCORE; i<LACTATE_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 4: TYP WŁÓKIEN MIĘŚNIOWYCH
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  💪 TYP WŁÓKIEN MIĘŚNIOWYCH                                               ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Triathlon wymaga głównie włókien wolnokurczliwych (Typ I) - wytrzymałość."
echo ""

# ACTN3 - kluczowy gen!
echo "Analizuję ACTN3 (typ włókien)..."
GT=$(check_snp "11:66560624")
ACTN3_STATUS=""
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            ACTN3_STATUS="RR"
            echo "  🔴 ACTN3 R577X: R/R - Przewaga włókien szybkich (Typ II)"
            echo "     → Lepszy dla sprintu, słabszy dla długich dystansów"
            ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 1))
            WEAKNESSES+=("ACTN3 R/R: Genetycznie predysponowany do siły/sprintu, nie ultra-wytrzymałości")
            RECOMMENDATIONS+=("Skup się na krótszych formatach: Sprint, Olympic triathlon")
            RECOMMENDATIONS+=("Dla Ironman: Więcej objętości treningowej niż osoby X/X")
            SWIM_SCORE=$((SWIM_SCORE + 2))  # R/R może być dobre dla pływania (siła!)
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            ACTN3_STATUS="RX"
            echo "  ✅ ACTN3 R577X: R/X - Mieszanka włókien - IDEALNE dla triathlonu!"
            echo "     → Dobra siła (pływanie) + wytrzymałość (rower, bieg)"
            ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 3))
            STRENGTHS+=("ACTN3 R/X: Optymalny genotyp dla triathlonu! Siła + wytrzymałość")
            SWIM_SCORE=$((SWIM_SCORE + 2))
            BIKE_SCORE=$((BIKE_SCORE + 2))
            RUN_SCORE=$((RUN_SCORE + 2))
            ;;
        "1/1"|"1|1")
            ACTN3_STATUS="XX"
            echo "  ✅ ACTN3 R577X: X/X - Przewaga włókien wolnych (Typ I)"
            echo "     → Świetne dla długich dystansów, słabsze dla sprintu"
            ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 3))
            STRENGTHS+=("ACTN3 X/X: Genetycznie stworzony do ultra-wytrzymałości!")
            BIKE_SCORE=$((BIKE_SCORE + 3))
            RUN_SCORE=$((RUN_SCORE + 3))
            RECOMMENDATIONS+=("ACTN3 X/X: Half/Full Ironman to Twoja naturalna dyscyplina!")
            ;;
    esac
fi

# ACE I/D
echo "Analizuję ACE (układ renina-angiotensyna)..."
GT=$(check_snp "17:63488529")
ACE_STATUS=""
if [ -n "$GT" ]; then
    case $GT in
        "1/1"|"1|1")
            ACE_STATUS="II"
            echo "  ✅ ACE I/D: I/I - Genotyp wytrzymałościowy"
            ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 3))
            STRENGTHS+=("ACE I/I: Świetna wydolność tlenowa, lepsza kapilaryzacja")
            BIKE_SCORE=$((BIKE_SCORE + 2))
            RUN_SCORE=$((RUN_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            ACE_STATUS="ID"
            echo "  🟡 ACE I/D: I/D - Zbalansowany"
            ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 2))
            ;;
        "0/0"|"0|0")
            ACE_STATUS="DD"
            echo "  🔴 ACE I/D: D/D - Genotyp siłowy"
            ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 1))
            WEAKNESSES+=("ACE D/D: Lepszy dla siły/sprintu niż ultra-wytrzymałości")
            ;;
    esac
fi

# ACVR1B - hipertrofia mięśni
echo "Analizuję ACVR1B (odpowiedź mięśni na trening)..."
GT=$(check_snp "12:52372642")
if [ -n "$GT" ]; then
    if [[ "$GT" != "0/0" ]] && [[ "$GT" != "0|0" ]]; then
        echo "  ✅ ACVR1B rs2854464: Lepsza odpowiedź na trening wytrzymałościowy"
        ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 2))
    else
        echo "  🟡 ACVR1B rs2854464: Standardowa odpowiedź"
        ENDURANCE_FIBER_SCORE=$((ENDURANCE_FIBER_SCORE + 1))
    fi
fi

echo ""
echo "   Endurance Fiber Score: $ENDURANCE_FIBER_SCORE / $ENDURANCE_FIBER_MAX"
printf "   ["
for ((i=0; i<ENDURANCE_FIBER_SCORE; i++)); do printf "█"; done
for ((i=ENDURANCE_FIBER_SCORE; i<ENDURANCE_FIBER_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 5: ODPORNOŚĆ NA KONTUZJE
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🦴 ODPORNOŚĆ NA KONTUZJE                                                 ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Triathlon = ogromne obciążenie ścięgien i stawów. Kontuzja = koniec sezonu."
echo ""

# COL5A1 - kolagen typ V (ścięgna, więzadła)
echo "Analizuję COL5A1 (elastyczność ścięgien)..."
GT=$(check_snp "9:137684151")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ COL5A1 rs12722: C/C - Sztywniejsze, ale mocniejsze ścięgna"
            INJURY_RESISTANCE_SCORE=$((INJURY_RESISTANCE_SCORE + 2))
            RUN_SCORE=$((RUN_SCORE + 1))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  🟡 COL5A1 rs12722: C/T - Pośredni"
            INJURY_RESISTANCE_SCORE=$((INJURY_RESISTANCE_SCORE + 1))
            ;;
        "1/1"|"1|1")
            echo "  🔴 COL5A1 rs12722: T/T - Elastyczne ścięgna, wyższe ryzyko kontuzji"
            WEAKNESSES+=("COL5A1 T/T: Wyższe ryzyko urazu ścięgna Achillesa!")
            RECOMMENDATIONS+=("PREWENCJA: Ćwiczenia ekscentryczne na Achilles, wolniejsza progresja biegu")
            ;;
    esac
fi

# COL1A1 - kolagen typ I
echo "Analizuję COL1A1 (jakość kolagenu)..."
GT=$(check_snp "17:50201632")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ COL1A1 rs1800012: G/G - Mocniejszy kolagen"
            INJURY_RESISTANCE_SCORE=$((INJURY_RESISTANCE_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0"|"1/1"|"1|1")
            echo "  🔴 COL1A1 rs1800012: T carrier - Słabszy kolagen"
            INJURY_RESISTANCE_SCORE=$((INJURY_RESISTANCE_SCORE + 1))
            RECOMMENDATIONS+=("KOLAGEN: Suplementacja 15g + wit. C przed treningiem")
            ;;
    esac
fi

# MMP3 - degradacja macierzy
echo "Analizuję MMP3 (degradacja macierzy)..."
GT=$(check_snp "11:102733807")
if [ -n "$GT" ]; then
    if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
        echo "  🔴 MMP3 5A/5A: Wyższe ryzyko tendinopatii Achillesa"
        WEAKNESSES+=("MMP3 5A/5A: Zwiększone ryzyko tendinopatii - monitoruj Achilles!")
    else
        echo "  🟡 MMP3: Standardowe ryzyko"
        INJURY_RESISTANCE_SCORE=$((INJURY_RESISTANCE_SCORE + 2))
    fi
fi

# GDF5 - rozwój chrząstki
echo "Analizuję GDF5 (zdrowie chrząstki)..."
GT=$(check_snp "20:35437977")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ GDF5 rs143383: G/G - Lepsza chrząstka"
            INJURY_RESISTANCE_SCORE=$((INJURY_RESISTANCE_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0"|"1/1"|"1|1")
            echo "  🔴 GDF5 rs143383: A carrier - Wyższe ryzyko osteoartrozy"
            WEAKNESSES+=("GDF5: Wyższe ryzyko problemów z chrząstką - dbaj o stawy!")
            RECOMMENDATIONS+=("STAWY: Glucosamine 1500mg + Chondroitin 1200mg + UC-II 40mg/dzień")
            ;;
    esac
fi

echo ""
echo "   Injury Resistance Score: $INJURY_RESISTANCE_SCORE / $INJURY_RESISTANCE_MAX"
printf "   ["
for ((i=0; i<INJURY_RESISTANCE_SCORE; i++)); do printf "█"; done
for ((i=INJURY_RESISTANCE_SCORE; i<INJURY_RESISTANCE_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 6: REGENERACJA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🔄 REGENERACJA                                                           ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Triathlon = 3 dyscypliny = potrójna kumulacja zmęczenia. Regeneracja kluczowa!"
echo ""

# IL6 - zapalenie
echo "Analizuję IL6 (odpowiedź zapalna)..."
GT=$(check_snp "7:22727026")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ IL6 rs1800795: G/G - Niska odpowiedź zapalna"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 3))
            STRENGTHS+=("IL6 G/G: Szybka regeneracja, mniejszy DOMS po ciężkich sesjach")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  🟡 IL6 rs1800795: G/C - Średnia odpowiedź"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
            ;;
        "1/1"|"1|1")
            echo "  🔴 IL6 rs1800795: C/C - Wysoka odpowiedź zapalna"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            WEAKNESSES+=("IL6 C/C: Wolniejsza regeneracja, więcej DOMS")
            RECOMMENDATIONS+=("REGENERACJA: Minimum 48-72h między ciężkimi sesjami")
            RECOMMENDATIONS+=("ANTYOKSYDANTY: Omega-3 (3-4g), kurkumina, tart cherry juice")
            ;;
    esac
fi

# TNF - czynnik martwicy nowotworów
echo "Analizuję TNF (czynnik zapalny)..."
GT=$(check_snp "6:31575254")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ TNF rs1800629: G/G - Niższa produkcja TNF-α"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
            ;;
        *)
            echo "  🟡 TNF rs1800629: A carrier - Wyższa produkcja TNF-α"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            RECOMMENDATIONS+=("TNF: Więcej snu (8-9h) i aktywna regeneracja (foam rolling, masaż)")
            ;;
    esac
fi

# CRP - białko C-reaktywne
echo "Analizuję CRP (marker zapalenia)..."
GT=$(check_snp "1:159712267")
if [ -n "$GT" ]; then
    if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
        echo "  🔴 CRP rs1205: T/T - Wyższy bazowy poziom CRP"
        RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
        RECOMMENDATIONS+=("CRP: Regularne badania CRP do monitorowania przeciążenia")
    else
        echo "  ✅ CRP rs1205: C carrier - Niższy bazowy poziom"
        RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
    fi
fi

# SOD2 - stres oksydacyjny
echo "Analizuję SOD2 (ochrona antyoksydacyjna)..."
GT=$(check_snp "6:160113872")
if [ -n "$GT" ]; then
    case $GT in
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  ✅ SOD2 Val16Ala: Val/Ala - Optymalna ochrona"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
            ;;
        "0/0"|"0|0")
            echo "  🟡 SOD2 Val16Ala: Val/Val - Niższa aktywność w mitochondriach"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            RECOMMENDATIONS+=("SOD2 Val/Val: NAC 600-1200mg/dzień dla wsparcia antyoksydacyjnego")
            ;;
        "1/1"|"1|1")
            echo "  🟡 SOD2 Val16Ala: Ala/Ala - Wyższa aktywność, ale potrzebuje manganu"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            RECOMMENDATIONS+=("SOD2 Ala/Ala: Upewnij się, że masz dość manganu w diecie")
            ;;
    esac
fi

echo ""
echo "   Recovery Score: $RECOVERY_SCORE / $RECOVERY_MAX"
printf "   ["
for ((i=0; i<RECOVERY_SCORE; i++)); do printf "█"; done
for ((i=RECOVERY_SCORE; i<RECOVERY_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 7: MENTALNOŚĆ I ODPORNOŚĆ PSYCHICZNA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🧠 MENTALNOŚĆ I ODPORNOŚĆ PSYCHICZNA                                     ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Ironman to 8-17h cierpienia. Głowa decyduje o wszystkim!"
echo ""

# COMT - kluczowy dla wydajności pod presją
echo "Analizuję COMT (odporność na stres)..."
GT=$(check_snp "22:19951271")
COMT_STATUS=""
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            COMT_STATUS="warrior"
            echo "  ⚔️  COMT Val158Met: Val/Val - WARRIOR"
            echo "     → Świetny pod presją, szybki klirens dopaminy"
            MENTAL_SCORE=$((MENTAL_SCORE + 3))
            STRENGTHS+=("COMT Warrior: Świetna wydajność w wyścigu, odporny na stres startowy")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            COMT_STATUS="balanced"
            echo "  ⚖️  COMT Val158Met: Val/Met - Zbalansowany"
            MENTAL_SCORE=$((MENTAL_SCORE + 2))
            ;;
        "1/1"|"1|1")
            COMT_STATUS="worrier"
            echo "  🧘 COMT Val158Met: Met/Met - WORRIER"
            echo "     → Lepsze planowanie, ale wrażliwszy na stres"
            MENTAL_SCORE=$((MENTAL_SCORE + 1))
            WEAKNESSES+=("COMT Worrier: Może 'choke' pod presją - potrzebujesz więcej przygotowania mentalnego")
            RECOMMENDATIONS+=("MENTAL: Wizualizacja, rutyny przedstartowe, praca z psychologiem sportu")
            RECOMMENDATIONS+=("COMT Met/Met: Ogranicz kofeinę w dniu wyścigu - może zwiększyć niepokój")
            ;;
    esac
fi

# BDNF - plastyczność mózgu
echo "Analizuję BDNF (plastyczność mózgu)..."
GT=$(check_snp "11:27679916")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ BDNF Val66Val: Wysoki BDNF - dobra neuroplastyczność"
            MENTAL_SCORE=$((MENTAL_SCORE + 2))
            ;;
        *)
            echo "  🟡 BDNF Met carrier: Niższy BDNF"
            MENTAL_SCORE=$((MENTAL_SCORE + 1))
            RECOMMENDATIONS+=("BDNF Met: Regularny trening aerobowy PODNOSI BDNF - to Twoja naturalna terapia!")
            ;;
    esac
fi

# DRD2 - motywacja i nagroda
echo "Analizuję DRD2 (układ nagrody)..."
GT=$(check_snp "11:113400106")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ DRD2 Taq1A: C/C - Normalna gęstość receptorów dopaminy"
            MENTAL_SCORE=$((MENTAL_SCORE + 2))
            ;;
        *)
            echo "  🟡 DRD2 Taq1A: T carrier - Mniej receptorów dopaminy"
            MENTAL_SCORE=$((MENTAL_SCORE + 1))
            FINDINGS+=("DRD2 T carrier: Możesz potrzebować większych celów/wyzwań dla motywacji")
            RECOMMENDATIONS+=("CELE: Stawiaj ambitne, mierzalne cele - potrzebujesz większej stymulacji")
            ;;
    esac
fi

# 5-HTTLPR (SLC6A4) - serotonina i nastrój
echo "Analizuję SLC6A4 (transporter serotoniny)..."
GT=$(check_snp "17:30194319")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ✅ SLC6A4: L/L - Stabilniejszy nastrój"
            MENTAL_SCORE=$((MENTAL_SCORE + 2))
            ;;
        *)
            echo "  🟡 SLC6A4: S carrier - Większa wrażliwość emocjonalna"
            MENTAL_SCORE=$((MENTAL_SCORE + 1))
            RECOMMENDATIONS+=("NASTRÓJ: Rytuały poranne, regularny sen, ekspozycja na światło - stabilizują nastrój")
            ;;
    esac
fi

echo ""
echo "   Mental Score: $MENTAL_SCORE / $MENTAL_MAX"
printf "   ["
for ((i=0; i<MENTAL_SCORE; i++)); do printf "█"; done
for ((i=MENTAL_SCORE; i<MENTAL_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 8: TERMOREGULACJA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🌡️ TERMOREGULACJA                                                        ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Ironman Kona = 35°C + wilgotność. Termoregulacja może wygrać lub przegrać wyścig!"
echo ""

# UCP2 - termogeneza
echo "Analizuję UCP2 (termogeneza)..."
GT=$(check_snp "11:73976705")
if [ -n "$GT" ]; then
    if [[ "$GT" != "0/0" ]] && [[ "$GT" != "0|0" ]]; then
        echo "  ✅ UCP2 rs660339: Wariant - Lepsza regulacja temperatury"
        THERMOREGULATION_SCORE=$((THERMOREGULATION_SCORE + 2))
    else
        echo "  🟡 UCP2 rs660339: Standardowy"
        THERMOREGULATION_SCORE=$((THERMOREGULATION_SCORE + 1))
    fi
fi

# UCP3 - termogeneza mięśniowa
echo "Analizuję UCP3 (termogeneza mięśniowa)..."
GT=$(check_snp "11:74009381")
if [ -n "$GT" ]; then
    if [[ "$GT" != "0/0" ]] && [[ "$GT" != "0|0" ]]; then
        echo "  ✅ UCP3 rs1800849: Wariant - Lepsza dysypacja ciepła"
        THERMOREGULATION_SCORE=$((THERMOREGULATION_SCORE + 2))
        STRENGTHS+=("UCP3: Lepsza tolerancja ciepła - przewaga w gorących warunkach")
    else
        echo "  🟡 UCP3 rs1800849: Standardowy"
        THERMOREGULATION_SCORE=$((THERMOREGULATION_SCORE + 1))
    fi
fi

# TRPM8 - wrażliwość na zimno
echo "Analizuję TRPM8 (wrażliwość na zimno)..."
GT=$(check_snp "2:234065133")
if [ -n "$GT" ]; then
    if [[ "$GT" != "0/0" ]] && [[ "$GT" != "0|0" ]]; then
        echo "  ❄️  TRPM8 rs11563208: Wyższa wrażliwość na zimno"
        RECOMMENDATIONS+=("ZIMNO: Możesz być wrażliwszy na zimną wodę - rozważ grubszy pianka")
    else
        echo "  ✅ TRPM8 rs11563208: Normalna tolerancja zimna"
        THERMOREGULATION_SCORE=$((THERMOREGULATION_SCORE + 2))
    fi
fi

echo ""
echo "   Thermoregulation Score: $THERMOREGULATION_SCORE / $THERMO_MAX"
printf "   ["
for ((i=0; i<THERMOREGULATION_SCORE; i++)); do printf "█"; done
for ((i=THERMOREGULATION_SCORE; i<THERMO_MAX; i++)); do printf "░"; done
printf "]\n"
echo ""

#===============================================================================
# SEKCJA 9: WYNIKI PER DYSCYPLINA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  📊 PREDYSPOZYCJE DO POSZCZEGÓLNYCH DYSCYPLIN                             ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Dodaj dodatkowe punkty na podstawie ogólnych wyników
SWIM_SCORE=$((SWIM_SCORE + VO2MAX_SCORE / 3))
BIKE_SCORE=$((BIKE_SCORE + VO2MAX_SCORE / 3 + FAT_METABOLISM_SCORE / 2 + LACTATE_SCORE / 2))
RUN_SCORE=$((RUN_SCORE + VO2MAX_SCORE / 3 + INJURY_RESISTANCE_SCORE / 2 + RECOVERY_SCORE / 3))

SWIM_MAX=15
BIKE_MAX=20
RUN_MAX=20

echo "   🏊 PŁYWANIE:  $SWIM_SCORE / $SWIM_MAX"
printf "      ["
for ((i=0; i<SWIM_SCORE && i<SWIM_MAX; i++)); do printf "█"; done
for ((i=SWIM_SCORE; i<SWIM_MAX; i++)); do printf "░"; done
printf "]\n"

echo "   🚴 KOLARSTWO: $BIKE_SCORE / $BIKE_MAX"
printf "      ["
for ((i=0; i<BIKE_SCORE && i<BIKE_MAX; i++)); do printf "█"; done
for ((i=BIKE_SCORE; i<BIKE_MAX; i++)); do printf "░"; done
printf "]\n"

echo "   🏃 BIEGANIE:  $RUN_SCORE / $RUN_MAX"
printf "      ["
for ((i=0; i<RUN_SCORE && i<RUN_MAX; i++)); do printf "█"; done
for ((i=RUN_SCORE; i<RUN_MAX; i++)); do printf "░"; done
printf "]\n"

echo ""

# Określ najsilniejszą i najsłabszą dyscyplinę
STRONGEST=""
WEAKEST=""
MAX_SCORE=0
MIN_SCORE=999

if [ $SWIM_SCORE -gt $MAX_SCORE ]; then MAX_SCORE=$SWIM_SCORE; STRONGEST="🏊 Pływanie"; fi
if [ $BIKE_SCORE -gt $MAX_SCORE ]; then MAX_SCORE=$BIKE_SCORE; STRONGEST="🚴 Kolarstwo"; fi
if [ $RUN_SCORE -gt $MAX_SCORE ]; then MAX_SCORE=$RUN_SCORE; STRONGEST="🏃 Bieganie"; fi

if [ $SWIM_SCORE -lt $MIN_SCORE ]; then MIN_SCORE=$SWIM_SCORE; WEAKEST="🏊 Pływanie"; fi
if [ $BIKE_SCORE -lt $MIN_SCORE ]; then MIN_SCORE=$BIKE_SCORE; WEAKEST="🚴 Kolarstwo"; fi
if [ $RUN_SCORE -lt $MIN_SCORE ]; then MIN_SCORE=$RUN_SCORE; WEAKEST="🏃 Bieganie"; fi

echo "   Twoja najsilniejsza dyscyplina: $STRONGEST"
echo "   Dyscyplina do poprawy: $WEAKEST"
echo ""

#===============================================================================
# SEKCJA 10: CAŁKOWITY WYNIK I REKOMENDOWANY FORMAT
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                       ║"
echo "║                         CAŁKOWITY WYNIK TRIATLONOWY                                   ║"
echo "║                                                                                       ║"
echo "╚═══════════════════════════════════════════════════════════════════════════════════════╝"
echo ""

TOTAL_SCORE=$((VO2MAX_SCORE + FAT_METABOLISM_SCORE + LACTATE_SCORE + ENDURANCE_FIBER_SCORE + INJURY_RESISTANCE_SCORE + RECOVERY_SCORE + MENTAL_SCORE + THERMOREGULATION_SCORE))
TOTAL_MAX=$((VO2MAX_MAX + FAT_METABOLISM_MAX + LACTATE_MAX + ENDURANCE_FIBER_MAX + INJURY_RESISTANCE_MAX + RECOVERY_MAX + MENTAL_MAX + THERMO_MAX))

PERCENT=$((TOTAL_SCORE * 100 / TOTAL_MAX))

echo "   ┌─────────────────────────────────────────────────────────────────────┐"
echo "   │                                                                     │"
printf "   │   CAŁKOWITY WYNIK: %3d / %2d punktów (%d%%)                         │\n" $TOTAL_SCORE $TOTAL_MAX $PERCENT
echo "   │                                                                     │"
printf "   │   ["
for ((i=0; i<TOTAL_SCORE*50/TOTAL_MAX; i++)); do printf "█"; done
for ((i=TOTAL_SCORE*50/TOTAL_MAX; i<50; i++)); do printf "░"; done
printf "]   │\n"
echo "   │                                                                     │"
echo "   └─────────────────────────────────────────────────────────────────────┘"
echo ""

# Określ rekomendowany format
echo "   📋 REKOMENDOWANY FORMAT TRIATHLONU:"
echo ""

if [ $PERCENT -ge 75 ]; then
    echo "   🏆 ELITARNY POTENCJAŁ TRIATLONOWY!"
    if [[ "$ACTN3_STATUS" == "XX" ]] || [[ "$ACE_STATUS" == "II" ]]; then
        echo "   → Twoja genetyka wskazuje na IRONMAN / Half-Ironman"
        echo "   → Masz naturalną przewagę w ultra-wytrzymałości"
    else
        echo "   → Świetne predyspozycje do wszystkich formatów"
        echo "   → Twoja wszechstronność pozwala na Olympic lub 70.3"
    fi
elif [ $PERCENT -ge 60 ]; then
    echo "   ✅ DOBRY POTENCJAŁ TRIATLONOWY"
    if [[ "$ACTN3_STATUS" == "RR" ]]; then
        echo "   → Sprint / Olympic - lepiej wykorzystasz genetykę siłową"
    elif [ $FAT_METABOLISM_SCORE -ge 7 ]; then
        echo "   → Half-Ironman / Ironman - świetny metabolizm tłuszczów"
    else
        echo "   → Olympic distance - dobry balans wymagań"
    fi
elif [ $PERCENT -ge 45 ]; then
    echo "   🟡 ŚREDNI POTENCJAŁ - ALE MOŻLIWY SUKCES Z TRENINGIEM!"
    echo "   → Sprint / Olympic - skup się na słabych stronach"
    echo "   → Dłuższe dystanse wymagają więcej pracy"
else
    echo "   🔴 GENETYKA NIE SPRZYJA TRIATHLONOWI"
    echo "   → Ale genetyka to nie wszystko! Trening może wiele nadrobić"
    echo "   → Zacznij od Sprint distance i buduj stopniowo"
fi

echo ""

#===============================================================================
# SEKCJA 11: MOCNE I SŁABE STRONY
#===============================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "💪 TWOJE GENETYCZNE MOCNE STRONY:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
for strength in "${STRENGTHS[@]}"; do
    echo "   ✅ $strength"
done
if [ ${#STRENGTHS[@]} -eq 0 ]; then
    echo "   (Analiza nie wykryła wyróżniających się mocnych stron)"
fi
echo ""

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "⚠️  OBSZARY DO PRACY:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
for weakness in "${WEAKNESSES[@]}"; do
    echo "   ⚠️  $weakness"
done
if [ ${#WEAKNESSES[@]} -eq 0 ]; then
    echo "   (Brak znaczących słabości genetycznych)"
fi
echo ""

#===============================================================================
# SEKCJA 12: SPERSONALIZOWANE REKOMENDACJE
#===============================================================================

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📋 SPERSONALIZOWANE REKOMENDACJE TRENINGOWE:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
for rec in "${RECOMMENDATIONS[@]}"; do
    echo "   → $rec"
    echo ""
done

# Dodaj ogólne rekomendacje
echo "   → OBJĘTOŚĆ: Bazuj plan na $([ $RECOVERY_SCORE -ge 7 ] && echo "wysokiej" || echo "umiarkowanej") objętości"
echo ""
echo "   → INTENSYWNOŚĆ: $([ "$COMT_STATUS" == "warrior" ] && echo "Dobrze tolerujesz intensywne sesje" || echo "Ostrożnie z intensywnością - planuj ją strategicznie")"
echo ""
echo "   → REGENERACJA: Planuj $([ $RECOVERY_SCORE -ge 7 ] && echo "1 dzień" || echo "2 dni") odpoczynku po długich/ciężkich sesjach"
echo ""

#===============================================================================
# ZAPIS RAPORTU
#===============================================================================

REPORT="$OUTDIR/TRIATHLON_REPORT_$(date +%Y%m%d).txt"

{
echo "================================================================================"
echo "              RAPORT GENETYCZNEJ ANALIZY TRIATLONOWEJ"
echo "              Wygenerowano: $(date)"
echo "================================================================================"
echo ""
echo "PROFIL GENETYCZNY:"
echo "  ACTN3: $ACTN3_STATUS"
echo "  ACE: $ACE_STATUS"
echo "  COMT: $COMT_STATUS"
echo ""
echo "WYNIKI:"
echo "  VO2max potential:    $VO2MAX_SCORE / $VO2MAX_MAX"
echo "  Fat Metabolism:      $FAT_METABOLISM_SCORE / $FAT_METABOLISM_MAX"
echo "  Lactate Threshold:   $LACTATE_SCORE / $LACTATE_MAX"
echo "  Endurance Fibers:    $ENDURANCE_FIBER_SCORE / $ENDURANCE_FIBER_MAX"
echo "  Injury Resistance:   $INJURY_RESISTANCE_SCORE / $INJURY_RESISTANCE_MAX"
echo "  Recovery:            $RECOVERY_SCORE / $RECOVERY_MAX"
echo "  Mental:              $MENTAL_SCORE / $MENTAL_MAX"
echo "  Thermoregulation:    $THERMOREGULATION_SCORE / $THERMO_MAX"
echo ""
echo "  CAŁKOWITY:           $TOTAL_SCORE / $TOTAL_MAX ($PERCENT%)"
echo ""
echo "PER DYSCYPLINA:"
echo "  Pływanie:  $SWIM_SCORE"
echo "  Kolarstwo: $BIKE_SCORE"
echo "  Bieganie:  $RUN_SCORE"
echo ""
echo "MOCNE STRONY:"
for s in "${STRENGTHS[@]}"; do echo "  - $s"; done
echo ""
echo "SŁABE STRONY:"
for w in "${WEAKNESSES[@]}"; do echo "  - $w"; done
echo ""
echo "REKOMENDACJE:"
for r in "${RECOMMENDATIONS[@]}"; do echo "  - $r"; done
echo ""
echo "================================================================================"
} > "$REPORT"

echo ""
echo "📁 Pełny raport zapisany: $REPORT"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🏊🚴🏃 TRIATHLON GENETIC ANALYSIS COMPLETE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
