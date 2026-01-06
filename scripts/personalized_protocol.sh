#!/bin/bash
#===============================================================================
#
#    ██████╗ ███████╗██████╗ ███████╗ ██████╗ ███╗   ██╗ █████╗ ██╗     ██╗███████╗███████╗██████╗ 
#    ██╔══██╗██╔════╝██╔══██╗██╔════╝██╔═══██╗████╗  ██║██╔══██╗██║     ██║╚══███╔╝██╔════╝██╔══██╗
#    ██████╔╝█████╗  ██████╔╝███████╗██║   ██║██╔██╗ ██║███████║██║     ██║  ███╔╝ █████╗  ██║  ██║
#    ██╔═══╝ ██╔══╝  ██╔══██╗╚════██║██║   ██║██║╚██╗██║██╔══██║██║     ██║ ███╔╝  ██╔══╝  ██║  ██║
#    ██║     ███████╗██║  ██║███████║╚██████╔╝██║ ╚████║██║  ██║███████╗██║███████╗███████╗██████╔╝
#    ╚═╝     ╚══════╝╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝╚══════╝╚══════╝╚═════╝ 
#
#    ██████╗ ██████╗  ██████╗ ████████╗ ██████╗  ██████╗ ██████╗ ██╗     
#    ██╔══██╗██╔══██╗██╔═══██╗╚══██╔══╝██╔═══██╗██╔════╝██╔═══██╗██║     
#    ██████╔╝██████╔╝██║   ██║   ██║   ██║   ██║██║     ██║   ██║██║     
#    ██╔═══╝ ██╔══██╗██║   ██║   ██║   ██║   ██║██║     ██║   ██║██║     
#    ██║     ██║  ██║╚██████╔╝   ██║   ╚██████╔╝╚██████╗╚██████╔╝███████╗
#    ╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝  ╚═════╝ ╚═════╝ ╚══════╝
#
#    SPERSONALIZOWANY PROTOKOL TRENINGOWY, SUPLEMENTACYJNY I ZYWIENIOWY
#    Na podstawie analizy genetycznej
#
#    HELIXIGHT - Precision Genetic Fitness
#
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"
OUTDIR="personalized_protocol"
mkdir -p "$OUTDIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                       ║"
echo "║            SPERSONALIZOWANY PROTOKOL GENETYCZNY                                       ║"
echo "║            Trening • Suplementacja • Zywienie • Styl zycia                            ║"
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

# Inicjalizacja zmiennych profilowych
POWER_SCORE=0
ENDURANCE_SCORE=0
RECOVERY_SPEED="normal"
INJURY_RISK="normal"
STRESS_TYPE="balanced"
CAFFEINE_METABOLISM="normal"
FAT_SENSITIVITY="normal"
CARB_SENSITIVITY="normal"
INFLAMMATION_LEVEL="normal"

# Rekomendacje (będą wypełniane)
declare -a TRAINING_RECS
declare -a SUPPLEMENT_RECS
declare -a NUTRITION_RECS
declare -a AVOID_RECS
declare -a LIFESTYLE_RECS
declare -a HORMONE_NOTES

#===============================================================================
# SEKCJA 1: ANALIZA GENOW TRENINGOWYCH
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  💪 ANALIZA PROFILU TRENINGOWEGO                                          ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# ACTN3 - kluczowy gen typu wlokien
echo "Analizuje ACTN3 (typ wlokien miesniowych)..."
GT=$(check_snp "11:66560624")
ACTN3_STATUS=""
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            ACTN3_STATUS="RR"
            POWER_SCORE=$((POWER_SCORE + 3))
            echo "  ACTN3: R/R → 🏋️ ELITARNA genetyka silowa!"
            TRAINING_RECS+=("SILA: Twoj genotyp ACTN3 R/R daje przewage w sile eksplozywnej. Priorytet: ciezkie treningi silowe 3-4x/tydzien.")
            TRAINING_RECS+=("SPRINT: Masz predyspozycje do sprintu i skoków. Wlacz trening plyometryczny.")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            ACTN3_STATUS="RX"
            POWER_SCORE=$((POWER_SCORE + 2))
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            echo "  ACTN3: R/X → ⚖️ Dobra mieszanka sila/wytrzymalosc"
            TRAINING_RECS+=("HYBRYDA: Twoj genotyp ACTN3 R/X pozwala na skuteczny trening zarowno silowy jak i wytrzymalosciowy.")
            ;;
        "1/1"|"1|1")
            ACTN3_STATUS="XX"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 3))
            echo "  ACTN3: X/X → 🏃 Genetyka wytrzymalosciowa"
            TRAINING_RECS+=("WYTRZYMALOSC: Twoj genotyp ACTN3 X/X daje przewage w sportach wytrzymalosciowych. Priorytet: bieganie, kolarstwo, plywanie.")
            TRAINING_RECS+=("OBJETOSC: Mozesz tolerowac wieksza objetosc treningowa niz osoby R/R.")
            ;;
    esac
fi

# ACE I/D - wytrzymalosc vs sila
echo "Analizuje ACE (uklad renina-angiotensyna)..."
GT=$(check_snp "17:63488529")
if [ -n "$GT" ]; then
    case $GT in
        "1/1"|"1|1")
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 2))
            echo "  ACE: I/I → Lepsza wydajnosc tlenowa"
            TRAINING_RECS+=("TLEN: ACE I/I = lepsza kapilaryzacja. Wlacz dlugie sesje cardio (60-90 min) do programu.")
            ;;
        "0/0"|"0|0")
            POWER_SCORE=$((POWER_SCORE + 2))
            echo "  ACE: D/D → Lepsza sila i hipertrofia"
            TRAINING_RECS+=("HIPERTROFIA: ACE D/D = lepsza odpowiedz hipertroficzna. Skup sie na progresji ciezaru.")
            ;;
    esac
fi

# PPARGC1A - biogeneza mitochondriow
echo "Analizuje PPARGC1A (mitochondria)..."
GT=$(check_snp "4:23814039")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            echo "  PPARGC1A: G/G → Dobra biogeneza mitochondriow"
            ;;
        "1/1"|"1|1")
            echo "  PPARGC1A: A/A → Slabsza adaptacja wytrzymalosciowa"
            TRAINING_RECS+=("MITOCHONDRIA: PPARGC1A A/A - potrzebujesz wiecej czasu na adaptacje wytrzymalosciowa. Stopniowo zwiększaj objetosc.")
            SUPPLEMENT_RECS+=("CoQ10: 100-200mg/dzien - wspiera funkcje mitochondriow (szczegolnie wazne przy Twoim PPARGC1A)")
            ;;
    esac
fi

# IL6 - regeneracja i zapalenie
echo "Analizuje IL6 (regeneracja)..."
GT=$(check_snp "7:22727026")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            RECOVERY_SPEED="fast"
            INFLAMMATION_LEVEL="low"
            echo "  IL6: G/G → Szybka regeneracja, niskie zapalenie"
            TRAINING_RECS+=("CZESTOTLIWOSC: IL6 G/G = szybsza regeneracja. Mozesz trenowac te same partie co 48h.")
            ;;
        "1/1"|"1|1")
            RECOVERY_SPEED="slow"
            INFLAMMATION_LEVEL="high"
            echo "  IL6: C/C → Wolniejsza regeneracja, wiecej DOMS"
            TRAINING_RECS+=("REGENERACJA: IL6 C/C = wolniejsza regeneracja. Potrzebujesz 72h+ miedzy ciezkimi sesjami.")
            TRAINING_RECS+=("DELOAD: Planuj deload co 3-4 tygodnie (nie 4-6 jak standardowo).")
            SUPPLEMENT_RECS+=("OMEGA-3: 3-4g EPA+DHA/dzien - kluczowe dla redukcji zapalenia przy Twoim IL6")
            SUPPLEMENT_RECS+=("KURKUMINA: 500-1000mg/dzien z piperyną - naturalne działanie przeciwzapalne")
            ;;
    esac
fi

# SOD2 - stres oksydacyjny
echo "Analizuje SOD2 (antyoksydanty)..."
GT=$(check_snp "6:160113872")
if [ -n "$GT" ]; then
    case $GT in
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  SOD2: Val/Ala → Optymalna ochrona antyoksydacyjna"
            ;;
        "0/0"|"0|0")
            echo "  SOD2: Val/Val → Nizsza aktywnosc w mitochondriach"
            SUPPLEMENT_RECS+=("NAC: 600-1200mg/dzien - wspiera glutation przy Twoim SOD2 Val/Val")
            ;;
        "1/1"|"1|1")
            echo "  SOD2: Ala/Ala → Wyzsza aktywnosc, ale wrazliwosc na niedobory"
            SUPPLEMENT_RECS+=("MANGAN: Upewnij sie, ze dieta dostarcza mangan (orzechy, nasiona) - kofaktor SOD2")
            ;;
    esac
fi

# COL1A1 i COL5A1 - ryzyko kontuzji
echo "Analizuje geny kolagenu (kontuzje)..."
GT=$(check_snp "17:50201632")
COL_RISK=0
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    COL_RISK=$((COL_RISK + 1))
fi

GT=$(check_snp "9:137684151")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    COL_RISK=$((COL_RISK + 1))
fi

if [ $COL_RISK -ge 1 ]; then
    INJURY_RISK="elevated"
    echo "  Kolagen: Warianty ryzyka kontuzji wykryte"
    TRAINING_RECS+=("ROZGRZEWKA: Minimum 15 minut przed kazda sesja! Twoje sciegna sa bardziej podatne na urazy.")
    TRAINING_RECS+=("MOBILNOSC: Codzienna praca nad mobilnoscia i elastycznoscia (yoga, stretching).")
    TRAINING_RECS+=("PROGRESJA: Wolniejsza progresja obciążeń (max +5% tygodniowo).")
    SUPPLEMENT_RECS+=("KOLAGEN: 10-15g hydrolizowanego kolagenu + 50mg wit. C - 30 min przed treningiem")
    SUPPLEMENT_RECS+=("GLICYNA: 3-5g/dzien - budulec kolagenu")
fi

# COMT - odpornosc na stres
echo "Analizuje COMT (stres i dopamina)..."
GT=$(check_snp "22:19951271")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            STRESS_TYPE="warrior"
            echo "  COMT: Val/Val → 'Warrior' - swietny pod presja"
            TRAINING_RECS+=("ZAWODY: COMT Warrior = swietna wydajnosc pod presja. Startuj w zawodach!")
            LIFESTYLE_RECS+=("STRES: Dobrze radzisz sobie ze stresem, ale mozesz potrzebowac wiecej stymulacji.")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            STRESS_TYPE="balanced"
            echo "  COMT: Val/Met → Zrownowazony profil stresowy"
            ;;
        "1/1"|"1|1")
            STRESS_TYPE="worrier"
            echo "  COMT: Met/Met → 'Worrier' - wrazliwy na stres"
            TRAINING_RECS+=("PRZYGOTOWANIE: COMT Worrier = potrzebujesz wiecej przygotowania mentalnego przed zawodami.")
            SUPPLEMENT_RECS+=("MAGNEZ: 400-600mg glycynatu magnezu wieczorem - kluczowy przy Met/Met")
            SUPPLEMENT_RECS+=("L-TEANINA: 200mg przed stresujacymi sytuacjami - wspiera GABA")
            SUPPLEMENT_RECS+=("ASHWAGANDHA: 300-600mg KSM-66 - obniza kortyzol, szczegolnie wazne dla Ciebie")
            LIFESTYLE_RECS+=("MEDYTACJA: Codzienne 10-20 min medytacji/oddychania - kluczowe przy Twoim COMT")
            AVOID_RECS+=("Ogranicz kofeine - Met/Met = wolniejszy klirens katecholamin")
            ;;
    esac
fi

# BDNF - plastycznosc mozgu i uczenie sie ruchowe
echo "Analizuje BDNF (plastycznosc mozgu)..."
GT=$(check_snp "11:27679916")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  BDNF: Val/Val → Szybkie uczenie sie nowych ruchow"
            TRAINING_RECS+=("TECHNIKA: BDNF Val/Val = szybko uczysz sie nowych ruchow. Eksperymentuj z roznymi cwiczeniami.")
            ;;
        "0/1"|"1/1")
            echo "  BDNF: Met carrier → Lepsza pamiec dlugoterminowa, ale wolniejsze uczenie"
            TRAINING_RECS+=("POWTORZENIA: BDNF Met = potrzebujesz wiecej powtorzen do opanowania techniki. Cierpliwosc!")
            TRAINING_RECS+=("CARDIO: Regularne cardio KRYTYCZNE - podnosi BDNF i kompensuje wariant Met")
            ;;
    esac
fi

echo ""

#===============================================================================
# SEKCJA 2: ANALIZA METABOLIZMU I ZYWIENIA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🥗 ANALIZA METABOLIZMU I ZYWIENIA                                        ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# FTO - otylosc
echo "Analizuje FTO (metabolizm energetyczny)..."
GT=$(check_snp "16:53820527")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  FTO: T/T → Niskie ryzyko otylosci genetycznej"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  FTO: T/A → Srednie ryzyko (+1.5 kg sredniej masy)"
            NUTRITION_RECS+=("KALORIE: FTO heterozygota - monitoruj kalorie, latwo przytyjesz.")
            ;;
        "1/1"|"1|1")
            echo "  FTO: A/A → Podwyzszone ryzyko otylosci genetycznej"
            NUTRITION_RECS+=("KALORIE: FTO A/A = genetyczna tendencja do tycia. Scisle monitoruj kalorie!")
            NUTRITION_RECS+=("BIALKO: Wyzsze bialko (2g/kg) pomoze z sytoscia przy FTO A/A")
            TRAINING_RECS+=("CWICZENIA NIWELUJA FTO: Regularna aktywnosc (5x/tydzien) calkowicie znosi efekt FTO!")
            ;;
    esac
fi

# TCF7L2 - weglowodany i cukrzyca
echo "Analizuje TCF7L2 (metabolizm weglowodanow)..."
GT=$(check_snp "10:114758349")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    CARB_SENSITIVITY="high"
    echo "  TCF7L2: Wariant ryzyka → Gorsza tolerancja weglowodanow"
    NUTRITION_RECS+=("WEGLOWODANY: TCF7L2 wariant = ogranicz weglowodany proste. Max 150-200g/dzien.")
    NUTRITION_RECS+=("TIMING: Jedz wegle glownie wokol treningu (przed i po).")
    NUTRITION_RECS+=("GI: Wybieraj zrodla niskoglikemiczne (bataty, ryz brazowy, owsianka).")
    AVOID_RECS+=("Unikaj cukrow prostych, bialego pieczywa, slodyczy - Twoj TCF7L2 to nie toleruje")
fi

# APOA2 - tluszcze nasycone
echo "Analizuje APOA2 (metabolizm tluszczow)..."
GT=$(check_snp "1:161222292")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    FAT_SENSITIVITY="high"
    echo "  APOA2: Wariant → Wrazliwosc na tluszcze nasycone"
    NUTRITION_RECS+=("TLUSZCZE NASYCONE: APOA2 wariant = ogranicz do <20g/dzien.")
    NUTRITION_RECS+=("ZAMIANY: Maslo → oliwa, ser → awokado, czerwone mieso → ryby/drob.")
    AVOID_RECS+=("Ogranicz: maslo, smalec, tlusty ser, czerwone mieso - Twoj APOA2 zle to metabolizuje")
fi

# CYP1A2 - kofeina
echo "Analizuje CYP1A2 (metabolizm kofeiny)..."
GT=$(check_snp "15:75041917")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            CAFFEINE_METABOLISM="fast"
            echo "  CYP1A2: A/A → Szybki metabolizer kofeiny"
            NUTRITION_RECS+=("KOFEINA: CYP1A2 szybki = kawa moze byc korzystna! 3-4 filizanki OK.")
            TRAINING_RECS+=("PRE-WORKOUT: 200-400mg kofeiny 30-60 min przed treningiem - dziala swietnie u Ciebie.")
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            CAFFEINE_METABOLISM="medium"
            echo "  CYP1A2: A/C → Sredni metabolizer"
            NUTRITION_RECS+=("KOFEINA: Max 2-3 kawy dziennie, ostatnia przed 14:00.")
            ;;
        "1/1"|"1|1")
            CAFFEINE_METABOLISM="slow"
            echo "  CYP1A2: C/C → Wolny metabolizer kofeiny!"
            NUTRITION_RECS+=("KOFEINA: CYP1A2 wolny = OGRANICZ kawe! Max 1-2 dziennie.")
            AVOID_RECS+=("Unikaj kofeiny po 12:00 - zaburzy sen przy Twoim CYP1A2")
            AVOID_RECS+=("Pre-workout z kofeina moze powodowac palpitacje i niepokój")
            TRAINING_RECS+=("PRE-WORKOUT: Uzyj wariantu BEZ kofeiny lub z niska dawka (50-100mg max).")
            ;;
    esac
fi

# MTHFR - metabolizm folianu
echo "Analizuje MTHFR (metabolizm folianu i metylacja)..."
GT=$(check_snp "1:11856378")
MTHFR_STATUS=""
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            MTHFR_STATUS="normal"
            echo "  MTHFR C677T: C/C → Normalny metabolizm folianu"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            MTHFR_STATUS="hetero"
            echo "  MTHFR C677T: C/T → Obnizony metabolizm (65% aktywnosci)"
            SUPPLEMENT_RECS+=("FOLIAN: L-metylofolian 400-800mcg/dzien (NIE kwas foliowy!)")
            SUPPLEMENT_RECS+=("B12: Metylokobalamina 1000mcg/dzien")
            ;;
        "1/1"|"1|1")
            MTHFR_STATUS="homo"
            echo "  MTHFR C677T: T/T → Znacznie obnizony metabolizm (30% aktywnosci)!"
            SUPPLEMENT_RECS+=("FOLIAN: L-metylofolian 800-1000mcg/dzien - KRYTYCZNE przy T/T!")
            SUPPLEMENT_RECS+=("B12: Metylokobalamina 1000-2000mcg/dzien")
            SUPPLEMENT_RECS+=("B6: P-5-P (aktywna forma) 25-50mg/dzien")
            SUPPLEMENT_RECS+=("BETAINA (TMG): 500-1000mg/dzien - wspiera alternatywna sciezke metylacji")
            AVOID_RECS+=("UNIKAJ kwasu foliowego (fortyfikowana zywnosc, tanie suplementy) - nie przetworzysz go!")
            NUTRITION_RECS+=("ZIELONE LISCIASTE: Duzo szpinaku, jarmuzu - naturalne zrodlo folianu")
            ;;
    esac
fi

# VDR - witamina D
echo "Analizuje VDR (receptor witaminy D)..."
GT=$(check_snp "12:48272895")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  VDR: Wariant → Wieksze zapotrzebowanie na wit. D"
    SUPPLEMENT_RECS+=("WITAMINA D3: 4000-5000 IU/dzien (wiecej niz standardowe 1000-2000)")
    SUPPLEMENT_RECS+=("WITAMINA K2: 100-200mcg MK-7 - zawsze z D3 dla prawidlowego metabolizmu wapnia")
fi

# FADS1 - omega-3
echo "Analizuje FADS1 (konwersja omega-3)..."
GT=$(check_snp "11:61567029")
if [ -n "$GT" ]; then
    if [[ "$GT" == "1/1" ]] || [[ "$GT" == "1|1" ]]; then
        echo "  FADS1: T/T → Slaba konwersja ALA do EPA/DHA"
        SUPPLEMENT_RECS+=("OMEGA-3: EPA+DHA bezposrednio 2-3g/dzien (nie z lnu/chia - nie przekonwertujesz)")
        NUTRITION_RECS+=("RYBY: Jedz tłuste ryby 3-4x/tydzien (losos, makrela, sardynki)")
        AVOID_RECS+=("Nie polegaj na roslinnych omega-3 (len, chia) - Twoj FADS1 nie przekonwertuje ich!")
    fi
fi

# BCMO1 - witamina A
echo "Analizuje BCMO1 (konwersja beta-karotenu)..."
GT=$(check_snp "16:81264597")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  BCMO1: Wariant → Slaba konwersja beta-karotenu do wit. A"
    SUPPLEMENT_RECS+=("WITAMINA A: Retinol 2500-5000 IU/dzien LUB jedz watrobke 1x/tydzien")
    NUTRITION_RECS+=("RETINOL: Jedz zrodla gotowej wit. A (watrobka, jaja, maslo) nie marchewke")
fi

# FUT2 - witamina B12
echo "Analizuje FUT2 (wchłanianie B12)..."
GT=$(check_snp "19:49206674")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  FUT2: Non-secretor → Nizsze wchlanianie B12"
    SUPPLEMENT_RECS+=("B12: Wyzsza dawka metylokobalaminy 1000-2000mcg/dzien lub zastrzyki")
fi

# ADH1B/ALDH2 - alkohol
echo "Analizuje metabolizm alkoholu..."
GT=$(check_snp "12:112241766")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  ALDH2: Wariant → Asian flush, nietolerancja alkoholu!"
    AVOID_RECS+=("ALKOHOL: Calkowicie unikaj! ALDH2 wariant = toksyczne aldehydy + ryzyko raka przeluku")
fi

GT=$(check_snp "4:99318162")
if [ -n "$GT" ] && [[ "$GT" == "0/0" ]]; then
    echo "  ADH1B: Wolny metabolizm alkoholu"
    AVOID_RECS+=("ALKOHOL: Ogranicz - wolny metabolizm = wieksze ryzyko uzaleznienia i uszkodzen")
fi

# Laktozaecho "Analizuje LCT (tolerancja laktozy)..."
GT=$(check_snp "2:136608646")
if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  LCT: C/C → Nietolerancja laktozy"
            AVOID_RECS+=("LAKTOZA: Unikaj mleka i produktow mlecznych lub uzywaj laktazy")
            NUTRITION_RECS+=("NABIAŁ: Wybieraj produkty bezlaktozowe lub fermentowane (jogurt, kefir, sery dojrzale)")
            ;;
        *)
            echo "  LCT: Tolerancja laktozy OK"
            ;;
    esac
fi

echo ""

#===============================================================================
# SEKCJA 3: HORMONY I UKŁAD HORMONALNY
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  🧪 ANALIZA UKLADU HORMONALNEGO                                           ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "UWAGA: To NIE sa rekomendacje medyczne. Hormonoterapia wymaga nadzoru lekarza!"
echo ""

# SHBG - globulina wiazaca hormony
echo "Analizuje SHBG (transport hormonow)..."
GT=$(check_snp "17:7534706")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  SHBG: Wariant → Moze wplywac na poziom wolnego testosteronu"
    HORMONE_NOTES+=("SHBG: Wariant moze oznaczac wyzszy SHBG = mniej wolnego testosteronu. Badaj FT nie tylko TT.")
fi

# AR - receptor androgenowy (CAG repeats - trudne do analizy z SNP)
echo "  AR (receptor androgenowy): Wymaga analizy STR (powtorzen CAG)"
HORMONE_NOTES+=("AR: Dlugosc CAG wplywa na wrazliwosc na androgeny. Krotsze = wieksza wrazliwosc.")

# CYP19A1 - aromataza
echo "Analizuje CYP19A1 (aromataza)..."
GT=$(check_snp "15:51502956")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  CYP19A1: Wariant → Moze wplywac na konwersje T→E2"
    HORMONE_NOTES+=("AROMATAZA: Wariant CYP19A1 - monitoruj stosunek T/E2. Moze wymagac modulacji.")
    NUTRITION_RECS+=("KRUCYFEROWE: Brokuly, kalafior, kapusta - naturalne modulatory aromatazy (DIM, I3C)")
fi

# SRD5A2 - 5-alfa reduktaza
echo "Analizuje SRD5A2 (5-alfa reduktaza)..."
GT=$(check_snp "2:31749strona")
# Pozycja przykladowa
HORMONE_NOTES+=("5AR: Aktywnosc 5-alfa reduktazy wplywa na konwersje T→DHT (lysienie, prostata, sila).")

# CYP17A1 - synteza androgenow
echo "Analizuje CYP17A1 (synteza androgenow)..."
GT=$(check_snp "10:102830531")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    HORMONE_NOTES+=("CYP17A1: Wariant moze wplywac na synteze androgenow nadnerczowych (DHEA, androstendion).")
fi

echo ""
echo "Notatki hormonalne:"
for note in "${HORMONE_NOTES[@]}"; do
    echo "  • $note"
done

echo ""

#===============================================================================
# SEKCJA 4: SEN I REGENERACJA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  😴 SEN I REGENERACJA                                                     ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Chronotyp
MORNING_SCORE=0
EVENING_SCORE=0

GT=$(check_snp "2:239186871")  # PER2
if [[ "$GT" == "1/1" ]]; then
    MORNING_SCORE=$((MORNING_SCORE + 2))
fi

GT=$(check_snp "4:56294068")  # CLOCK
if [[ "$GT" != "0/0" ]] && [ -n "$GT" ]; then
    EVENING_SCORE=$((EVENING_SCORE + 2))
fi

if [ $MORNING_SCORE -gt $EVENING_SCORE ]; then
    echo "  Chronotyp: 🐦 SKOWRONEK (ranny ptaszek)"
    LIFESTYLE_RECS+=("SEN: Chodz spac 21:00-22:00, wstawaj 5:00-6:00. Trenuj rano (6-10).")
    TRAINING_RECS+=("TIMING: Najlepsza wydajnosc treningowa: 8:00-12:00")
elif [ $EVENING_SCORE -gt $MORNING_SCORE ]; then
    echo "  Chronotyp: 🦉 SOWA (nocny marek)"
    LIFESTYLE_RECS+=("SEN: Chodz spac 23:00-00:00, wstawaj 7:00-8:00. Trenuj popoludniu/wieczorem.")
    TRAINING_RECS+=("TIMING: Najlepsza wydajnosc treningowa: 16:00-20:00")
else
    echo "  Chronotyp: 🕐 Posredni"
    LIFESTYLE_RECS+=("SEN: Elastyczny rytm - dostosuj do stylu zycia, ale utrzymuj stalą porę.")
fi

# ADA - sen głęboki
GT=$(check_snp "20:44619933")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]]; then
    echo "  ADA: Wariant → Wiecej snu glebokiego"
    LIFESTYLE_RECS+=("SEN: Twoj wariant ADA = potrzebujesz wiecej snu (8-9h). Nie skracaj!")
fi

echo ""

#===============================================================================
# SEKCJA 5: PODSUMOWANIE PROTOKOLU
#===============================================================================

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════════════════╗"
echo "║                                                                                       ║"
echo "║                    TWOJ SPERSONALIZOWANY PROTOKOL                                     ║"
echo "║                                                                                       ║"
echo "╚═══════════════════════════════════════════════════════════════════════════════════════╝"
echo ""

# PROFIL
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📊 TWOJ PROFIL GENETYCZNY:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "   Typ atletyczny:      $([ $POWER_SCORE -gt $ENDURANCE_SCORE ] && echo "💪 SILOWY" || ([ $ENDURANCE_SCORE -gt $POWER_SCORE ] && echo "🏃 WYTRZYMALOSCIOWY" || echo "⚖️ MIESZANY"))"
echo "   Regeneracja:         $([ "$RECOVERY_SPEED" == "fast" ] && echo "⚡ Szybka" || ([ "$RECOVERY_SPEED" == "slow" ] && echo "🐢 Wolna" || echo "🔄 Normalna"))"
echo "   Ryzyko kontuzji:     $([ "$INJURY_RISK" == "elevated" ] && echo "⚠️ Podwyzszone" || echo "✅ Normalne")"
echo "   Typ stresowy:        $([ "$STRESS_TYPE" == "warrior" ] && echo "⚔️ Warrior" || ([ "$STRESS_TYPE" == "worrier" ] && echo "🧘 Worrier" || echo "⚖️ Balanced"))"
echo "   Metabolizm kofeiny:  $([ "$CAFFEINE_METABOLISM" == "fast" ] && echo "☕ Szybki" || ([ "$CAFFEINE_METABOLISM" == "slow" ] && echo "🐢 Wolny" || echo "🔄 Sredni"))"
echo "   MTHFR:               $([ "$MTHFR_STATUS" == "homo" ] && echo "🔴 T/T (wymaga uwagi!)" || ([ "$MTHFR_STATUS" == "hetero" ] && echo "🟡 C/T" || echo "✅ C/C"))"
echo ""

# TRENING
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🏋️ REKOMENDACJE TRENINGOWE:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
for rec in "${TRAINING_RECS[@]}"; do
    echo "   ✓ $rec"
    echo ""
done

# SUPLEMENTACJA
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "💊 REKOMENDOWANA SUPLEMENTACJA:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Podstawowa suplementacja dla kazdego
echo "   BAZA (dla kazdego):"
echo "   ├─ Witamina D3: 2000-4000 IU/dzien"
echo "   ├─ Omega-3: 2-3g EPA+DHA/dzien"
echo "   ├─ Magnez: 300-400mg wieczorem"
echo "   └─ Cynk: 15-30mg/dzien (jesli trenujesz ciezko)"
echo ""
echo "   SPERSONALIZOWANE (na podstawie Twoich genow):"
for rec in "${SUPPLEMENT_RECS[@]}"; do
    echo "   ✓ $rec"
done
echo ""

# ZYWIENIE
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🥗 REKOMENDACJE ZYWIENIOWE:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
for rec in "${NUTRITION_RECS[@]}"; do
    echo "   ✓ $rec"
done
echo ""

# CZEGO UNIKAC
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "⛔ CZEGO UNIKAC:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
for rec in "${AVOID_RECS[@]}"; do
    echo "   ✗ $rec"
done
echo ""

# STYL ZYCIA
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🌟 STYL ZYCIA:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
for rec in "${LIFESTYLE_RECS[@]}"; do
    echo "   ✓ $rec"
done
echo ""

#===============================================================================
# ZAPIS DO PLIKU
#===============================================================================

REPORT="$OUTDIR/PROTOCOL_$(date +%Y%m%d).txt"

{
echo "================================================================================"
echo "          SPERSONALIZOWANY PROTOKOL GENETYCZNY"
echo "          Wygenerowano: $(date)"
echo "================================================================================"
echo ""
echo "PROFIL:"
echo "  Typ atletyczny: $([ $POWER_SCORE -gt $ENDURANCE_SCORE ] && echo "SILOWY" || ([ $ENDURANCE_SCORE -gt $POWER_SCORE ] && echo "WYTRZYMALOSCIOWY" || echo "MIESZANY"))"
echo "  ACTN3: $ACTN3_STATUS"
echo "  Regeneracja: $RECOVERY_SPEED"
echo "  COMT: $STRESS_TYPE"
echo "  CYP1A2 (kofeina): $CAFFEINE_METABOLISM"
echo "  MTHFR: $MTHFR_STATUS"
echo ""
echo "TRENING:"
for rec in "${TRAINING_RECS[@]}"; do
    echo "  - $rec"
done
echo ""
echo "SUPLEMENTACJA:"
for rec in "${SUPPLEMENT_RECS[@]}"; do
    echo "  - $rec"
done
echo ""
echo "ZYWIENIE:"
for rec in "${NUTRITION_RECS[@]}"; do
    echo "  - $rec"
done
echo ""
echo "UNIKAC:"
for rec in "${AVOID_RECS[@]}"; do
    echo "  - $rec"
done
echo ""
echo "================================================================================"
echo "UWAGA: To nie sa rekomendacje medyczne. Skonsultuj z lekarzem przed zmianami."
echo "================================================================================"
} > "$REPORT"

echo ""
echo "📁 Pelny raport zapisany: $REPORT"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "⚠️  DISCLAIMER: To są rekomendacje informacyjne na podstawie genetyki."
echo "    Nie zastępują porady lekarza, dietetyka ani trenera."
echo "    Hormonoterapia wymaga ZAWSZE nadzoru specjalisty!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
