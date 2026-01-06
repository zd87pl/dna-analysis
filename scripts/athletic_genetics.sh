#!/bin/bash
#===============================================================================
# ROZSZERZONY PROFIL GENETYKI SPORTOWEJ
# Kompleksowa analiza potencjalu atletycznego
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║           ROZSZERZONY PROFIL GENETYKI SPORTOWEJ                           ║"
echo "║           Kompleksowa analiza potencjalu atletycznego                     ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

#===============================================================================
# KATEGORIA 1: SILA I MOC
#===============================================================================

echo "💪 SILA I MOC (Power)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

POWER_SCORE=0

# ACTN3 - kluczowy gen sily
result=$(bcftools query -r 11:66560624 -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1)
GT=$(echo $result | awk '{print $3}')

if [ -n "$GT" ]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ACTN3 (rs1815739): R/R → 🏆 ELITARNA genetyka silowa!"
            echo "    • Pelna funkcja alfa-aktyniny-3 w wloknach szybkokurczliwych"
            echo "    • Przewaga w: sprint, skoki, rzuty, podnoszenie ciezarow"
            echo "    • ~18% populacji ma ten genotyp"
            POWER_SCORE=$((POWER_SCORE + 3))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  ACTN3 (rs1815739): R/X → Dobra mieszanka sila/wytrzymalosc"
            echo "    • Dobry dla sportow mieszanych (pilka nozna, koszykowka)"
            POWER_SCORE=$((POWER_SCORE + 2))
            ;;
        "1/1"|"1|1")
            echo "  ACTN3 (rs1815739): X/X → Genetyka wytrzymalosciowa"
            echo "    • Lepsza wydajnosc w sportach wytrzymalosciowych"
            echo "    • ~18% populacji ma ten genotyp"
            POWER_SCORE=$((POWER_SCORE + 0))
            ;;
    esac
fi

# AGT - angiotensynogen (sila miesni)
result=$(bcftools query -r 1:230845794 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "1/1"|"1|1")
            echo "  AGT (rs699): T/T → Wyzsza sila miesniowa"
            POWER_SCORE=$((POWER_SCORE + 1))
            ;;
        "0/1")
            echo "  AGT (rs699): M/T → Srednia predyspozycja"
            ;;
    esac
fi

# IL6 - odpowiedz na trening silowy
result=$(bcftools query -r 7:22727026 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  IL6 (rs1800795): G/G → Lepsza odpowiedz na trening silowy"
            POWER_SCORE=$((POWER_SCORE + 1))
            ;;
    esac
fi

# MSTN - miostatyna (hamulec wzrostu miesni)
result=$(bcftools query -r 2:190379711 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  MSTN: Wariant → Potencjalnie lepsza hipertrofia miesni"
    POWER_SCORE=$((POWER_SCORE + 1))
fi

echo ""
echo "  WYNIK SILY: $POWER_SCORE / 6"
echo ""

#===============================================================================
# KATEGORIA 2: WYTRZYMALOSC
#===============================================================================

echo "🏃 WYTRZYMALOSC (Endurance)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

ENDURANCE_SCORE=0

# ACE - enzym konwertujacy angiotensyne
result=$(bcftools query -r 17:63488529 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "1/1"|"1|1")
            echo "  ACE (I/D): I/I → Lepsza wytrzymalosc, wydajnosc tlenowa"
            echo "    • Przewaga w: biegi dlugie, kolarstwo, plywanie"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 2))
            ;;
        "0/1")
            echo "  ACE (I/D): I/D → Zrownowazony profil"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
        "0/0"|"0|0")
            echo "  ACE (I/D): D/D → Lepsza sila i moc"
            ;;
    esac
fi

# PPARGC1A (PGC-1alpha) - biogeneza mitochondriow
result=$(bcftools query -r 4:23814039 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  PPARGC1A (rs8192678): G/G → Lepsza adaptacja wytrzymalosciowa"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 2))
            ;;
        "0/1")
            echo "  PPARGC1A (rs8192678): G/A → Srednia odpowiedz"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
    esac
fi

# VEGFA - angiogeneza
result=$(bcftools query -r 6:43770613 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  VEGFA (rs2010963): G/G → Lepsza angiogeneza"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
    esac
fi

# NRF2 - odpowiedz antyoksydacyjna
result=$(bcftools query -r 2:178095031 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  NFE2L2: Wariant → Lepsza ochrona antyoksydacyjna"
    ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
fi

# PPARA - metabolizm tluszczow jako paliwa
result=$(bcftools query -r 22:46615880 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  PPARA (rs4253778): G/G → Lepsze spalanie tluszczow"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
    esac
fi

echo ""
echo "  WYNIK WYTRZYMALOSCI: $ENDURANCE_SCORE / 7"
echo ""

#===============================================================================
# KATEGORIA 3: REGENERACJA
#===============================================================================

echo "🔄 REGENERACJA (Recovery)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

RECOVERY_SCORE=0

# IL6 - zapalenie po treningu
result=$(bcftools query -r 7:22727026 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  IL6 (rs1800795): G/G → Szybsza regeneracja, mniej zapalenia"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
            ;;
        "0/1")
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
        "1/1"|"1|1")
            echo "  IL6 (rs1800795): C/C → Wolniejsza regeneracja, wiecej DOMS"
            ;;
    esac
fi

# CRP - C-reaktywne bialko
result=$(bcftools query -r 1:159682233 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ] && [ "$result" == "0/0" ]; then
    echo "  CRP (rs1205): Nizszy poziom zapalenia"
    RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
fi

# TNF - czynnik martwicy nowotworu
result=$(bcftools query -r 6:31575254 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  TNF (rs1800629): G/G → Nizsza odpowiedz zapalna"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
    esac
fi

# SOD2 - dysmutaza ponadtlenkowa
result=$(bcftools query -r 6:160113872 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/1")
            echo "  SOD2 (rs4880): Heterozygota → Optymalna ochrona antyoksydacyjna"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
            ;;
        "0/0"|"1/1")
            echo "  SOD2 (rs4880): Homozygota → Mniej optymalna ochrona"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
    esac
fi

# COL1A1 - kolagen (regeneracja sciegien)
result=$(bcftools query -r 17:50201632 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  COL1A1 (rs1800012): G/G → Lepszy kolagen, mniej kontuzji"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
    esac
fi

echo ""
echo "  WYNIK REGENERACJI: $RECOVERY_SCORE / 7"
echo ""

#===============================================================================
# KATEGORIA 4: RYZYKO KONTUZJI
#===============================================================================

echo "🩹 RYZYKO KONTUZJI"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

INJURY_RISK=0

# COL5A1 - ACL i urazy sciegien
result=$(bcftools query -r 9:137684151 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "1/1"|"1|1")
            echo "  COL5A1 (rs12722): T/T → Podwyzszone ryzyko urazow ACL"
            INJURY_RISK=$((INJURY_RISK + 2))
            ;;
        "0/1")
            echo "  COL5A1 (rs12722): C/T → Umiarkowane ryzyko urazow"
            INJURY_RISK=$((INJURY_RISK + 1))
            ;;
        "0/0"|"0|0")
            echo "  COL5A1 (rs12722): C/C → Nizsze ryzyko urazow sciegien"
            ;;
    esac
fi

# GDF5 - rozwoj chrzastki
result=$(bcftools query -r 20:34025756 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "1/1"|"1|1")
            echo "  GDF5 (rs143383): A/A → Podwyzszone ryzyko osteoarthritis"
            INJURY_RISK=$((INJURY_RISK + 1))
            ;;
    esac
fi

# MMP3 - degradacja macierzy
result=$(bcftools query -r 11:102731092 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  MMP3: Wariant → Wieksze ryzyko urazow sciegna Achillesa"
    INJURY_RISK=$((INJURY_RISK + 1))
fi

echo ""
if [ $INJURY_RISK -le 1 ]; then
    echo "  RYZYKO KONTUZJI: NISKIE ✅"
elif [ $INJURY_RISK -le 2 ]; then
    echo "  RYZYKO KONTUZJI: SREDNIE 🟡"
else
    echo "  RYZYKO KONTUZJI: PODWYZSZONE 🔴"
    echo "  → Rekomendacja: Dluzsze rozgrzewki, praca nad elastycznoscia"
fi
echo ""

#===============================================================================
# KATEGORIA 5: PSYCHOLOGIA SPORTOWA
#===============================================================================

echo "🧠 PSYCHOLOGIA SPORTOWA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# COMT - odpornosc na stres rywalizacji
result=$(bcftools query -r 22:19951271 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  COMT (rs4680): Val/Val → 'Warrior'"
            echo "    • Lepsza wydajnosc pod presja"
            echo "    • Nizszy poziom dopaminy bazowo = lepszy podczas stresu"
            echo "    • Dobry dla: zawody, mecze, wazne starty"
            ;;
        "0/1"|"0|1")
            echo "  COMT (rs4680): Val/Met → Zrownowazony profil"
            echo "    • Dobra elastycznosc poznawcza"
            ;;
        "1/1"|"1|1")
            echo "  COMT (rs4680): Met/Met → 'Worrier'"
            echo "    • Lepsza koncentracja i precyzja w spokojnych warunkach"
            echo "    • Moze potrzebowac technik radzenia ze stresem przed zawodami"
            echo "    • Dobry dla: trening techniczny, sporty precyzyjne"
            ;;
    esac
fi

# DRD4 - poszukiwanie wrazn / ryzyko
result=$(bcftools query -r 11:637269 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ] && [ "$result" != "0/0" ]; then
    echo "  DRD4: Wariant 7R → Wieksze poszukiwanie wrazn"
    echo "    • Dobry dla: sporty ekstremalne, wysokie ryzyko"
fi

# BDNF - uczenie sie ruchowe
result=$(bcftools query -r 11:27679916 -f '[%GT]\n' "$VCF" 2>/dev/null | head -1)
if [ -n "$result" ]; then
    case $result in
        "0/0"|"0|0")
            echo "  BDNF (rs6265): Val/Val → Szybsze uczenie sie nowych ruchow"
            ;;
        "0/1"|"1/1")
            echo "  BDNF (rs6265): Met carrier → Lepsza pamiec ruchowa dlugoterminowa"
            ;;
    esac
fi

echo ""

#===============================================================================
# PODSUMOWANIE PROFILU
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    PODSUMOWANIE PROFILU SPORTOWEGO                        ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Oblicz calkowite wyniki
POWER_PCT=$(echo "scale=0; $POWER_SCORE * 100 / 6" | bc)
ENDURANCE_PCT=$(echo "scale=0; $ENDURANCE_SCORE * 100 / 7" | bc)
RECOVERY_PCT=$(echo "scale=0; $RECOVERY_SCORE * 100 / 7" | bc)

echo "WYNIKI:"
echo ""
printf "  Sila i Moc:      "
for i in $(seq 1 10); do
    if [ $i -le $(($POWER_PCT / 10)) ]; then
        printf "█"
    else
        printf "░"
    fi
done
echo " $POWER_PCT%"

printf "  Wytrzymalosc:    "
for i in $(seq 1 10); do
    if [ $i -le $(($ENDURANCE_PCT / 10)) ]; then
        printf "█"
    else
        printf "░"
    fi
done
echo " $ENDURANCE_PCT%"

printf "  Regeneracja:     "
for i in $(seq 1 10); do
    if [ $i -le $(($RECOVERY_PCT / 10)) ]; then
        printf "█"
    else
        printf "░"
    fi
done
echo " $RECOVERY_PCT%"

echo ""

# Rekomendacje sportow
echo "REKOMENDOWANE DYSCYPLINY:"
echo ""

if [ $POWER_SCORE -ge 4 ] && [ $ENDURANCE_SCORE -le 3 ]; then
    echo "  🥇 SPORTY SILOWE I SZYBKOSCIOWE:"
    echo "     • Sprint (100m, 200m)"
    echo "     • Podnoszenie ciezarow"
    echo "     • Skoki"
    echo "     • Rzuty"
    echo "     • CrossFit (komponenty silowe)"
    echo "     • Futbol amerykanski"
elif [ $ENDURANCE_SCORE -ge 4 ] && [ $POWER_SCORE -le 3 ]; then
    echo "  🥇 SPORTY WYTRZYMALOSCIOWE:"
    echo "     • Biegi dlugie (5km+)"
    echo "     • Kolarstwo szosowe"
    echo "     • Plywanie dlugodystansowe"
    echo "     • Triathlon"
    echo "     • Biegi gorskie"
elif [ $POWER_SCORE -ge 3 ] && [ $ENDURANCE_SCORE -ge 3 ]; then
    echo "  🥇 SPORTY MIESZANE (SILA + WYTRZYMALOSC):"
    echo "     • Pilka nozna"
    echo "     • Koszykowka"
    echo "     • Tenis"
    echo "     • MMA / Judo / Zapasy"
    echo "     • CrossFit"
    echo "     • Biegi sredniodystansowe (800m-1500m)"
fi

echo ""
echo "WSKAZOWKI TRENINGOWE:"
echo ""

if [ $RECOVERY_SCORE -le 3 ]; then
    echo "  ⚠️ Regeneracja: Wolniejsza - zaplanuj wiecej dni odpoczynku"
    echo "     • 48-72h miedzy intensywnymi treningami"
    echo "     • Techniki regeneracji: sen, masaz, sauna"
fi

if [ $INJURY_RISK -ge 2 ]; then
    echo "  ⚠️ Kontuzje: Podwyzszone ryzyko - ostroznosc!"
    echo "     • Dluzsze rozgrzewki (15-20 min)"
    echo "     • Praca nad elastycznoscia i stabilizacja"
    echo "     • Unikaj gwaltownych zmian obciazenia"
fi

echo ""
