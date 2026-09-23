#!/bin/bash
#===============================================================================
# HELIXIGHT LOCAL ANALYSIS TOOLKIT
# Kompleksowa analiza VCF na macOS
#===============================================================================

# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"
set -eo pipefail

VCF="${1:-}"
hx_require_arg "$VCF" "$(basename "$0") <variants.vcf.gz> [reads.bam] [output_dir]"
BAM="${2:-}"   # opcjonalny: używany tylko w podpowiedzi dla mtDNA
OUTPUT_DIR="${3:-./helixight_analysis}"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║          HELIXIGHT LOCAL ANALYSIS TOOLKIT                         ║"
echo "║          Analiza VCF + BAM na macOS                               ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "VCF: $VCF"
echo "BAM: ${BAM:-<brak>}"
echo "Output: $OUTPUT_DIR"
echo ""

# Współrzędne GRCh38 (scripts/lib/GRCh38_COORDINATES.md). Biblioteka normalizuje
# nazwy kontigów i buduje zindeksowaną kopię roboczą; dalej używamy tylko jej.
hx_prepare_vcf "$VCF"
VCF="$HX_VCF"

# snp_info REGION [FALLBACK]: "REF>ALT GT" pierwszego SNP w regionie (genotyp
# z fazą znormalizowaną); FALLBACK (domyślnie N/A) gdy brak rekordu / wywołania.
snp_info() {
    local gt alleles
    gt=$(hx_gt "$1")
    if [[ -z $gt ]]; then
        printf '%s\n' "${2:-N/A}"
        return 0
    fi
    alleles=$(hx_query "$1" '%REF>%ALT\n' -e 'GT="mis" || TYPE!="snp"')
    printf '%s %s\n' "${alleles%%$'\n'*}" "$gt"
}

mkdir -p "$OUTPUT_DIR"/{pharmacogenomics,ancestry,traits,health,reports}

#===============================================================================
# MENU ANALIZ
#===============================================================================

echo "Dostępne analizy:"
echo ""
echo "  1) Rozszerzona farmakogenomika (CYP haplotypy)"
echo "  2) Ancestry & pochodzenie (haplogrupy)"
echo "  3) Cechy fizyczne (kolor oczu, włosów, etc.)"
echo "  4) Ryzyko chorób (rozszerzone SNP)"
echo "  5) mtDNA - haplogrupa matczyna"
echo "  6) Y-DNA - haplogrupa ojcowska (jeśli mężczyzna)"
echo "  7) Runs of Homozygosity (ROH)"
echo "  8) Statystyki genomu"
echo "  9) WSZYSTKIE ANALIZY"
echo "  0) Wyjście"
echo ""

read -p "Wybierz analizę [0-9]: " choice

#===============================================================================
# 1. ROZSZERZONA FARMAKOGENOMIKA
#===============================================================================

run_pharmacogenomics() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  ANALIZA FARMAKOGENOMICZNA"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/pharmacogenomics/pharma_report.txt"
    
    {
        echo "RAPORT FARMAKOGENOMICZNY"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        # CYP2D6 - kluczowe warianty
        echo "=== CYP2D6 (metabolizm ~25% leków) ==="
        # GRCh38: przesunięcie GRCh37->38 w tym regionie (-396002) wynika z kanonicznych
        # pozycji rs3892097/rs1065852 w scripts/lib/GRCh38_COORDINATES.md
        echo "Lokalizacja: chr22:42126498-42130881 (GRCh38)"
        CYP2D6_REGION=$(hx_query 22:42126498-42130881 '  %CHROM:%POS %REF>%ALT (QUAL=%QUAL)\n')  # TODO(build): unverified GRCh38
        echo "${CYP2D6_REGION:-  Brak wariantów w regionie}"
        echo ""
        
        # Kluczowe SNP CYP2D6 (GRCh38)
        echo "Kluczowe SNP CYP2D6:"
        for snp in "22:42128945 rs3892097 (*4)" "22:42130692 rs1065852 (*10)" "22:42127941 rs16947 (*2)"; do
            read -r pos label <<<"$snp"
            echo "  $label $pos: $(snp_info "$pos")"
        done
        echo ""
        
        # CYP2C19
        echo "=== CYP2C19 (klopidogrel, IPP, SSRI) ==="
        echo "rs4244285 (*2): $(snp_info 10:94781859)"
        echo "rs12248560 (*17): $(snp_info 10:94761900)"
        echo ""
        
        # CYP2C9
        echo "=== CYP2C9 (warfaryna, NLPZ) ==="
        echo "rs1799853 (*2): $(snp_info 10:94942290)"
        echo "rs1057910 (*3): $(snp_info 10:94981296)"
        echo ""
        
        # CYP3A4/CYP3A5
        echo "=== CYP3A4/CYP3A5 (statyny, immunosupresja) ==="
        echo "CYP3A5*3 rs776746: $(snp_info 7:99672916)"
        echo "CYP3A4*22 rs35599367: $(snp_info 7:99768693)"
        echo ""
        
        # VKORC1
        echo "=== VKORC1 (warfaryna) ==="
        echo "rs9923231: $(snp_info 16:31096368)"
        echo ""
        
        # SLCO1B1
        echo "=== SLCO1B1 (statyny - miopatia) ==="
        echo "rs4149056 (*5): $(snp_info 12:21178615)"
        echo ""
        
        # DPYD
        echo "=== DPYD (chemioterapia 5-FU) ==="
        echo "rs3918290 (*2A): $(snp_info 1:97915614)"
        echo "rs55886062 (*13): $(snp_info 1:98348885)"  # TODO(build): unverified GRCh38
        echo ""
        
        # TPMT
        echo "=== TPMT (tiopuryny - azatiopryna) ==="
        echo "rs1800462 (*2): $(snp_info 6:18143724)"
        echo "rs1800460 (*3B): $(snp_info 6:18138997)"
        echo "rs1142345 (*3C): $(snp_info 6:18130687)"
        echo ""
        
        # HLA - reakcje na leki
        echo "=== HLA (ciężkie reakcje skórne) ==="
        echo "Region HLA-B: chr6:31321649-31324989"
        HLAB_COUNT=$(hx_query 6:31321649-31324989 '%POS\n' | wc -l | tr -d ' ')  # TODO(build): unverified GRCh38
        echo "  Znaleziono $HLAB_COUNT wariantów w regionie HLA-B"
        echo ""
        
    } > "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# 2. ANCESTRY - POCHODZENIE
#===============================================================================

run_ancestry() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  ANALIZA POCHODZENIA (Ancestry Informative Markers)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/ancestry/ancestry_report.txt"
    
    {
        echo "RAPORT POCHODZENIA"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        # AIMs - Ancestry Informative Markers
        echo "=== MARKERY INFORMACYJNE POCHODZENIA ==="
        echo ""
        
        # Europejskie markery
        echo "Markery europejskie:"
        echo "  rs1426654 (SLC24A5 - jasna skóra): $(snp_info 15:48426484)"  # TODO(build): unverified GRCh38
        echo "  rs16891982 (SLC45A2 - jasna skóra): $(snp_info 5:33951693)"  # TODO(build): unverified GRCh38
        echo "  rs12913832 (HERC2 - niebieskie oczy): $(snp_info 15:28120472)"
        echo ""
        
        # Azjatyckie markery
        echo "Markery azjatyckie:"
        echo "  rs3827760 (EDAR - grube włosy): $(snp_info 2:109513601)"  # TODO(build): unverified GRCh38
        echo "  rs671 (ALDH2 - Asian flush): $(snp_info 12:111803962)"
        echo ""
        
        # Afrykańskie markery
        echo "Markery afrykańskie:"
        echo "  rs2814778 (DARC - malaria): $(snp_info 1:159174683)"  # TODO(build): unverified GRCh38
        echo ""
        
        # Tolerancja laktozy (różnice populacyjne)
        echo "Adaptacje populacyjne:"
        echo "  rs4988235 (LCT - tolerancja laktozy): $(snp_info 2:135851076)"
        echo "  rs182549 (LCT - alt marker): $(snp_info 2:136616754)"  # TODO(build): unverified GRCh38
        echo ""
        
    } > "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# 3. CECHY FIZYCZNE
#===============================================================================

run_traits() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  CECHY FIZYCZNE (Physical Traits)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/traits/traits_report.txt"
    
    {
        echo "RAPORT CECH FIZYCZNYCH"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        # KOLOR OCZU
        echo "=== KOLOR OCZU ==="
        echo "rs12913832 (HERC2 - główny): $(snp_info 15:28120472)"
        echo "  AA = brązowe oczy (dominujące)"
        echo "  AG = brązowe/zielone"
        echo "  GG = niebieskie oczy"
        echo ""
        echo "rs1800407 (OCA2): $(snp_info 15:28230318)"  # TODO(build): unverified GRCh38
        echo "rs12896399 (SLC24A4): $(snp_info 14:92773663)"  # TODO(build): unverified GRCh38
        echo ""
        
        # KOLOR WŁOSÓW
        echo "=== KOLOR WŁOSÓW ==="
        echo "rs12821256 (KITLG - blond): $(snp_info 12:89328335)"  # TODO(build): unverified GRCh38
        echo "rs1805007 (MC1R - rude): $(snp_info 16:89919436)"
        echo "rs1805008 (MC1R - rude): $(snp_info 16:89986091)"  # TODO(build): unverified GRCh38
        echo "rs1805009 (MC1R - rude): $(snp_info 16:89986144)"  # TODO(build): unverified GRCh38
        echo ""
        
        # ŁYSIENIE
        echo "=== ŁYSIENIE (androgenowe) ==="
        echo "rs2180439 (chr20): $(snp_info 20:21985252)"  # TODO(build): unverified GRCh38
        echo "rs6152 (AR - receptor androgenowy): $(snp_info X:66943551)"  # TODO(build): unverified GRCh38
        echo ""
        
        # PIEGI
        echo "=== PIEGI ==="
        echo "rs4911414 (IRF4): $(snp_info 6:396321)"
        echo "rs1805007 (MC1R): już powyżej"
        echo ""
        
        # GORZKI SMAK
        echo "=== PERCEPCJA SMAKU ==="
        echo "rs713598 (TAS2R38 - gorzki smak): $(snp_info 7:141972804)"
        echo "rs1726866 (TAS2R38): $(snp_info 7:141672705)"  # TODO(build): unverified GRCh38
        echo "  PAV/PAV = silnie czujesz gorzki smak"
        echo "  AVI/AVI = nie czujesz gorzkiego"
        echo ""
        
        # KICHANIE NA SŁOŃCE
        echo "=== INNE CECHY ==="
        # rs10427255: pozycja GRCh38 niezweryfikowana (~2:145.37 Mb); 2:135851076 to rs4988235 (LCT)
        echo "rs10427255 (kichanie na słońce - photic sneeze): pozycja niezweryfikowana - pominięto"
        echo "rs17822931 (ABCC11 - typ woskowiny): $(snp_info 16:48224287)"
        echo "  CC/CT = mokra woskowina (europejski typ)"
        echo "  TT = sucha woskowina (azjatycki typ)"
        echo ""
        
        # KOFEINA
        echo "=== METABOLIZM KOFEINY ==="
        echo "rs762551 (CYP1A2): $(snp_info 15:74749576)"
        echo "  AA = szybki metabolizer (kawa OK)"
        echo "  AC/CC = wolny metabolizer (ogranicz kawę)"
        echo ""
        
        # ALKOHOL
        echo "=== METABOLIZM ALKOHOLU ==="
        echo "rs1229984 (ADH1B): $(snp_info 4:99318162)"
        echo "rs671 (ALDH2 - Asian flush): $(snp_info 12:111803962)"
        echo ""
        
        # WZROST
        echo "=== WZROST (wybrane warianty) ==="
        echo "rs1042725 (HMGA2): $(snp_info 12:66220748)"  # TODO(build): unverified GRCh38
        echo "rs143383 (GDF5): $(snp_info 20:35437977)"
        echo ""
        
    } > "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# 4. ROZSZERZONE RYZYKO CHORÓB
#===============================================================================

run_health_risks() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  ROZSZERZONE RYZYKO CHORÓB"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/health/health_risks_extended.txt"
    
    {
        echo "ROZSZERZONY RAPORT RYZYKA CHORÓB"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        # CHOROBA PARKINSONA
        echo "=== CHOROBA PARKINSONA ==="
        echo "rs34637584 (LRRK2 G2019S): $(snp_info 12:40734202)"  # TODO(build): unverified GRCh38
        echo "rs76763715 (GBA N370S): $(snp_info 1:155235843)"
        echo ""
        
        # CHOROBA ALZHEIMERA
        echo "=== CHOROBA ALZHEIMERA ==="
        echo "rs429358 (APOE): $(snp_info 19:44908684)"
        echo "rs7412 (APOE): $(snp_info 19:44908822)"
        echo "rs75932628 (TREM2): $(snp_info 6:41129252)"  # TODO(build): unverified GRCh38
        echo ""
        
        # STWARDNIENIE ROZSIANE (MS)
        echo "=== STWARDNIENIE ROZSIANE ==="
        echo "rs3135388 (HLA-DRB1*15:01 proxy): $(snp_info 6:32557515)"  # TODO(build): unverified GRCh38
        echo "rs6897932 (IL7R): $(snp_info 5:35874575)"  # TODO(build): unverified GRCh38
        echo ""
        
        # CELIAKIA
        echo "=== CELIAKIA ==="
        echo "rs2187668 (HLA-DQ2.5): $(snp_info 6:32713862)"  # TODO(build): unverified GRCh38
        echo "rs7454108 (HLA-DQ8): $(snp_info 6:32681277)"  # TODO(build): unverified GRCh38
        echo ""
        
        # JASKRA
        echo "=== JASKRA ==="
        echo "rs10483727 (SIX1/SIX6): $(snp_info 14:60974593)"  # TODO(build): unverified GRCh38
        echo "rs4236601 (CAV1/CAV2): $(snp_info 7:116164078)"  # TODO(build): unverified GRCh38
        echo ""
        
        # AMD (zwyrodnienie plamki)
        echo "=== ZWYRODNIENIE PLAMKI (AMD) ==="
        echo "rs1061170 (CFH Y402H): $(snp_info 1:196659237)"  # TODO(build): unverified GRCh38
        echo "rs10490924 (ARMS2): $(snp_info 10:124214448)"  # TODO(build): unverified GRCh38
        echo ""
        
        # ŁUSZCZYCA
        echo "=== ŁUSZCZYCA ==="
        echo "rs12191877 (HLA-C): $(snp_info 6:31271836)"  # TODO(build): unverified GRCh38
        echo "rs4406273 (IL12B): $(snp_info 5:158818745)"  # TODO(build): unverified GRCh38
        echo ""
        
        # MIGRENY
        echo "=== MIGRENY ==="
        echo "rs2651899 (PRDM16): $(snp_info 1:3094240)"  # TODO(build): unverified GRCh38
        echo "rs10166942 (TRPM8): $(snp_info 2:234824723)"  # TODO(build): unverified GRCh38
        echo ""
        
        # ASTMA
        echo "=== ASTMA ==="
        echo "rs8076131 (ORMDL3): $(snp_info 17:38122679)"  # TODO(build): unverified GRCh38
        echo "rs2305480 (GSDMB): $(snp_info 17:38072472)"  # TODO(build): unverified GRCh38
        echo ""
        
        # DEPRESJA
        echo "=== DEPRESJA (wybrane markery) ==="
        echo "rs6265 (BDNF Val66Met): $(snp_info 11:27658369)"
        echo "rs4680 (COMT Val158Met): $(snp_info 22:19963748)"
        echo "rs25531 (SLC6A4 - serotonin transporter): $(snp_info 17:28564346)"  # TODO(build): unverified GRCh38
        echo ""
        
        # OSTEOPOROZA
        echo "=== OSTEOPOROZA ==="
        echo "rs2282679 (GC/DBP): $(snp_info 4:72618323)"  # TODO(build): unverified GRCh38
        echo "rs4988235 (LCT - wchłanianie wapnia): $(snp_info 2:135851076)"
        echo ""
        
        # TARCZYCA
        echo "=== CHOROBY TARCZYCY ==="
        echo "rs965513 (FOXE1 - rak tarczycy): $(snp_info 9:100556109)"  # TODO(build): unverified GRCh38
        echo "rs3184504 (SH2B3 - Hashimoto): $(snp_info 12:111884608)"  # TODO(build): unverified GRCh38
        echo ""
        
    } > "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# 5. mtDNA ANALIZA
#===============================================================================

run_mtdna() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  ANALIZA mtDNA (haplogrupa matczyna)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/ancestry/mtdna_report.txt"
    
    {
        echo "RAPORT mtDNA"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        # Sprawdź czy mamy MT w VCF
        MT_COUNT=$(hx_query MT '%POS\n' | wc -l | tr -d ' ')
        
        if [ "$MT_COUNT" -gt "0" ]; then
            echo "Znaleziono $MT_COUNT wariantów mtDNA"
            echo ""
            
            # Wypisz wszystkie warianty MT
            echo "=== WARIANTY mtDNA ==="
            hx_query MT 'm.%POS %REF>%ALT\n'
            echo ""
            
            # Kluczowe pozycje dla haplogrup
            echo "=== KLUCZOWE POZYCJE HAPLOGRUP ==="
            echo "Pozycja 263 (A>G = większość nie-afrykańskich): $(snp_info MT:263-263 'brak wariantu = A')"
            echo "Pozycja 7028 (C>T = haplogrupa H): $(snp_info MT:7028-7028 'brak wariantu')"
            echo "Pozycja 12705 (C>T = haplogrupa R): $(snp_info MT:12705-12705 'brak wariantu')"
            echo "Pozycja 10398 (A>G = haplogrupy I, J, K): $(snp_info MT:10398-10398 'brak wariantu')"
            echo ""
            
            # Warianty patogenne mtDNA
            echo "=== WARIANTY PATOGENNE mtDNA (sprawdź) ==="
            echo "m.3243 (MELAS): $(snp_info MT:3243-3243 'OK - brak')"
            echo "m.8344 (MERRF): $(snp_info MT:8344-8344 'OK - brak')"
            echo "m.11778 (LHON): $(snp_info MT:11778-11778 'OK - brak')"
            echo "m.14484 (LHON): $(snp_info MT:14484-14484 'OK - brak')"
            echo ""
            
            echo "TIP: Wklej warianty mtDNA do https://haplogrep.i-med.ac.at/"
            echo "     żeby uzyskać dokładną haplogrupę"
            
        else
            echo "Brak wariantów MT w VCF."
            echo ""
            echo "Spróbuj wyekstrahować z BAM:"
            echo "  bcftools mpileup -r MT -f ref.fa $BAM | bcftools call -mv > mt_variants.vcf"
        fi
        
    } > "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# 6. Y-DNA ANALIZA
#===============================================================================

run_ydna() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  ANALIZA Y-DNA (haplogrupa ojcowska)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/ancestry/ydna_report.txt"
    
    {
        echo "RAPORT Y-DNA"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        # Sprawdź czy mamy Y w VCF
        Y_COUNT=$(hx_query Y '%POS\n' | wc -l | tr -d ' ')
        
        if [ "$Y_COUNT" -gt "100" ]; then
            echo "Znaleziono $Y_COUNT wariantów Y-DNA"
            echo "Płeć: MĘŻCZYZNA"
            echo ""
            
            # Kluczowe SNP dla głównych haplogrup
            echo "=== GŁÓWNE HAPLOGRUPY Y-DNA ==="
            
            # Haplogrupa R (Europa Zachodnia)
            echo ""
            echo "Haplogrupa R (Europa Zachodnia, częsta w Polsce):"
            echo "  M207 (R): $(snp_info Y:14969634)"  # TODO(build): unverified GRCh38
            echo "  M343 (R1b): $(snp_info Y:22801040)"  # TODO(build): unverified GRCh38
            echo "  M269 (R1b1a1a2): $(snp_info Y:22739367)"  # TODO(build): unverified GRCh38
            echo "  M198/M417 (R1a): $(snp_info Y:16867800)"  # TODO(build): unverified GRCh38
            
            # Haplogrupa I (Północna Europa)
            echo ""
            echo "Haplogrupa I (Skandynawia, Bałkany):"
            echo "  M170 (I): $(snp_info Y:15597896)"  # TODO(build): unverified GRCh38
            echo "  M253 (I1): $(snp_info Y:13532401)"  # TODO(build): unverified GRCh38
            echo "  M438 (I2): $(snp_info Y:8579240)"  # TODO(build): unverified GRCh38
            
            # Haplogrupa E (Afryka, Bliski Wschód)
            echo ""
            echo "Haplogrupa E (Afryka, Bliski Wschód):"
            echo "  M96 (E): $(snp_info Y:21867774)"  # TODO(build): unverified GRCh38
            
            # Haplogrupa J (Bliski Wschód)
            echo ""
            echo "Haplogrupa J (Bliski Wschód, Żydzi, Arabowie):"
            echo "  M304 (J): $(snp_info Y:14969634)"  # TODO(build): unverified GRCh38
            
            # Haplogrupa N (Finlandia, Syberia)
            echo ""
            echo "Haplogrupa N (Finlandia, Syberia, Bałtyk):"
            echo "  M231 (N): $(snp_info Y:8116869)"  # TODO(build): unverified GRCh38
            
            echo ""
            echo "TIP: Użyj narzędzia yhaplo lub wklej dane do:"
            echo "     https://www.yfull.com/"
            
        else
            echo "Brak lub mało wariantów Y-DNA."
            echo "Możliwe przyczyny:"
            echo "  - Próbka od kobiety (XX)"
            echo "  - Niskie pokrycie chromosomu Y"
            echo ""
            
            # Sprawdź stosunek X:Y
            X_COUNT=$(hx_query X '%POS\n' | wc -l | tr -d ' ')
            echo "Warianty X: $X_COUNT"
            echo "Warianty Y: $Y_COUNT"
            
            if [ "$X_COUNT" -gt "1000" ] && [ "$Y_COUNT" -lt "100" ]; then
                echo ""
                echo "Wniosek: Prawdopodobnie KOBIETA (XX)"
            fi
        fi
        
    } > "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# 7. RUNS OF HOMOZYGOSITY (ROH)
#===============================================================================

run_roh() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  RUNS OF HOMOZYGOSITY (ROH) - Pokrewieństwo rodziców"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/ancestry/roh_report.txt"
    
    {
        echo "RAPORT ROH (Runs of Homozygosity)"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        echo "Obliczanie ROH (to może potrwać)..."
        
    } > "$OUT"
    
    # Uruchom bcftools roh
    bcftools roh -G30 --AF-dflt 0.4 "$VCF" 2>"$OUTPUT_DIR/ancestry/roh.log" | \
        grep "^RG" > "$OUTPUT_DIR/ancestry/roh_regions.txt" || true
    
    {
        echo ""
        
        # Analiza wyników
        if [ -s "$OUTPUT_DIR/ancestry/roh_regions.txt" ]; then
            TOTAL_ROH=$(awk '{sum+=$6} END {print sum}' "$OUTPUT_DIR/ancestry/roh_regions.txt" || echo "0")
            NUM_REGIONS=$(wc -l < "$OUTPUT_DIR/ancestry/roh_regions.txt" | tr -d ' ')
            
            echo "=== WYNIKI ROH ==="
            echo "Liczba regionów ROH: $NUM_REGIONS"
            echo "Całkowita długość ROH: $TOTAL_ROH bp"
            echo ""
            
            # Współczynnik inbredu (przybliżony)
            # F ≈ ROH_total / 3,000,000,000 (wielkość genomu)
            if [ -n "$TOTAL_ROH" ] && [ "$TOTAL_ROH" -gt "0" ]; then
                F_COEF=$(awk -v t="$TOTAL_ROH" 'BEGIN {printf "%.6f", t / 3000000000}')
                echo "Przybliżony współczynnik inbredu (F): $F_COEF"
                echo ""
                
                echo "=== INTERPRETACJA ==="
                echo "F < 0.01 (ROH < 30 Mb): Brak bliskiego pokrewieństwa rodziców"
                echo "F = 0.01-0.03: Możliwe dalekie pokrewieństwo"
                echo "F = 0.03-0.06: Kuzynostwo 2-3 stopnia"
                echo "F > 0.06: Bliskie pokrewieństwo rodziców"
            fi
            
            echo ""
            echo "=== NAJWIĘKSZE REGIONY ROH ==="
            # No `head` in the middle: under pipefail its early exit would kill
            # sort with SIGPIPE (rc 141) once the input exceeds the pipe buffer.
            sort -k6 -rn "$OUTPUT_DIR/ancestry/roh_regions.txt" | \
                awk 'NR<=10 {printf "chr%s:%s-%s (%d bp)\n", $3, $4, $5, $6}'
            
        else
            echo "Nie znaleziono znaczących regionów ROH."
            echo "To sugeruje brak bliskiego pokrewieństwa rodziców."
        fi
        
    } >> "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# 8. STATYSTYKI GENOMU
#===============================================================================

run_stats() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  STATYSTYKI GENOMU"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    OUT="$OUTPUT_DIR/reports/genome_stats.txt"
    
    {
        echo "STATYSTYKI GENOMU"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        echo "=== PODSUMOWANIE WARIANTÓW ==="
        TOTAL=$(bcftools view -H "$VCF" | wc -l | tr -d ' ')
        SNPS=$(bcftools view -v snps -H "$VCF" | wc -l | tr -d ' ')
        INDELS=$(bcftools view -v indels -H "$VCF" | wc -l | tr -d ' ')
        
        echo "Wszystkie warianty: $TOTAL"
        echo "SNPs: $SNPS"
        echo "Indele: $INDELS"
        echo ""
        
        # Heterozygotyczność
        # hx_count zwraca zawsze jedną liczbę (grep -c || echo dawało "0\n0")
        GT_LIST="$OUTPUT_DIR/reports/genotypes.tmp"
        bcftools query -f '[%GT]\n' "$VCF" > "$GT_LIST"
        HET=$(hx_count "$GT_LIST" '^0[/|]1$')
        HOM_ALT=$(hx_count "$GT_LIST" '^1[/|]1$')
        rm -f "$GT_LIST"
        
        echo "=== ZYGOTYCZNOŚĆ ==="
        echo "Heterozygoty (0/1): $HET"
        echo "Homozygoty alt (1/1): $HOM_ALT"
        
        if [ "$TOTAL" -gt "0" ]; then
            HET_RATIO=$(awk -v h="$HET" -v t="$TOTAL" 'BEGIN {printf "%.4f", h / t}')
            echo "Stosunek het/total: $HET_RATIO"
        fi
        echo ""
        
        # Ti/Tv ratio
        echo "=== Ti/Tv RATIO ==="
        # A single-file stats run emits exactly one TSTV line; no `exit` in awk,
        # which would SIGPIPE bcftools under pipefail on whole-genome input.
        bcftools stats "$VCF" | \
            awk '/^TSTV\t0\t/ {printf "Transitions: %s\nTransversions: %s\nTi/Tv: %s\n", $3, $4, $5}'
        echo "(Oczekiwane dla WGS: ~2.0-2.1)"
        echo ""
        
        # Warianty per chromosom
        echo "=== WARIANTY PER CHROMOSOM ==="
        for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
            count=$(hx_query "$chr" '%POS\n' | wc -l | tr -d ' ')
            printf "chr%-2s: %'d\n" "$chr" "$count"
        done
        echo ""
        
        # Jakość wariantów
        echo "=== ROZKŁAD JAKOŚCI ==="
        echo "QUAL < 20:  $(bcftools view -H "$VCF" | awk '$6 < 20' | wc -l | tr -d ' ')"
        echo "QUAL 20-50: $(bcftools view -H "$VCF" | awk '$6 >= 20 && $6 < 50' | wc -l | tr -d ' ')"
        echo "QUAL 50-100: $(bcftools view -H "$VCF" | awk '$6 >= 50 && $6 < 100' | wc -l | tr -d ' ')"
        echo "QUAL >= 100: $(bcftools view -H "$VCF" | awk '$6 >= 100' | wc -l | tr -d ' ')"
        echo ""
        
    } > "$OUT"
    
    echo "✓ Zapisano: $OUT"
    cat "$OUT"
}

#===============================================================================
# WYKONANIE
#===============================================================================

case $choice in
    1) run_pharmacogenomics ;;
    2) run_ancestry ;;
    3) run_traits ;;
    4) run_health_risks ;;
    5) run_mtdna ;;
    6) run_ydna ;;
    7) run_roh ;;
    8) run_stats ;;
    9) 
        run_pharmacogenomics
        run_ancestry
        run_traits
        run_health_risks
        run_mtdna
        run_ydna
        run_roh
        run_stats
        ;;
    0) echo "Wyjście." ;;
    *) echo "Nieprawidłowy wybór" ;;
esac

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Wszystkie raporty zapisane w: $OUTPUT_DIR"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
hx_report_coverage
