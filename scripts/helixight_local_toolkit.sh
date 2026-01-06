#!/bin/bash
#===============================================================================
# HELIXIGHT LOCAL ANALYSIS TOOLKIT
# Kompleksowa analiza VCF na macOS
#===============================================================================

set -e

VCF="${1:-saryd_variants.vcf.gz}"
BAM="${2:-saryd.bam}"
OUTPUT_DIR="${3:-./helixight_analysis}"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║          HELIXIGHT LOCAL ANALYSIS TOOLKIT                         ║"
echo "║          Analiza VCF + BAM na macOS                               ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "VCF: $VCF"
echo "BAM: $BAM"
echo "Output: $OUTPUT_DIR"
echo ""

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

run_all=false
[ "$choice" == "9" ] && run_all=true

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
        echo "Lokalizacja: chr22:42522500-42526883"
        bcftools view -r 22:42522500-42526883 "$VCF" 2>/dev/null | grep -v "^#" | \
            awk '{print "  " $1 ":" $2 " " $4 ">" $5 " (QUAL=" $6 ")"}' || echo "  Brak wariantów w regionie"
        echo ""
        
        # Kluczowe SNP CYP2D6
        echo "Kluczowe SNP CYP2D6:"
        for snp in "22:42524947" "22:42525772" "22:42526694"; do
            result=$(bcftools view -r "$snp" "$VCF" 2>/dev/null | grep -v "^#" | head -1)
            if [ -n "$result" ]; then
                echo "  $snp: $(echo $result | awk '{print $4 ">" $5}')"
            fi
        done
        echo ""
        
        # CYP2C19
        echo "=== CYP2C19 (klopidogrel, IPP, SSRI) ==="
        echo "rs4244285 (*2): $(bcftools view -r 10:96541616 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs12248560 (*17): $(bcftools view -r 10:96521657 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # CYP2C9
        echo "=== CYP2C9 (warfaryna, NLPZ) ==="
        echo "rs1799853 (*2): $(bcftools view -r 10:96702047 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs1057910 (*3): $(bcftools view -r 10:96741053 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # CYP3A4/CYP3A5
        echo "=== CYP3A4/CYP3A5 (statyny, immunosupresja) ==="
        echo "CYP3A5*3 rs776746: $(bcftools view -r 7:99672916 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "CYP3A4*22 rs35599367: $(bcftools view -r 7:99361466 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # VKORC1
        echo "=== VKORC1 (warfaryna) ==="
        echo "rs9923231: $(bcftools view -r 16:31107689 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # SLCO1B1
        echo "=== SLCO1B1 (statyny - miopatia) ==="
        echo "rs4149056 (*5): $(bcftools view -r 12:21331549 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # DPYD
        echo "=== DPYD (chemioterapia 5-FU) ==="
        echo "rs3918290 (*2A): $(bcftools view -r 1:97915614 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs55886062 (*13): $(bcftools view -r 1:98348885 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # TPMT
        echo "=== TPMT (tiopuryny - azatiopryna) ==="
        echo "rs1800460 (*3B): $(bcftools view -r 6:18130918 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs1142345 (*3C): $(bcftools view -r 6:18139228 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # HLA - reakcje na leki
        echo "=== HLA (ciężkie reakcje skórne) ==="
        echo "Region HLA-B: chr6:31321649-31324989"
        bcftools view -r 6:31321649-31324989 "$VCF" 2>/dev/null | grep -v "^#" | wc -l | \
            xargs -I {} echo "  Znaleziono {} wariantów w regionie HLA-B"
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
        echo "  rs1426654 (SLC24A5 - jasna skóra): $(bcftools view -r 15:48426484 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  rs16891982 (SLC45A2 - jasna skóra): $(bcftools view -r 5:33951693 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  rs12913832 (HERC2 - niebieskie oczy): $(bcftools view -r 15:28365618 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # Azjatyckie markery
        echo "Markery azjatyckie:"
        echo "  rs3827760 (EDAR - grube włosy): $(bcftools view -r 2:109513601 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  rs671 (ALDH2 - Asian flush): $(bcftools view -r 12:112241766 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # Afrykańskie markery
        echo "Markery afrykańskie:"
        echo "  rs2814778 (DARC - malaria): $(bcftools view -r 1:159174683 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # Tolerancja laktozy (różnice populacyjne)
        echo "Adaptacje populacyjne:"
        echo "  rs4988235 (LCT - tolerancja laktozy): $(bcftools view -r 2:136608646 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  rs182549 (LCT - alt marker): $(bcftools view -r 2:136616754 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
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
        echo "rs12913832 (HERC2 - główny): $(bcftools view -r 15:28365618 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  AA = brązowe oczy (dominujące)"
        echo "  AG = brązowe/zielone"
        echo "  GG = niebieskie oczy"
        echo ""
        echo "rs1800407 (OCA2): $(bcftools view -r 15:28230318 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs12896399 (SLC24A4): $(bcftools view -r 14:92773663 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # KOLOR WŁOSÓW
        echo "=== KOLOR WŁOSÓW ==="
        echo "rs12821256 (KITLG - blond): $(bcftools view -r 12:89328335 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs1805007 (MC1R - rude): $(bcftools view -r 16:89985844 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs1805008 (MC1R - rude): $(bcftools view -r 16:89986091 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs1805009 (MC1R - rude): $(bcftools view -r 16:89986144 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # ŁYSIENIE
        echo "=== ŁYSIENIE (androgenowe) ==="
        echo "rs2180439 (chr20): $(bcftools view -r 20:21985252 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs6152 (AR - receptor androgenowy): $(bcftools view -r X:66943551 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # PIEGI
        echo "=== PIEGI ==="
        echo "rs4911414 (IRF4): $(bcftools view -r 6:396321 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs1805007 (MC1R): już powyżej"
        echo ""
        
        # GORZKI SMAK
        echo "=== PERCEPCJA SMAKU ==="
        echo "rs713598 (TAS2R38 - gorzki smak): $(bcftools view -r 7:141672604 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs1726866 (TAS2R38): $(bcftools view -r 7:141672705 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  PAV/PAV = silnie czujesz gorzki smak"
        echo "  AVI/AVI = nie czujesz gorzkiego"
        echo ""
        
        # KICHANIE NA SŁOŃCE
        echo "=== INNE CECHY ==="
        echo "rs10427255 (kichanie na słońce - photic sneeze): $(bcftools view -r 2:135851076 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs17822931 (ABCC11 - typ woskowiny): $(bcftools view -r 16:48258198 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  CC = sucha woskowina (azjatycki typ)"
        echo "  CT/TT = mokra woskowina (europejski typ)"
        echo ""
        
        # KOFEINA
        echo "=== METABOLIZM KOFEINY ==="
        echo "rs762551 (CYP1A2): $(bcftools view -r 15:75041917 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "  AA = szybki metabolizer (kawa OK)"
        echo "  AC/CC = wolny metabolizer (ogranicz kawę)"
        echo ""
        
        # ALKOHOL
        echo "=== METABOLIZM ALKOHOLU ==="
        echo "rs1229984 (ADH1B): $(bcftools view -r 4:100239319 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs671 (ALDH2 - Asian flush): $(bcftools view -r 12:112241766 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # WZROST
        echo "=== WZROST (wybrane warianty) ==="
        echo "rs1042725 (HMGA2): $(bcftools view -r 12:66220748 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs6060373 (GDF5): $(bcftools view -r 20:34025756 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
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
        echo "rs34637584 (LRRK2 G2019S): $(bcftools view -r 12:40734202 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs76763715 (GBA N370S): $(bcftools view -r 1:155205634 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # CHOROBA ALZHEIMERA
        echo "=== CHOROBA ALZHEIMERA ==="
        echo "rs429358 (APOE): $(bcftools view -r 19:45411941 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs7412 (APOE): $(bcftools view -r 19:45412079 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs75932628 (TREM2): $(bcftools view -r 6:41129252 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # STWARDNIENIE ROZSIANE (MS)
        echo "=== STWARDNIENIE ROZSIANE ==="
        echo "rs3135388 (HLA-DRB1*15:01 proxy): $(bcftools view -r 6:32557515 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs6897932 (IL7R): $(bcftools view -r 5:35874575 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # CELIAKIA
        echo "=== CELIAKIA ==="
        echo "rs2187668 (HLA-DQ2.5): $(bcftools view -r 6:32713862 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs7454108 (HLA-DQ8): $(bcftools view -r 6:32681277 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # JASKRA
        echo "=== JASKRA ==="
        echo "rs10483727 (SIX1/SIX6): $(bcftools view -r 14:60974593 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs4236601 (CAV1/CAV2): $(bcftools view -r 7:116164078 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # AMD (zwyrodnienie plamki)
        echo "=== ZWYRODNIENIE PLAMKI (AMD) ==="
        echo "rs1061170 (CFH Y402H): $(bcftools view -r 1:196659237 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs10490924 (ARMS2): $(bcftools view -r 10:124214448 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # ŁUSZCZYCA
        echo "=== ŁUSZCZYCA ==="
        echo "rs12191877 (HLA-C): $(bcftools view -r 6:31271836 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs4406273 (IL12B): $(bcftools view -r 5:158818745 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # MIGRENY
        echo "=== MIGRENY ==="
        echo "rs2651899 (PRDM16): $(bcftools view -r 1:3094240 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs10166942 (TRPM8): $(bcftools view -r 2:234824723 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # ASTMA
        echo "=== ASTMA ==="
        echo "rs8076131 (ORMDL3): $(bcftools view -r 17:38122679 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs2305480 (GSDMB): $(bcftools view -r 17:38072472 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # DEPRESJA
        echo "=== DEPRESJA (wybrane markery) ==="
        echo "rs6265 (BDNF Val66Met): $(bcftools view -r 11:27679916 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs4680 (COMT Val158Met): $(bcftools view -r 22:19951271 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs25531 (SLC6A4 - serotonin transporter): $(bcftools view -r 17:28564346 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # OSTEOPOROZA
        echo "=== OSTEOPOROZA ==="
        echo "rs2282679 (GC/DBP): $(bcftools view -r 4:72618323 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs4988235 (LCT - wchłanianie wapnia): $(bcftools view -r 2:136608646 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo ""
        
        # TARCZYCA
        echo "=== CHOROBY TARCZYCY ==="
        echo "rs965513 (FOXE1 - rak tarczycy): $(bcftools view -r 9:100556109 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
        echo "rs3184504 (SH2B3 - Hashimoto): $(bcftools view -r 12:111884608 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
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
        MT_COUNT=$(bcftools view -r MT "$VCF" 2>/dev/null | grep -v "^#" | wc -l || echo "0")
        MT_COUNT=$(echo $MT_COUNT | tr -d ' ')
        
        if [ "$MT_COUNT" -gt "0" ]; then
            echo "Znaleziono $MT_COUNT wariantów mtDNA"
            echo ""
            
            # Wypisz wszystkie warianty MT
            echo "=== WARIANTY mtDNA ==="
            bcftools view -r MT "$VCF" 2>/dev/null | grep -v "^#" | \
                awk '{printf "m.%s %s>%s\n", $2, $4, $5}'
            echo ""
            
            # Kluczowe pozycje dla haplogrup
            echo "=== KLUCZOWE POZYCJE HAPLOGRUP ==="
            echo "Pozycja 263 (A>G = większość nie-afrykańskich): $(bcftools view -r MT:263-263 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'brak wariantu = A')"
            echo "Pozycja 7028 (C>T = haplogrupa H): $(bcftools view -r MT:7028-7028 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'brak wariantu')"
            echo "Pozycja 12705 (C>T = haplogrupa R): $(bcftools view -r MT:12705-12705 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'brak wariantu')"
            echo "Pozycja 10398 (A>G = haplogrupy I, J, K): $(bcftools view -r MT:10398-10398 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'brak wariantu')"
            echo ""
            
            # Warianty patogenne mtDNA
            echo "=== WARIANTY PATOGENNE mtDNA (sprawdź) ==="
            echo "m.3243 (MELAS): $(bcftools view -r MT:3243-3243 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'OK - brak')"
            echo "m.8344 (MERRF): $(bcftools view -r MT:8344-8344 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'OK - brak')"
            echo "m.11778 (LHON): $(bcftools view -r MT:11778-11778 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'OK - brak')"
            echo "m.14484 (LHON): $(bcftools view -r MT:14484-14484 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'OK - brak')"
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
        Y_COUNT=$(bcftools view -r Y "$VCF" 2>/dev/null | grep -v "^#" | wc -l || echo "0")
        Y_COUNT=$(echo $Y_COUNT | tr -d ' ')
        
        if [ "$Y_COUNT" -gt "100" ]; then
            echo "Znaleziono $Y_COUNT wariantów Y-DNA"
            echo "Płeć: MĘŻCZYZNA"
            echo ""
            
            # Kluczowe SNP dla głównych haplogrup
            echo "=== GŁÓWNE HAPLOGRUPY Y-DNA ==="
            
            # Haplogrupa R (Europa Zachodnia)
            echo ""
            echo "Haplogrupa R (Europa Zachodnia, częsta w Polsce):"
            echo "  M207 (R): $(bcftools view -r Y:14969634 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            echo "  M343 (R1b): $(bcftools view -r Y:22801040 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            echo "  M269 (R1b1a1a2): $(bcftools view -r Y:22739367 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            echo "  M198/M417 (R1a): $(bcftools view -r Y:16867800 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            
            # Haplogrupa I (Północna Europa)
            echo ""
            echo "Haplogrupa I (Skandynawia, Bałkany):"
            echo "  M170 (I): $(bcftools view -r Y:15597896 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            echo "  M253 (I1): $(bcftools view -r Y:13532401 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            echo "  M438 (I2): $(bcftools view -r Y:8579240 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            
            # Haplogrupa E (Afryka, Bliski Wschód)
            echo ""
            echo "Haplogrupa E (Afryka, Bliski Wschód):"
            echo "  M96 (E): $(bcftools view -r Y:21867774 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            
            # Haplogrupa J (Bliski Wschód)
            echo ""
            echo "Haplogrupa J (Bliski Wschód, Żydzi, Arabowie):"
            echo "  M304 (J): $(bcftools view -r Y:14969634 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            
            # Haplogrupa N (Finlandia, Syberia)
            echo ""
            echo "Haplogrupa N (Finlandia, Syberia, Bałtyk):"
            echo "  M231 (N): $(bcftools view -r Y:8116869 "$VCF" 2>/dev/null | grep -v '^#' | awk '{print $4">"$5}' || echo 'N/A')"
            
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
            X_COUNT=$(bcftools view -r X "$VCF" 2>/dev/null | grep -v "^#" | wc -l || echo "0")
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
    bcftools roh -G30 --AF-dflt 0.4 "$VCF" 2>/dev/null | \
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
                F_COEF=$(echo "scale=6; $TOTAL_ROH / 3000000000" | bc)
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
            sort -k6 -rn "$OUTPUT_DIR/ancestry/roh_regions.txt" | head -10 | \
                awk '{printf "chr%s:%s-%s (%d bp)\n", $3, $4, $5, $6}'
            
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
        TOTAL=$(bcftools view -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
        SNPS=$(bcftools view -v snps -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
        INDELS=$(bcftools view -v indels -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
        
        echo "Wszystkie warianty: $TOTAL"
        echo "SNPs: $SNPS"
        echo "Indele: $INDELS"
        echo ""
        
        # Heterozygotyczność
        HET=$(bcftools view -H "$VCF" 2>/dev/null | grep -c "0/1" || echo "0")
        HOM_ALT=$(bcftools view -H "$VCF" 2>/dev/null | grep -c "1/1" || echo "0")
        
        echo "=== ZYGOTYCZNOŚĆ ==="
        echo "Heterozygoty (0/1): $HET"
        echo "Homozygoty alt (1/1): $HOM_ALT"
        
        if [ "$TOTAL" -gt "0" ]; then
            HET_RATIO=$(echo "scale=4; $HET / $TOTAL" | bc)
            echo "Stosunek het/total: $HET_RATIO"
        fi
        echo ""
        
        # Ti/Tv ratio
        echo "=== Ti/Tv RATIO ==="
        bcftools stats "$VCF" 2>/dev/null | grep "TSTV" | head -1 | \
            awk '{printf "Transitions: %s\nTransversions: %s\nTi/Tv: %s\n", $3, $4, $5}'
        echo "(Oczekiwane dla WGS: ~2.0-2.1)"
        echo ""
        
        # Warianty per chromosom
        echo "=== WARIANTY PER CHROMOSOM ==="
        for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
            count=$(bcftools view -r "$chr" -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
            printf "chr%-2s: %'d\n" "$chr" "$count"
        done
        echo ""
        
        # Jakość wariantów
        echo "=== ROZKŁAD JAKOŚCI ==="
        echo "QUAL < 20:  $(bcftools view -H "$VCF" 2>/dev/null | awk '$6 < 20' | wc -l | tr -d ' ')"
        echo "QUAL 20-50: $(bcftools view -H "$VCF" 2>/dev/null | awk '$6 >= 20 && $6 < 50' | wc -l | tr -d ' ')"
        echo "QUAL 50-100: $(bcftools view -H "$VCF" 2>/dev/null | awk '$6 >= 50 && $6 < 100' | wc -l | tr -d ' ')"
        echo "QUAL >= 100: $(bcftools view -H "$VCF" 2>/dev/null | awk '$6 >= 100' | wc -l | tr -d ' ')"
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
