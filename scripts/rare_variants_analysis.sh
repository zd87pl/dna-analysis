#!/bin/bash
#===============================================================================
#
#    ██████╗  █████╗ ██████╗ ███████╗    ██╗   ██╗ █████╗ ██████╗ ██╗ █████╗ ███╗   ██╗████████╗███████╗
#    ██╔══██╗██╔══██╗██╔══██╗██╔════╝    ██║   ██║██╔══██╗██╔══██╗██║██╔══██╗████╗  ██║╚══██╔══╝██╔════╝
#    ██████╔╝███████║██████╔╝█████╗      ██║   ██║███████║██████╔╝██║███████║██╔██╗ ██║   ██║   ███████╗
#    ██╔══██╗██╔══██║██╔══██╗██╔══╝      ╚██╗ ██╔╝██╔══██║██╔══██╗██║██╔══██║██║╚██╗██║   ██║   ╚════██║
#    ██║  ██║██║  ██║██║  ██║███████╗     ╚████╔╝ ██║  ██║██║  ██║██║██║  ██║██║ ╚████║   ██║   ███████║
#    ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝      ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝
#
#    ANALIZA RZADKOSCI WARIANTOW
#    Jak unikalne jest Twoje DNA?
#
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"
OUTDIR="rare_variants_analysis"
mkdir -p "$OUTDIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║            ANALIZA RZADKOSCI WARIANTOW                                    ║"
echo "║            Jak unikalne jest Twoje DNA?                                   ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Plik VCF: $VCF"
echo ""

#===============================================================================
# KROK 1: Podstawowe statystyki
#===============================================================================

echo "📊 PODSTAWOWE STATYSTYKI WARIANTOW"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Całkowita liczba wariantów
TOTAL=$(bcftools view -H "$VCF" 2>/dev/null | wc -l)
echo "  Calkowita liczba wariantow: $TOTAL"

# SNPs vs Indels
SNPS=$(bcftools view -v snps -H "$VCF" 2>/dev/null | wc -l)
INDELS=$(bcftools view -v indels -H "$VCF" 2>/dev/null | wc -l)
echo "  SNPs: $SNPS"
echo "  Indels: $INDELS"

# Heterozygoty vs Homozygoty
HET=$(bcftools view -H "$VCF" 2>/dev/null | grep -E "0/1|0\|1|1/0|1\|0" | wc -l)
HOM_ALT=$(bcftools view -H "$VCF" 2>/dev/null | grep -E "1/1|1\|1" | wc -l)
echo "  Heterozygoty: $HET"
echo "  Homozygoty ALT: $HOM_ALT"
echo ""

#===============================================================================
# KROK 2: Analiza Ti/Tv ratio
#===============================================================================

echo "📊 JAKOSC DANYCH - Ti/Tv RATIO"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Transitions: A<>G, C<>T
# Transversions: A<>C, A<>T, G<>C, G<>T

TRANSITIONS=$(bcftools view -v snps -H "$VCF" 2>/dev/null | \
    awk '{
        ref=$4; alt=$5;
        if ((ref=="A" && alt=="G") || (ref=="G" && alt=="A") ||
            (ref=="C" && alt=="T") || (ref=="T" && alt=="C"))
            print
    }' | wc -l)

TRANSVERSIONS=$(bcftools view -v snps -H "$VCF" 2>/dev/null | \
    awk '{
        ref=$4; alt=$5;
        if ((ref=="A" && alt=="C") || (ref=="C" && alt=="A") ||
            (ref=="A" && alt=="T") || (ref=="T" && alt=="A") ||
            (ref=="G" && alt=="C") || (ref=="C" && alt=="G") ||
            (ref=="G" && alt=="T") || (ref=="T" && alt=="G"))
            print
    }' | wc -l)

if [ $TRANSVERSIONS -gt 0 ]; then
    TITV_RATIO=$(echo "scale=2; $TRANSITIONS / $TRANSVERSIONS" | bc)
else
    TITV_RATIO="N/A"
fi

echo "  Transitions (Ti): $TRANSITIONS"
echo "  Transversions (Tv): $TRANSVERSIONS"
echo "  Ti/Tv ratio: $TITV_RATIO"
echo ""

if [ "$TITV_RATIO" != "N/A" ]; then
    TITV_NUM=$(echo "$TITV_RATIO" | bc)
    if (( $(echo "$TITV_NUM >= 2.0 && $TITV_NUM <= 2.2" | bc -l) )); then
        echo "  ✅ Ti/Tv w normie dla WGS (oczekiwane: 2.0-2.1)"
    elif (( $(echo "$TITV_NUM >= 2.8 && $TITV_NUM <= 3.3" | bc -l) )); then
        echo "  ✅ Ti/Tv w normie dla WES (oczekiwane: 2.8-3.3)"
    else
        echo "  ⚠️ Ti/Tv poza typowym zakresem - sprawdz jakosc danych"
    fi
fi
echo ""

#===============================================================================
# KROK 3: Warianty wg chromosomow
#===============================================================================

echo "📊 ROZKLAD WARIANTOW WG CHROMOSOMOW"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

printf "  %-6s %10s %10s\n" "Chr" "Warianty" "% calości"
echo "  ─────────────────────────────────"

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
    COUNT=$(bcftools view -r "$CHR" -H "$VCF" 2>/dev/null | wc -l)
    if [ $TOTAL -gt 0 ]; then
        PCT=$(echo "scale=1; $COUNT * 100 / $TOTAL" | bc)
    else
        PCT="0"
    fi
    if [ $COUNT -gt 0 ]; then
        printf "  %-6s %10s %9s%%\n" "$CHR" "$COUNT" "$PCT"
    fi
done
echo ""

#===============================================================================
# KROK 4: Szacowanie rzadkich wariantow (bez gnomAD)
#===============================================================================

echo "📊 SZACOWANIE RZADKOSCI WARIANTOW"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Bez pelnej bazy gnomAD, szacujemy rzadkosc na podstawie:"
echo "- Warianty w regionach kodujacych"
echo "- Warianty z przewidywanym duzym efektem"
echo "- Warianty bez rs ID (potencjalnie unikalne)"
echo ""

# Warianty bez rs ID (novel)
NOVEL=$(bcftools view -H "$VCF" 2>/dev/null | awk '$3=="." || $3==""' | wc -l)
WITH_RS=$(bcftools view -H "$VCF" 2>/dev/null | awk '$3!="." && $3!=""' | wc -l)

echo "  Warianty z rs ID (znane): $WITH_RS"
echo "  Warianty bez rs ID (potencjalnie rzadkie/nowe): $NOVEL"

if [ $TOTAL -gt 0 ]; then
    NOVEL_PCT=$(echo "scale=1; $NOVEL * 100 / $TOTAL" | bc)
    echo "  Procent potencjalnie rzadkich: ${NOVEL_PCT}%"
fi
echo ""

#===============================================================================
# KROK 5: Warianty w genach - regiony kodujace
#===============================================================================

echo "📊 ANALIZA REGIONOW KODUJACYCH"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Kluczowe geny - sprawdz ile wariantow
GENES_TO_CHECK="BRCA1 BRCA2 TP53 APC MLH1 MSH2 LDLR APOB PCSK9 CFTR DMD"

echo "  Warianty w wybranych genach klinicznych:"
echo "  ─────────────────────────────────────────"

# Zdefiniuj pozycje genow (GRCh38)
declare -A GENE_REGIONS=(
    ["BRCA1"]="17:43044295-43125364"
    ["BRCA2"]="13:32315086-32400266"
    ["TP53"]="17:7668421-7687490"
    ["APC"]="5:112707498-112846239"
    ["MLH1"]="3:36993332-37050918"
    ["MSH2"]="2:47403067-47709830"
    ["LDLR"]="19:11089362-11133830"
    ["CFTR"]="7:117465784-117715971"
    ["DMD"]="X:31097677-33339609"
    ["APOE"]="19:44905754-44909393"
    ["ATM"]="11:108222484-108369102"
    ["CHEK2"]="22:28687743-28742422"
)

for GENE in BRCA1 BRCA2 TP53 APC MLH1 MSH2 LDLR CFTR ATM CHEK2 APOE; do
    REGION="${GENE_REGIONS[$GENE]}"
    if [ -n "$REGION" ]; then
        COUNT=$(bcftools view -r "$REGION" -H "$VCF" 2>/dev/null | wc -l)
        printf "  %-10s %5s wariantow\n" "$GENE:" "$COUNT"
    fi
done
echo ""

#===============================================================================
# KROK 6: Typy wariantow wg efektu (na podstawie pozycji)
#===============================================================================

echo "📊 SZACOWANIE EFEKTU WARIANTOW"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Bez pelnej anotacji, szacujemy na podstawie typu wariantu:"
echo ""

# Frameshift candidates (indels nie podzielne przez 3 w regionach kodujacych)
INDELS_1=$(bcftools view -v indels -H "$VCF" 2>/dev/null | \
    awk '{ref=$4; alt=$5; diff=length(alt)-length(ref); if(diff%3!=0) print}' | wc -l)
echo "  Potencjalne frameshift (indel nie %3): $INDELS_1"

# Duże delecje (>50bp)
LARGE_DEL=$(bcftools view -v indels -H "$VCF" 2>/dev/null | \
    awk '{ref=$4; alt=$5; if(length(ref)-length(alt) > 50) print}' | wc -l)
echo "  Duze delecje (>50bp): $LARGE_DEL"

# Duże insercje (>50bp)
LARGE_INS=$(bcftools view -v indels -H "$VCF" 2>/dev/null | \
    awk '{ref=$4; alt=$5; if(length(alt)-length(ref) > 50) print}' | wc -l)
echo "  Duze insercje (>50bp): $LARGE_INS"

echo ""

#===============================================================================
# KROK 7: Porownanie z 1000 Genomes (jesli dostepne AF w VCF)
#===============================================================================

echo "📊 ALLELE FREQUENCY Z VCF (jesli dostepne)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Sprawdz czy VCF ma pole AF
HAS_AF=$(bcftools view -h "$VCF" 2>/dev/null | grep "INFO=<ID=AF" | head -1)

if [ -n "$HAS_AF" ]; then
    echo "  VCF zawiera informacje o czestosci alleli (AF)"
    echo ""
    
    # Rozklad AF
    echo "  Rozklad czestosci alleli:"
    echo "  ─────────────────────────────────────────"
    
    # Ultra-rare (AF < 0.001)
    ULTRA_RARE=$(bcftools view -H "$VCF" 2>/dev/null | \
        bcftools query -f '%INFO/AF\n' 2>/dev/null | \
        awk '$1!="." && $1<0.001' | wc -l)
    
    # Rare (0.001 <= AF < 0.01)
    RARE=$(bcftools view -H "$VCF" 2>/dev/null | \
        bcftools query -f '%INFO/AF\n' 2>/dev/null | \
        awk '$1!="." && $1>=0.001 && $1<0.01' | wc -l)
    
    # Low frequency (0.01 <= AF < 0.05)
    LOW_FREQ=$(bcftools view -H "$VCF" 2>/dev/null | \
        bcftools query -f '%INFO/AF\n' 2>/dev/null | \
        awk '$1!="." && $1>=0.01 && $1<0.05' | wc -l)
    
    # Common (AF >= 0.05)
    COMMON=$(bcftools view -H "$VCF" 2>/dev/null | \
        bcftools query -f '%INFO/AF\n' 2>/dev/null | \
        awk '$1!="." && $1>=0.05' | wc -l)
    
    printf "  %-25s %10s\n" "Ultra-rzadkie (AF<0.1%):" "$ULTRA_RARE"
    printf "  %-25s %10s\n" "Rzadkie (0.1-1%):" "$RARE"
    printf "  %-25s %10s\n" "Niskiej czestosci (1-5%):" "$LOW_FREQ"
    printf "  %-25s %10s\n" "Czeste (>5%):" "$COMMON"
else
    echo "  VCF nie zawiera informacji AF - uzywam szacunkow"
    echo ""
    echo "  Typowy rozklad dla europejskiego genomu:"
    echo "  ─────────────────────────────────────────"
    echo "  • ~150,000-200,000 rzadkich wariantow (AF<1%)"
    echo "  • ~50,000-100,000 ultra-rzadkich (AF<0.1%)"
    echo "  • ~100-500 prywatnych (tylko u Ciebie)"
fi
echo ""

#===============================================================================
# KROK 8: Warianty homozygotyczne rzadkie
#===============================================================================

echo "📊 HOMOZYGOTY - POTENCJALNIE KLINICZNE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Homozygoty w genach chorobowych sa szczegolnie wazne dla chorob recesywnych."
echo ""

# Znajdz homozygoty w genach klinicznych
echo "  Homozygoty ALT w kluczowych genach:"
echo "  ─────────────────────────────────────────"

for GENE in BRCA1 BRCA2 TP53 CFTR ATM MLH1 MSH2; do
    REGION="${GENE_REGIONS[$GENE]}"
    if [ -n "$REGION" ]; then
        HOM_COUNT=$(bcftools view -r "$REGION" -H "$VCF" 2>/dev/null | \
            grep -E "1/1|1\|1" | wc -l)
        if [ $HOM_COUNT -gt 0 ]; then
            printf "  ⚠️ %-10s %5s homozygot ALT\n" "$GENE:" "$HOM_COUNT"
        else
            printf "  ✅ %-10s %5s homozygot\n" "$GENE:" "$HOM_COUNT"
        fi
    fi
done
echo ""

#===============================================================================
# KROK 9: Podsumowanie
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                       PODSUMOWANIE ANALIZY                                ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "  📊 Calkowita liczba wariantow: $TOTAL"
echo "  🧬 SNPs / Indels: $SNPS / $INDELS"
echo "  📈 Ti/Tv ratio: $TITV_RATIO"
echo "  🔍 Warianty bez rs ID: $NOVEL (${NOVEL_PCT:-?}%)"
echo ""
echo "  💡 INTERPRETACJA:"
echo "  ─────────────────────────────────────────"
echo "  • Typowy europejski genom ma ~4-5 mln wariantow vs GRCh38"
echo "  • ~99.9% to warianty neutralne lub lagodne"
echo "  • ~20,000 wariantow zmienia bialka (missense)"
echo "  • ~100-500 to potencjalnie szkodliwe"
echo "  • ~50-100 to warianty zwiazane z chorobami recesywnymi (nosicielstwo)"
echo ""

# Zapisz raport
REPORT="$OUTDIR/rare_variants_report.txt"
cat << EOF > "$REPORT"
================================================================================
              RAPORT ANALIZY RZADKOSCI WARIANTOW
              Wygenerowano: $(date)
================================================================================

PODSTAWOWE STATYSTYKI:
- Calkowita liczba wariantow: $TOTAL
- SNPs: $SNPS
- Indels: $INDELS
- Ti/Tv ratio: $TITV_RATIO
- Heterozygoty: $HET
- Homozygoty ALT: $HOM_ALT

RZADKOSC:
- Warianty z rs ID: $WITH_RS
- Warianty bez rs ID (potencjalnie rzadkie): $NOVEL

================================================================================
EOF

echo "  📁 Raport zapisany: $REPORT"
echo ""
