#!/bin/bash
#===============================================================================
# ANOTACJA VCF BAZA dbSNP
# Dodaje rsID do wariantow
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"
OUTPUT="${VCF%.vcf.gz}_annotated.vcf.gz"

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    ANOTACJA VCF BAZA dbSNP                                ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Sprawdz czy mamy dbSNP
DBSNP="dbsnp156_grch38.vcf.gz"

if [ ! -f "$DBSNP" ]; then
    echo "📥 Pobieram dbSNP (to zajmie kilka minut, plik ~3GB)..."
    echo ""
    echo "UWAGA: Pelny dbSNP jest duzy. Mozesz tez uzyc mniejszej wersji:"
    echo "  - common_all: tylko czeste warianty (~600MB)"
    echo ""
    
    # Pobierz common_all (mniejszy, szybszy)
    wget -q --show-progress \
        "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz" \
        -O dbsnp_common.vcf.gz
    
    wget -q \
        "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz.tbi" \
        -O dbsnp_common.vcf.gz.tbi
    
    DBSNP="dbsnp_common.vcf.gz"
    echo "✅ dbSNP pobrany!"
fi

echo ""
echo "🔧 Anotuję VCF..."
echo "   Input:  $VCF"
echo "   dbSNP:  $DBSNP"
echo "   Output: $OUTPUT"
echo ""

# Anotuj
bcftools annotate \
    -a "$DBSNP" \
    -c ID \
    -O z \
    -o "$OUTPUT" \
    "$VCF"

# Indeksuj
bcftools index -t "$OUTPUT"

echo ""
echo "✅ Anotacja zakonczona!"
echo ""

# Statystyki
TOTAL=$(bcftools view -H "$OUTPUT" | wc -l)
WITH_RS=$(bcftools query -f '%ID\n' "$OUTPUT" | grep -c "^rs")
WITHOUT_RS=$((TOTAL - WITH_RS))

echo "📊 STATYSTYKI:"
echo "   Calkowite warianty:    $TOTAL"
echo "   Z rsID (znane):        $WITH_RS ($(echo "scale=1; $WITH_RS * 100 / $TOTAL" | bc)%)"
echo "   Bez rsID (potenc. nowe): $WITHOUT_RS ($(echo "scale=1; $WITHOUT_RS * 100 / $TOTAL" | bc)%)"
echo ""
echo "Plik wynikowy: $OUTPUT"
