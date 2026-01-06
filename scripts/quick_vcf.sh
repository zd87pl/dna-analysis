#!/bin/bash
#===============================================================================
# SZYBKIE GENEROWANIE VCF Z BAM
# Minimalna wersja - tylko VCF
#===============================================================================

BAM="$1"
OUTPUT="${2:-saryd_variants.vcf.gz}"
REFERENCE="${3:-}"
THREADS="${4:-8}"

if [ -z "$BAM" ]; then
    echo "Użycie: $0 <plik.bam> [output.vcf.gz] [reference.fa] [threads]"
    echo ""
    echo "Przykłady:"
    echo "  $0 saryd.bam"
    echo "  $0 saryd.bam my_variants.vcf.gz human_g1k_v37.fasta 16"
    exit 1
fi

echo "═══════════════════════════════════════════════════════════════"
echo "         SZYBKIE GENEROWANIE VCF"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "BAM:        $BAM"
echo "Output:     $OUTPUT"
echo "Reference:  ${REFERENCE:-BRAK (ograniczone możliwości)}"
echo "Threads:    $THREADS"
echo ""

# Sprawdź czy BAM jest zindeksowany
if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.*}.bai" ]; then
    echo "[1/4] Indeksowanie BAM..."
    samtools index -@ $THREADS "$BAM"
else
    echo "[1/4] Index BAM istnieje ✓"
fi

# Variant calling
echo "[2/4] Variant calling (to może potrwać 2-6h dla WGS)..."
echo "      Rozpoczęto: $(date)"

if [ -n "$REFERENCE" ] && [ -f "$REFERENCE" ]; then
    # Z referencją - pełny variant calling
    bcftools mpileup \
        -Ou \
        -f "$REFERENCE" \
        --threads "$THREADS" \
        --max-depth 250 \
        --min-MQ 20 \
        --min-BQ 20 \
        "$BAM" 2>/dev/null | \
    bcftools call \
        -mv \
        --threads "$THREADS" \
        -Oz \
        -o "$OUTPUT"
else
    echo "UWAGA: Brak referencji - używam trybu uproszczonego"
    echo ""
    echo "Pobierz referencję GRCh37:"
    echo "  wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
    echo "  gunzip human_g1k_v37.fasta.gz"
    echo "  samtools faidx human_g1k_v37.fasta"
    echo ""
    
    # Bez referencji - tylko statystyki
    echo "Generuję statystyki pokrycia zamiast VCF..."
    samtools depth "$BAM" | \
        awk 'BEGIN{OFS="\t"} {sum[$1]+=$3; cnt[$1]++} END{for(c in sum) print c, sum[c]/cnt[c]}' \
        > "${OUTPUT%.vcf.gz}_coverage_per_chrom.txt"
    
    echo "Zapisano: ${OUTPUT%.vcf.gz}_coverage_per_chrom.txt"
    exit 0
fi

echo "      Zakończono: $(date)"

# Indeksowanie VCF
echo "[3/4] Indeksowanie VCF..."
tabix -p vcf "$OUTPUT"

# Statystyki
echo "[4/4] Generowanie statystyk..."
bcftools stats "$OUTPUT" > "${OUTPUT%.vcf.gz}_stats.txt"

# Podsumowanie
echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "                      GOTOWE!"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "Wygenerowane pliki:"
echo "  • $OUTPUT"
echo "  • ${OUTPUT}.tbi"
echo "  • ${OUTPUT%.vcf.gz}_stats.txt"
echo ""

# Szybkie statystyki
TOTAL=$(bcftools view -H "$OUTPUT" | wc -l)
SNPS=$(bcftools view -v snps -H "$OUTPUT" | wc -l)
INDELS=$(bcftools view -v indels -H "$OUTPUT" | wc -l)

echo "Statystyki wariantów:"
echo "  • Wszystkie: $TOTAL"
echo "  • SNPs:      $SNPS"
echo "  • Indele:    $INDELS"
echo ""
echo "Podgląd pierwszych wariantów:"
bcftools view -H "$OUTPUT" | head -5
