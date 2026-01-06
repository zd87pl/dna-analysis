#!/bin/bash
#===============================================================================
# TURBO VCF - Parallel Variant Calling dla Apple Silicon
# Zoptymalizowany dla M4 MacBook Air z 24GB RAM
#===============================================================================

set -e

BAM="$1"
OUTPUT="${2:-variants.vcf.gz}"
REFERENCE="${3:-human_g1k_v37.fasta}"
THREADS="${4:-10}"  # M4 Air ma 10 rdzeni (4P + 6E)

# Kolory
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}"
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║         TURBO VCF - Parallel Variant Calling                      ║"
echo "║         Zoptymalizowany dla Apple Silicon M4                      ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Sprawdź argumenty
if [ -z "$BAM" ]; then
    echo "Użycie: $0 <plik.bam> [output.vcf.gz] [reference.fa] [threads]"
    echo ""
    echo "Przykład:"
    echo "  $0 saryd.bam saryd_variants.vcf.gz human_g1k_v37.fasta 10"
    exit 1
fi

# Katalog tymczasowy (użyj SSD!)
TEMP_DIR=$(mktemp -d)
echo -e "${BLUE}Katalog tymczasowy: $TEMP_DIR${NC}"

# Cleanup przy wyjściu
cleanup() {
    echo -e "\n${YELLOW}Sprzątanie plików tymczasowych...${NC}"
    rm -rf "$TEMP_DIR"
}
trap cleanup EXIT

#===============================================================================
# KROK 1: Sprawdzenie i przygotowanie
#===============================================================================

echo -e "\n${YELLOW}[1/5] Sprawdzanie plików...${NC}"

# Sprawdź czy BAM jest zindeksowany
if [ ! -f "${BAM}.bai" ] && [ ! -f "${BAM%.*}.bai" ]; then
    echo "  Indeksowanie BAM..."
    samtools index -@ $THREADS "$BAM"
fi
echo -e "  ${GREEN}✓${NC} BAM zindeksowany"

# Sprawdź referencję
if [ ! -f "$REFERENCE" ]; then
    echo -e "  ${RED}✗${NC} Brak pliku referencyjnego: $REFERENCE"
    exit 1
fi
if [ ! -f "${REFERENCE}.fai" ]; then
    echo "  Indeksowanie referencji..."
    samtools faidx "$REFERENCE"
fi
echo -e "  ${GREEN}✓${NC} Referencja gotowa"

# Pobierz listę chromosomów
CHROMOSOMES=$(samtools idxstats "$BAM" | cut -f1 | grep -v '^\*$' | head -25)
NUM_CHROM=$(echo "$CHROMOSOMES" | wc -l | tr -d ' ')
echo -e "  ${GREEN}✓${NC} Znaleziono $NUM_CHROM chromosomów"

#===============================================================================
# KROK 2: Parallel Variant Calling
#===============================================================================

echo -e "\n${YELLOW}[2/5] Parallel variant calling ($NUM_CHROM chromosomów)...${NC}"
echo "  Startuje: $(date '+%H:%M:%S')"

# Funkcja do przetwarzania jednego chromosomu
process_chromosome() {
    local chrom=$1
    local bam=$2
    local ref=$3
    local temp_dir=$4
    local output_vcf="${temp_dir}/${chrom}.vcf.gz"
    
    # Variant calling dla chromosomu
    bcftools mpileup \
        -Ou \
        -f "$ref" \
        -r "$chrom" \
        --max-depth 200 \
        --min-MQ 20 \
        --min-BQ 20 \
        "$bam" 2>/dev/null | \
    bcftools call \
        -mv \
        -Oz \
        -o "$output_vcf" 2>/dev/null
    
    # Indeksuj
    tabix -p vcf "$output_vcf" 2>/dev/null
    
    echo "  ✓ $chrom"
}

export -f process_chromosome

# Uruchom równolegle (max 8 procesów - balans CPU/RAM)
# Dla M4 Air z 24GB: 8 procesów * ~2GB = 16GB używane
MAX_PARALLEL=8

echo "$CHROMOSOMES" | xargs -P $MAX_PARALLEL -I {} bash -c \
    "process_chromosome {} '$BAM' '$REFERENCE' '$TEMP_DIR'"

echo -e "  ${GREEN}✓${NC} Wszystkie chromosomy przetworzone"
echo "  Zakończono: $(date '+%H:%M:%S')"

#===============================================================================
# KROK 3: Łączenie VCF
#===============================================================================

echo -e "\n${YELLOW}[3/5] Łączenie plików VCF...${NC}"

# Lista plików do połączenia (w kolejności chromosomów)
VCF_LIST="${TEMP_DIR}/vcf_list.txt"
> "$VCF_LIST"

for chrom in $CHROMOSOMES; do
    vcf_file="${TEMP_DIR}/${chrom}.vcf.gz"
    if [ -f "$vcf_file" ]; then
        echo "$vcf_file" >> "$VCF_LIST"
    fi
done

# Połącz wszystkie VCF
bcftools concat \
    -f "$VCF_LIST" \
    -Oz \
    -o "$OUTPUT" \
    --threads $THREADS

echo -e "  ${GREEN}✓${NC} VCF połączony"

#===============================================================================
# KROK 4: Indeksowanie końcowe
#===============================================================================

echo -e "\n${YELLOW}[4/5] Indeksowanie końcowego VCF...${NC}"
tabix -p vcf "$OUTPUT"
echo -e "  ${GREEN}✓${NC} Indeks utworzony"

#===============================================================================
# KROK 5: Statystyki
#===============================================================================

echo -e "\n${YELLOW}[5/5] Generowanie statystyk...${NC}"

STATS_FILE="${OUTPUT%.vcf.gz}_stats.txt"
bcftools stats "$OUTPUT" > "$STATS_FILE"

# Szybkie podsumowanie
TOTAL=$(bcftools view -H "$OUTPUT" | wc -l | tr -d ' ')
SNPS=$(bcftools view -v snps -H "$OUTPUT" | wc -l | tr -d ' ')
INDELS=$(bcftools view -v indels -H "$OUTPUT" | wc -l | tr -d ' ')

echo -e "${GREEN}"
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║                         GOTOWE!                                   ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo "Pliki wynikowe:"
echo "  • $OUTPUT"
echo "  • ${OUTPUT}.tbi"
echo "  • $STATS_FILE"
echo ""
echo "Statystyki:"
echo "  • Wszystkie warianty: $TOTAL"
echo "  • SNPs:               $SNPS"
echo "  • Indele:             $INDELS"
echo ""
echo "Podgląd pierwszych wariantów:"
bcftools view -H "$OUTPUT" | head -3
