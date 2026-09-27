#!/bin/bash
#===============================================================================
# TURBO VCF - Parallel Variant Calling dla Apple Silicon
# Zoptymalizowany dla M4 MacBook Air z 24GB RAM
#===============================================================================

set -e

BAM="${1:-}"
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

# Sprawdź argumenty (BAM jest wymagany)
usage() {
    echo "Użycie: $0 <plik.bam> [output.vcf.gz] [reference.fa] [threads]" >&2
    echo "" >&2
    echo "Przykład:" >&2
    echo "  $0 sample.bam sample_variants.vcf.gz GRCh38.fa 10" >&2
    exit 1
}
[ -n "$BAM" ] || usage
if [ ! -f "$BAM" ]; then
    echo -e "${RED}✗${NC} Brak pliku BAM: $BAM" >&2
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
    samtools index -@ "$THREADS" "$BAM"
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

# Pobierz listę chromosomów. Nazwy kontigów pochodzą z nagłówka BAM (dane
# wejściowe!): dopuszczamy tylko bezpieczne znaki i nigdy nie wklejamy ich do
# łańcucha poleceń (patrz xargs poniżej).
ALL_CONTIGS=$(samtools idxstats "$BAM" | cut -f1 | grep -v '^\*$' || true)
BAD_CONTIGS=$(echo "$ALL_CONTIGS" | grep -v '^$' | grep -Ev '^[0-9A-Za-z_.-]+$' || true)
if [ -n "$BAD_CONTIGS" ]; then
    echo -e "  ${RED}✗${NC} Niedozwolone nazwy kontigów w nagłówku BAM (dozwolone: [0-9A-Za-z_.-]):" >&2
    echo "$BAD_CONTIGS" | sed 's/^/      /' >&2
    exit 1
fi
CHROMOSOMES=$(echo "$ALL_CONTIGS" | head -25)
# grep -c prints 0 but exits 1 on empty input; without `|| true` set -e would
# abort here before the diagnostic below.
NUM_CHROM=$(echo "$CHROMOSOMES" | grep -c . || true)
if [ "$NUM_CHROM" -eq 0 ]; then
    echo -e "  ${RED}✗${NC} Nie znaleziono kontigów w nagłówku BAM" >&2
    exit 1
fi
echo -e "  ${GREEN}✓${NC} Znaleziono $NUM_CHROM chromosomów"

# Każdy kontig BAM musi istnieć w referencji: bcftools mpileup kończy się
# kodem 0 i pustym VCF dla kontigu, którego referencja nie zawiera.
MISSING_CONTIGS=$(echo "$CHROMOSOMES" | grep -vxF -f <(cut -f1 "${REFERENCE}.fai") || true)
if [ -n "$MISSING_CONTIGS" ]; then
    echo -e "  ${RED}✗${NC} Kontigi BAM nieobecne w referencji ${REFERENCE}:" >&2
    echo "$MISSING_CONTIGS" | sed 's/^/      /' >&2
    echo "  Sprawdź prefiks 'chr' i wersję genomu (BAM i referencja muszą się zgadzać)." >&2
    exit 1
fi

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

    # Variant calling dla chromosomu. Błędy nie są wyciszane, a każdy
    # niezerowy status (pipefail) przerywa proces potomny, co xargs propaguje.
    bcftools mpileup \
        -Ou \
        -f "$ref" \
        -r "$chrom" \
        --max-depth 200 \
        --min-MQ 20 \
        --min-BQ 20 \
        "$bam" | \
    bcftools call \
        -mv \
        -Oz \
        -o "$output_vcf"

    # Indeksuj
    tabix -p vcf "$output_vcf"

    echo "  ✓ $chrom"
}

export -f process_chromosome

# Uruchom równolegle (max 8 procesów - balans CPU/RAM)
# Dla M4 Air z 24GB: 8 procesów * ~2GB = 16GB używane
MAX_PARALLEL=8

# Nazwa kontigu i ścieżki są przekazywane jako argumenty pozycyjne, nigdy
# wklejane do łańcucha poleceń (kontig "chr2;cmd" nie wykona cmd).
if ! echo "$CHROMOSOMES" | xargs -P "$MAX_PARALLEL" -I {} bash -c \
    'set -euo pipefail; process_chromosome "$1" "$2" "$3" "$4"' _ {} "$BAM" "$REFERENCE" "$TEMP_DIR"; then
    echo -e "  ${RED}✗${NC} Variant calling nie powiódł się dla co najmniej jednego chromosomu" >&2
    exit 1
fi

echo -e "  ${GREEN}✓${NC} Wszystkie chromosomy przetworzone"
echo "  Zakończono: $(date '+%H:%M:%S')"

#===============================================================================
# KROK 3: Łączenie VCF
#===============================================================================

echo -e "\n${YELLOW}[3/5] Łączenie plików VCF...${NC}"

# Lista plików do połączenia (w kolejności chromosomów)
VCF_LIST="${TEMP_DIR}/vcf_list.txt"
: > "$VCF_LIST"

# Brak pliku (lub indeksu) dla chromosomu to błąd, nie cichy brak danych.
for chrom in $CHROMOSOMES; do
    vcf_file="${TEMP_DIR}/${chrom}.vcf.gz"
    if [ ! -s "$vcf_file" ] || [ ! -s "${vcf_file}.tbi" ]; then
        echo -e "  ${RED}✗${NC} Brak wyniku dla chromosomu $chrom ($vcf_file)" >&2
        exit 1
    fi
    echo "$vcf_file" >> "$VCF_LIST"
done

# Połącz wszystkie VCF
bcftools concat \
    -f "$VCF_LIST" \
    -Oz \
    -o "$OUTPUT" \
    --threads "$THREADS"

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
