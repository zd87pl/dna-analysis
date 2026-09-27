#!/bin/bash
#===============================================================================
# ULTRA VCF - Maximum Speed Variant Calling
# Wszystkie optymalizacje dla M4 MacBook Air 24GB
#===============================================================================

set -e

BAM="${1:-}"
OUTPUT="${2:-variants.vcf.gz}"
REFERENCE="${3:-human_g1k_v37.fasta}"

# M4 Air specs
TOTAL_CORES=10      # 4 Performance + 6 Efficiency
PARALLEL_JOBS=6     # Zostaw trochę dla systemu
MAX_DEPTH=150       # Mniej = szybciej (standardowo 250)

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║              ULTRA VCF - Maximum Speed Mode                       ║"
echo "║              M4 MacBook Air 24GB Optimized                        ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

if [ -z "$BAM" ] || [ -z "$REFERENCE" ]; then
    echo "Użycie: $0 <plik.bam> [output.vcf.gz] [reference.fa]" >&2
    exit 1
fi
if [ ! -f "$BAM" ]; then
    echo "BŁĄD: Brak pliku BAM: $BAM" >&2
    exit 1
fi
if [ ! -f "$REFERENCE" ]; then
    echo "BŁĄD: Brak pliku referencyjnego: $REFERENCE" >&2
    exit 1
fi
if [ ! -f "${REFERENCE}.fai" ]; then
    samtools faidx "$REFERENCE"
fi

TEMP_DIR=$(mktemp -d)
# Pojedyncze cudzysłowy: $TEMP_DIR rozwijany przy wyjściu, nie teraz
trap 'rm -rf "$TEMP_DIR"' EXIT

START_TIME=$(date +%s)

#===============================================================================
# Optymalizacja 1: RAM disk dla plików tymczasowych (jeśli dostępny)
#===============================================================================

# Na macOS można użyć /tmp który jest często na SSD
# Dla ekstra szybkości możesz stworzyć RAM disk:
# diskutil erasevolume HFS+ 'RAMDisk' `hdiutil attach -nomount ram://8388608`
# (4GB RAM disk)

if [ -d "/Volumes/RAMDisk" ]; then
    rm -rf "$TEMP_DIR"
    TEMP_DIR="/Volumes/RAMDisk/vcf_temp_$$"
    mkdir -p "$TEMP_DIR"
    trap 'rm -rf "$TEMP_DIR"' EXIT
    echo "🚀 Używam RAM disk dla max szybkości!"
fi

#===============================================================================
# Optymalizacja 2: Podział na regiony (nie chromosomy)
#===============================================================================

echo ""
echo "[1/4] Dzielenie genomu na regiony..."

# Podziel genom na ~50 regionów dla lepszego load balancing
# Chromosom 1 jest ~8x większy niż chr21, więc podział na równe części jest lepszy

generate_regions() {
    local ref=$1
    local num_regions=$2
    local output=$3
    
    # Pobierz długości chromosomów
    awk -v n="$num_regions" '
    BEGIN { total = 0 }
    { 
        chroms[NR] = $1
        lens[NR] = $2
        total += $2
    }
    END {
        region_size = int(total / n)
        region_id = 0
        
        for (i = 1; i <= NR; i++) {
            chrom = chroms[i]
            len = lens[i]
            start = 1
            
            while (start < len) {
                end = start + region_size - 1
                if (end > len) end = len
                
                printf "%s:%d-%d\n", chrom, start, end
                region_id++
                start = end + 1
            }
        }
    }
    ' "${ref}.fai" > "$output"
}

REGIONS_FILE="$TEMP_DIR/regions.txt"
generate_regions "$REFERENCE" 50 "$REGIONS_FILE"
NUM_REGIONS=$(wc -l < "$REGIONS_FILE" | tr -d ' ')
echo "  Utworzono $NUM_REGIONS regionów"

#===============================================================================
# Optymalizacja 3: Parallel z progress bar
#===============================================================================

echo ""
echo "[2/4] Variant calling ($PARALLEL_JOBS równoległych procesów)..."
echo "  Start: $(date '+%H:%M:%S')"

# Funkcja przetwarzania regionu
process_region() {
    local region=$1
    local idx=$2
    local bam=$3
    local ref=$4
    local temp_dir=$5
    
    # Bezpieczna nazwa pliku
    local safe_name
    safe_name=$(echo "$region" | tr ':' '_' | tr '-' '_')
    local output_vcf="${temp_dir}/${idx}_${safe_name}.vcf.gz"

    bcftools mpileup \
        -Ou \
        -f "$ref" \
        -r "$region" \
        --max-depth "$MAX_DEPTH" \
        --min-MQ 20 \
        --min-BQ 20 \
        --no-BAQ \
        "$bam" | \
    bcftools call \
        -mv \
        -Oz \
        -o "$output_vcf"
    
    # Indeksuj tylko jeśli plik istnieje i ma zawartość
    if [ -s "$output_vcf" ]; then
        tabix -f -p vcf "$output_vcf" 2>/dev/null || true
    fi
}

export -f process_region
export MAX_DEPTH

# Progress tracking
PROGRESS_FILE="$TEMP_DIR/progress"
echo 0 > "$PROGRESS_FILE"

# Uruchom równolegle z numeracją (dla poprawnej kolejności).
# Wszystkie wartości (ścieżki, region) idą jako argumenty pozycyjne, nigdy
# wklejane do łańcucha poleceń: -L 1 dokleja "idx region" po stałych
# argumentach, więc $1..$5 = BAM REF TEMP_DIR PROGRESS_FILE NUM_REGIONS,
# $6 = idx, $7 = region.
if ! cat -n "$REGIONS_FILE" | \
    xargs -P "$PARALLEL_JOBS" -L 1 bash -c '
        set -euo pipefail
        bam=$1; ref=$2; temp_dir=$3; progress_file=$4; num_regions=$5
        idx=$6
        region=$7
        process_region "$region" "$idx" "$bam" "$ref" "$temp_dir"

        # Update progress
        current=$(cat "$progress_file")
        echo $((current + 1)) > "$progress_file"
        printf "\r  Postęp: %d/%d regionów" $((current + 1)) "$num_regions"
    ' _ "$BAM" "$REFERENCE" "$TEMP_DIR" "$PROGRESS_FILE" "$NUM_REGIONS"; then
    echo "" >&2
    echo "BŁĄD: Variant calling nie powiódł się dla co najmniej jednego regionu" >&2
    exit 1
fi

echo ""
echo "  Zakończono: $(date '+%H:%M:%S')"

#===============================================================================
# Optymalizacja 4: Szybkie łączenie
#===============================================================================

echo ""
echo "[3/4] Łączenie VCF (concat)..."

# Lista plików posortowana po numerze
VCF_LIST="$TEMP_DIR/vcf_list.txt"
ls "$TEMP_DIR"/*.vcf.gz 2>/dev/null | sort -t'_' -k1 -n > "$VCF_LIST" || true

if [ ! -s "$VCF_LIST" ]; then
    echo "BŁĄD: Brak plików VCF do połączenia!"
    exit 1
fi

# Połącz
bcftools concat \
    -f "$VCF_LIST" \
    -Oz \
    -o "$OUTPUT" \
    --threads "$TOTAL_CORES" \
    -a  # Allow overlaps (dla regionów które mogą się nakładać)

#===============================================================================
# Finalizacja
#===============================================================================

echo ""
echo "[4/4] Finalizacja..."

tabix -p vcf "$OUTPUT"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
MINUTES=$((DURATION / 60))
SECONDS=$((DURATION % 60))

# Statystyki
TOTAL=$(bcftools view -H "$OUTPUT" 2>/dev/null | wc -l | tr -d ' ')
SNPS=$(bcftools view -v snps -H "$OUTPUT" 2>/dev/null | wc -l | tr -d ' ')
INDELS=$(bcftools view -v indels -H "$OUTPUT" 2>/dev/null | wc -l | tr -d ' ')

echo ""
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║                      ZAKOŃCZONO!                                  ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "⏱️  Czas wykonania: ${MINUTES}m ${SECONDS}s"
echo ""
echo "📊 Statystyki:"
echo "   • Wszystkie warianty: $TOTAL"
echo "   • SNPs:               $SNPS"
echo "   • Indele:             $INDELS"
echo ""
echo "📁 Pliki:"
echo "   • $OUTPUT"
echo "   • ${OUTPUT}.tbi"
