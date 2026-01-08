#!/bin/bash
#===============================================================================
# QUICK VCF GENERATION FROM BAM
# Minimal version - VCF output only
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

SCRIPT_NAME="$(basename "$0")"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

error_exit() {
    echo -e "${RED}Error: $1${NC}" >&2
    exit 1
}

warn() {
    echo -e "${YELLOW}Warning: $1${NC}" >&2
}

success() {
    echo -e "${GREEN}$1${NC}"
}

# Check for required dependencies
check_dependencies() {
    local missing=()
    for cmd in bcftools samtools tabix; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        error_exit "Missing required dependencies: ${missing[*]}\nPlease run: ./install.sh"
    fi
}

# Validate BAM file
validate_bam() {
    local bam="$1"

    if [[ ! -f "$bam" ]]; then
        error_exit "BAM file not found: $bam"
    fi

    # Use samtools to validate
    if ! samtools quickcheck "$bam" 2>/dev/null; then
        error_exit "File does not appear to be a valid BAM: $bam"
    fi
}

# Validate reference file
validate_reference() {
    local ref="$1"

    if [[ ! -f "$ref" ]]; then
        error_exit "Reference file not found: $ref"
    fi

    # Check for index
    if [[ ! -f "${ref}.fai" ]]; then
        warn "Reference index not found, creating..."
        if ! samtools faidx "$ref"; then
            error_exit "Failed to index reference file"
        fi
    fi
}

usage() {
    echo "Usage: $SCRIPT_NAME <bam_file> [output.vcf.gz] [reference.fa] [threads]"
    echo ""
    echo "Arguments:"
    echo "  bam_file      Path to BAM file (required)"
    echo "  output        Output VCF filename (default: <bam_name>_variants.vcf.gz)"
    echo "  reference     Reference genome FASTA (required for variant calling)"
    echo "  threads       Number of threads (default: 8)"
    echo ""
    echo "Examples:"
    echo "  $SCRIPT_NAME sample.bam"
    echo "  $SCRIPT_NAME sample.bam my_variants.vcf.gz GRCh38.fa 16"
    exit 1
}

#===============================================================================
# MAIN SCRIPT
#===============================================================================

# Check arguments
if [[ $# -lt 1 ]]; then
    usage
fi

BAM="$1"
THREADS="${4:-8}"

# Validate threads is a number
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
    THREADS=8
fi

# Check dependencies
check_dependencies

# Validate BAM
validate_bam "$BAM"

# Set output name
if [[ -n "${2:-}" ]]; then
    OUTPUT="$2"
else
    OUTPUT="$(basename "$BAM" .bam)_variants.vcf.gz"
fi

# Handle reference
REFERENCE="${3:-}"

echo "═══════════════════════════════════════════════════════════════"
echo "         QUICK VCF GENERATION"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "BAM:        $BAM"
echo "Output:     $OUTPUT"
echo "Reference:  ${REFERENCE:-NOT PROVIDED (limited functionality)}"
echo "Threads:    $THREADS"
echo ""

# Check if BAM is indexed
if [[ ! -f "${BAM}.bai" ]] && [[ ! -f "${BAM%.*}.bai" ]]; then
    echo "[1/4] Indexing BAM..."
    if ! samtools index -@ "$THREADS" "$BAM"; then
        error_exit "Failed to index BAM file"
    fi
else
    success "[1/4] BAM index exists ✓"
fi

# Variant calling
echo "[2/4] Variant calling (this may take 2-6h for WGS)..."
echo "      Started: $(date)"

if [[ -n "$REFERENCE" ]] && [[ -f "$REFERENCE" ]]; then
    # Validate reference
    validate_reference "$REFERENCE"

    # Full variant calling with reference
    if ! bcftools mpileup \
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
        -o "$OUTPUT"; then
        error_exit "Variant calling failed"
    fi
else
    warn "No reference provided - using simplified mode"
    echo ""
    echo "To download a reference genome:"
    echo "  wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
    echo "  gunzip human_g1k_v37.fasta.gz"
    echo "  samtools faidx human_g1k_v37.fasta"
    echo ""

    # Without reference - just generate coverage stats
    echo "Generating coverage statistics instead of VCF..."
    COVERAGE_FILE="${OUTPUT%.vcf.gz}_coverage_per_chrom.txt"

    if ! samtools depth "$BAM" | \
        awk 'BEGIN{OFS="\t"} {sum[$1]+=$3; cnt[$1]++} END{for(c in sum) print c, sum[c]/cnt[c]}' \
        > "$COVERAGE_FILE"; then
        error_exit "Coverage calculation failed"
    fi

    success "Saved: $COVERAGE_FILE"
    exit 0
fi

echo "      Finished: $(date)"

# Index VCF
echo "[3/4] Indexing VCF..."
if ! tabix -p vcf "$OUTPUT"; then
    warn "Failed to create tabix index, trying bcftools..."
    bcftools index "$OUTPUT" || warn "Indexing failed"
fi

# Statistics
echo "[4/4] Generating statistics..."
STATS_FILE="${OUTPUT%.vcf.gz}_stats.txt"
bcftools stats "$OUTPUT" > "$STATS_FILE" 2>/dev/null || true

# Summary
echo ""
echo "═══════════════════════════════════════════════════════════════"
echo "                      DONE!"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "Generated files:"
echo "  • $OUTPUT"
if [[ -f "${OUTPUT}.tbi" ]]; then
    echo "  • ${OUTPUT}.tbi"
fi
if [[ -f "${OUTPUT}.csi" ]]; then
    echo "  • ${OUTPUT}.csi"
fi
if [[ -f "$STATS_FILE" ]]; then
    echo "  • $STATS_FILE"
fi
echo ""

# Quick statistics
TOTAL=$(bcftools view -H "$OUTPUT" 2>/dev/null | wc -l || echo "0")
SNPS=$(bcftools view -v snps -H "$OUTPUT" 2>/dev/null | wc -l || echo "0")
INDELS=$(bcftools view -v indels -H "$OUTPUT" 2>/dev/null | wc -l || echo "0")

echo "Variant statistics:"
echo "  • Total:   $TOTAL"
echo "  • SNPs:    $SNPS"
echo "  • Indels:  $INDELS"
echo ""
echo "Preview of first variants:"
bcftools view -H "$OUTPUT" 2>/dev/null | head -5 || true
echo ""
