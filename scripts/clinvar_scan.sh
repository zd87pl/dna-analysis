#!/bin/bash
#===============================================================================
# CLINVAR PATHOGENIC VARIANT SCANNER
# Finds all pathogenic variants in your genome
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

SCRIPT_NAME="$(basename "$0")"
OUTDIR="clinvar_analysis"
LOCK_FILE="/tmp/clinvar_scan.lock"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_INDEX_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"

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
    cleanup
    exit 1
}

warn() {
    echo -e "${YELLOW}Warning: $1${NC}" >&2
}

success() {
    echo -e "${GREEN}$1${NC}"
}

cleanup() {
    # Remove lock file if we own it
    if [[ -f "$LOCK_FILE" ]] && [[ "$(cat "$LOCK_FILE" 2>/dev/null)" == "$$" ]]; then
        rm -f "$LOCK_FILE"
    fi
}

trap cleanup EXIT

# Check for required dependencies
check_dependencies() {
    local missing=()
    for cmd in bcftools wget; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        error_exit "Missing required dependencies: ${missing[*]}\nPlease run: ./install.sh"
    fi
}

# Acquire lock to prevent race conditions
acquire_lock() {
    local max_wait=60
    local waited=0

    while [[ -f "$LOCK_FILE" ]]; do
        local lock_pid
        lock_pid=$(cat "$LOCK_FILE" 2>/dev/null || echo "")

        # Check if the process holding the lock is still running
        if [[ -n "$lock_pid" ]] && ! kill -0 "$lock_pid" 2>/dev/null; then
            # Stale lock, remove it
            rm -f "$LOCK_FILE"
            break
        fi

        if [[ $waited -ge $max_wait ]]; then
            error_exit "Could not acquire lock after ${max_wait}s. Another instance may be running."
        fi

        echo "Waiting for another instance to finish..."
        sleep 2
        ((waited+=2))
    done

    echo "$$" > "$LOCK_FILE"
}

# Validate VCF file
validate_vcf() {
    local vcf="$1"

    if [[ ! -f "$vcf" ]]; then
        error_exit "VCF file not found: $vcf"
    fi

    # Basic format check
    if [[ "$vcf" =~ \.gz$ ]]; then
        if ! zcat "$vcf" 2>/dev/null | head -1 | grep -q "^##fileformat=VCF"; then
            error_exit "File does not appear to be a valid VCF: $vcf"
        fi
    else
        if ! head -1 "$vcf" 2>/dev/null | grep -q "^##fileformat=VCF"; then
            error_exit "File does not appear to be a valid VCF: $vcf"
        fi
    fi
}

#===============================================================================
# MAIN SCRIPT
#===============================================================================

# Check arguments
if [[ $# -lt 1 ]]; then
    echo "Usage: $SCRIPT_NAME <vcf_file>"
    echo ""
    echo "Arguments:"
    echo "  vcf_file    Path to your VCF file (.vcf or .vcf.gz)"
    echo ""
    echo "Example:"
    echo "  $SCRIPT_NAME my_variants.vcf.gz"
    exit 1
fi

VCF="$1"

# Check dependencies
check_dependencies

# Validate input VCF
validate_vcf "$VCF"

# Create output directory
mkdir -p "$OUTDIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║            CLINVAR PATHOGENIC VARIANT SCANNER                             ║"
echo "║            Searching for pathogenic variants in your genome               ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "VCF File: $VCF"
echo "Results: $OUTDIR/"
echo ""

#===============================================================================
# STEP 1: Download ClinVar (if not present)
#===============================================================================

CLINVAR_VCF="$OUTDIR/clinvar_GRCh38.vcf.gz"

# Acquire lock before downloading to prevent race conditions
acquire_lock

if [[ ! -f "$CLINVAR_VCF" ]]; then
    echo "📥 Downloading ClinVar database (GRCh38)..."
    echo "   This may take a few minutes (~70MB)..."

    # Download with error checking
    if ! wget -q --show-progress -O "${CLINVAR_VCF}.tmp" "$CLINVAR_URL"; then
        rm -f "${CLINVAR_VCF}.tmp"
        error_exit "Failed to download ClinVar VCF"
    fi

    if ! wget -q -O "${CLINVAR_VCF}.tbi.tmp" "$CLINVAR_INDEX_URL"; then
        rm -f "${CLINVAR_VCF}.tmp" "${CLINVAR_VCF}.tbi.tmp"
        error_exit "Failed to download ClinVar index"
    fi

    # Verify downloaded files are valid
    if ! zcat "${CLINVAR_VCF}.tmp" 2>/dev/null | head -1 | grep -q "^##fileformat=VCF"; then
        rm -f "${CLINVAR_VCF}.tmp" "${CLINVAR_VCF}.tbi.tmp"
        error_exit "Downloaded ClinVar file appears to be invalid"
    fi

    # Move to final location atomically
    mv "${CLINVAR_VCF}.tmp" "$CLINVAR_VCF"
    mv "${CLINVAR_VCF}.tbi.tmp" "${CLINVAR_VCF}.tbi"

    success "✅ ClinVar downloaded successfully!"
else
    echo "✅ Using existing ClinVar database"
fi
echo ""

#===============================================================================
# STEP 2: Extract pathogenic variants from ClinVar
#===============================================================================

echo "🔬 Extracting Pathogenic/Likely_pathogenic variants from ClinVar..."

CLINVAR_PATH_VCF="$OUTDIR/clinvar_pathogenic_only.vcf.gz"

if [[ ! -f "$CLINVAR_PATH_VCF" ]]; then
    if ! bcftools view -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' \
        "$CLINVAR_VCF" -Oz -o "${CLINVAR_PATH_VCF}.tmp" 2>/dev/null; then
        rm -f "${CLINVAR_PATH_VCF}.tmp"
        error_exit "Failed to filter ClinVar for pathogenic variants"
    fi

    if ! bcftools index -f "${CLINVAR_PATH_VCF}.tmp" 2>/dev/null; then
        rm -f "${CLINVAR_PATH_VCF}.tmp"
        error_exit "Failed to index filtered ClinVar file"
    fi

    mv "${CLINVAR_PATH_VCF}.tmp" "$CLINVAR_PATH_VCF"
    mv "${CLINVAR_PATH_VCF}.tmp.csi" "${CLINVAR_PATH_VCF}.csi" 2>/dev/null || true

    PATHOGENIC_COUNT=$(bcftools view -H "$CLINVAR_PATH_VCF" 2>/dev/null | wc -l)
    echo "   ClinVar contains $PATHOGENIC_COUNT Pathogenic/Likely_pathogenic variants"
else
    echo "   Using existing pathogenic variants file"
fi
echo ""

#===============================================================================
# STEP 3: Find common variants (bcftools isec)
#===============================================================================

echo "🔍 Comparing your genome with ClinVar Pathogenic..."
echo ""

mkdir -p "$OUTDIR/isec_results"

# Ensure VCF is indexed
if [[ ! -f "${VCF}.tbi" ]] && [[ ! -f "${VCF}.csi" ]]; then
    echo "   Indexing your VCF..."
    if ! bcftools index -f "$VCF" 2>/dev/null; then
        error_exit "Failed to index VCF file. Ensure it's bgzip compressed."
    fi
fi

# Find intersection - variants that are BOTH in your VCF AND in ClinVar Pathogenic
if ! bcftools isec -p "$OUTDIR/isec_results" -n=2 -w1 \
    "$VCF" "$CLINVAR_PATH_VCF" 2>/dev/null; then
    error_exit "Failed to compare VCF with ClinVar"
fi

#===============================================================================
# STEP 4: Analyze results
#===============================================================================

RESULTS_FILE="$OUTDIR/isec_results/0000.vcf"
REPORT="$OUTDIR/PATHOGENIC_VARIANTS_FOUND.txt"

if [[ -f "$RESULTS_FILE" ]]; then
    FOUND_COUNT=$(grep -cv "^#" "$RESULTS_FILE" 2>/dev/null || echo "0")

    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    if [[ $FOUND_COUNT -eq 0 ]]; then
        success "✅ No pathogenic variants found in your genome!"
        echo ""
        echo "   This is GOOD news - none of the known pathogenic variants"
        echo "   from ClinVar were detected in your VCF."

        cat << EOF > "$REPORT"
================================================================================
                CLINVAR REPORT - PATHOGENIC VARIANTS
                Generated: $(date)
================================================================================

RESULT: NO pathogenic variants found

Your genome was compared against the ClinVar database containing all known
variants classified as Pathogenic or Likely_pathogenic.

No matches were found.

NOTE: This does not mean absence of ANY genetic risk.
- ClinVar does not contain all pathogenic variants
- Some variants may not yet be classified
- This analysis does not include CNV, STR, or structural variants

================================================================================
EOF
    else
        echo -e "${RED}🔴 FOUND $FOUND_COUNT variants present in ClinVar Pathogenic!${NC}"
        echo ""

        cat << EOF > "$REPORT"
================================================================================
                CLINVAR REPORT - PATHOGENIC VARIANTS
                Generated: $(date)
================================================================================

⚠️ FOUND $FOUND_COUNT VARIANTS IN CLINVAR PATHOGENIC

DETAILS BELOW:
================================================================================

EOF

        echo "   Details of found variants:"
        echo "   ─────────────────────────────────────────"

        # Process up to 50 variants
        local count=0
        while IFS=$'\t' read -r CHROM POS RSID REF ALT QUAL FILTER INFO FORMAT SAMPLE _; do
            [[ "$CHROM" == "#"* ]] && continue
            ((count++))
            [[ $count -gt 50 ]] && break

            GT="${SAMPLE%%:*}"

            # Get info from ClinVar
            CLINVAR_INFO=$(bcftools query -r "${CHROM}:${POS}" \
                -f '%ID\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\n' \
                "$CLINVAR_PATH_VCF" 2>/dev/null | head -1) || true

            CV_ID=$(echo "$CLINVAR_INFO" | cut -f1)
            CLNSIG=$(echo "$CLINVAR_INFO" | cut -f2)
            CLNDN=$(echo "$CLINVAR_INFO" | cut -f3 | sed 's/_/ /g')
            REVIEW=$(echo "$CLINVAR_INFO" | cut -f4)

            # Genotype interpretation
            case $GT in
                "0/1"|"0|1"|"1/0"|"1|0") GENO="HETEROZYGOUS" ;;
                "1/1"|"1|1") GENO="HOMOZYGOUS" ;;
                *) GENO="$GT" ;;
            esac

            echo ""
            echo "   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "   📍 chr${CHROM}:${POS} ${REF}>${ALT}"
            echo "      ClinVar ID: ${CV_ID:-N/A}"
            echo "      Classification: ${CLNSIG:-N/A}"
            echo "      Disease: ${CLNDN:-N/A}"
            echo "      Your genotype: $GENO"
            echo "      Review status: ${REVIEW:-N/A}"

            # Save to report
            cat << EOF >> "$REPORT"
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
VARIANT: chr${CHROM}:${POS} ${REF}>${ALT}
ClinVar ID: ${CV_ID:-N/A}
Classification: ${CLNSIG:-N/A}
Disease: ${CLNDN:-N/A}
Your genotype: $GENO
Review status: ${REVIEW:-N/A}
Link: https://www.ncbi.nlm.nih.gov/clinvar/?term=${CV_ID:-}

EOF
        done < "$RESULTS_FILE"

        if [[ $FOUND_COUNT -gt 50 ]]; then
            echo ""
            echo "   ... and $((FOUND_COUNT - 50)) more (see full report)"
        fi
    fi

    echo ""
else
    error_exit "Comparison failed - check if VCF is valid"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📁 Results saved:"
echo "   $REPORT"
echo "   $OUTDIR/isec_results/"
echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  ⚠️ NOTE: Results require interpretation by a clinical geneticist!        ║"
echo "║  Not all 'pathogenic' variants cause disease in everyone.                 ║"
echo "║  Penetrance, expression, and environmental factors matter.                ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
