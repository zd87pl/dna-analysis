#!/bin/bash
#===============================================================================
# CLINVAR PATHOGENIC VARIANT SCANNER
# Finds all pathogenic variants in your genome
#===============================================================================

set -euo pipefail

# Shared VCF helpers (hx_prepare_vcf, hx_validate_vcf, hx_report_coverage, ...)
# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"

#===============================================================================
# CONFIGURATION
#===============================================================================

SCRIPT_NAME="$(basename "$0")"
OUTDIR="clinvar_analysis"
# ClinVar cache: shared across runs (the frontend sets HELIXIGHT_REFERENCE_DIR;
# otherwise the ~70 MB download would be repeated in every per-run directory).
REF_DIR="${HELIXIGHT_REFERENCE_DIR:-$OUTDIR}"
CLINVAR_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
# Re-download when the cached copy is older than this many days.
CLINVAR_MAX_AGE_DAYS="${HELIXIGHT_CLINVAR_MAX_AGE_DAYS:-90}"
# Temporary download file (set while a download is in flight; removed by cleanup)
DOWNLOAD_TMP=""

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

# md5 hex digest of a file: md5sum (coreutils) or, on macOS, md5 -q.
file_md5() {
    if command -v md5sum &> /dev/null; then
        md5sum "$1" | cut -d' ' -f1
    else
        md5 -q "$1"
    fi
}

cleanup() {
    # Remove a partial download; the flock on fd 9 is released by the kernel
    # when the process exits, so no stale-lock handling is needed.
    if [[ -n "$DOWNLOAD_TMP" ]]; then
        rm -f "$DOWNLOAD_TMP" "$DOWNLOAD_TMP.md5"
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
    # flock (util-linux) and md5sum (coreutils) are not shipped by macOS:
    # md5sum falls back to `md5 -q` (see file_md5), a missing flock only
    # disables the download lock (warned below).
    if ! command -v md5sum &> /dev/null && ! command -v md5 &> /dev/null; then
        missing+=("md5sum")
    fi

    if [[ ${#missing[@]} -gt 0 ]]; then
        error_exit "Missing required dependencies: ${missing[*]}\nPlease run: ./install.sh"
    fi
}

# Validate VCF file (plain, gzip or BGZF; header must parse)
validate_vcf() {
    local vcf="$1"

    if [[ ! -f "$vcf" ]]; then
        error_exit "VCF file not found: $vcf"
    fi

    if ! hx_validate_vcf "$vcf"; then
        error_exit "File does not appear to be a valid VCF: $vcf"
    fi
}

# ##fileDate of a ClinVar VCF (provenance), "unknown" when absent
clinvar_release_date() {
    local line
    line=$(bcftools view -h "$1" 2>/dev/null | grep -m1 '^##fileDate' || true)
    line=${line#*=}
    printf '%s\n' "${line:-unknown}"
}

# True when FILE is missing or older than CLINVAR_MAX_AGE_DAYS days
clinvar_cache_stale() {
    local f="$1"
    [[ -s "$f" ]] || return 0
    [[ -n $(find "$f" -mtime +"$CLINVAR_MAX_AGE_DAYS" -print 2>/dev/null) ]]
}

# Download ClinVar into REF_DIR with md5 verification. Returns non-zero on
# failure (the caller decides whether a stale cached copy can be used).
download_clinvar() {
    local tmp md5_expected md5_actual

    echo "📥 Downloading ClinVar database (GRCh38)..."
    echo "   This may take a few minutes (~70MB)..."

    tmp=$(mktemp "$REF_DIR/clinvar_download.XXXXXX") || return 1
    DOWNLOAD_TMP="$tmp"

    # wget's own diagnostics stay on stderr so a blocked network is visible.
    if ! wget -nv --timeout=60 --tries=3 -c "$CLINVAR_URL" -O "$tmp"; then
        warn "Failed to download ClinVar VCF from $CLINVAR_URL"
        rm -f "$tmp"; DOWNLOAD_TMP=""
        return 1
    fi

    # Integrity: NCBI publishes <file>.md5 next to the VCF. Compare the first
    # 32-hex token in it with the local md5sum (works for both the
    # "MD5(file)= hash" and "hash  file" layouts).
    if wget -nv --timeout=60 --tries=3 "${CLINVAR_URL}.md5" -O "$tmp.md5" 2>/dev/null; then
        # (tr, not ${var,,}: bash 4 only, macOS ships bash 3.2)
        md5_expected=$(grep -oE '[0-9a-fA-F]{32}' "$tmp.md5" | head -n1 | tr '[:upper:]' '[:lower:]' || true)
        md5_actual=$(file_md5 "$tmp")
        if [[ -z "$md5_expected" ]] || [[ "$md5_expected" != "$md5_actual" ]]; then
            warn "ClinVar checksum mismatch (expected ${md5_expected:-?}, got $md5_actual) - discarding download"
            rm -f "$tmp" "$tmp.md5"; DOWNLOAD_TMP=""
            return 1
        fi
        echo "   md5 checksum verified"
    elif [[ "${HELIXIGHT_SKIP_CHECKSUM:-0}" == "1" ]]; then
        warn "ClinVar .md5 file unavailable - HELIXIGHT_SKIP_CHECKSUM=1 set, continuing WITHOUT integrity verification"
    else
        warn "ClinVar .md5 file unavailable at ${CLINVAR_URL}.md5 - refusing unverified download (set HELIXIGHT_SKIP_CHECKSUM=1 to bypass)"
        rm -f "$tmp" "$tmp.md5"; DOWNLOAD_TMP=""
        return 1
    fi

    # Verify the downloaded file is a readable VCF
    if ! hx_validate_vcf "$tmp"; then
        warn "Downloaded ClinVar file appears to be invalid"
        rm -f "$tmp" "$tmp.md5"; DOWNLOAD_TMP=""
        return 1
    fi

    # Move to final location atomically; drop the derived pathogenic-only
    # file so it is rebuilt from the new release.
    mv -f "$tmp" "$CLINVAR_VCF"
    rm -f "$tmp.md5" "$CLINVAR_PATH_VCF" "${CLINVAR_PATH_VCF}.csi" "${CLINVAR_PATH_VCF}.tbi"
    DOWNLOAD_TMP=""
    success "ClinVar downloaded successfully!"
    return 0
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
hx_require_arg "$VCF" "$SCRIPT_NAME <vcf_file>"

# Check dependencies
check_dependencies

# Validate input VCF
validate_vcf "$VCF"

# Build (or reuse) a bgzip+indexed, bare-contig, GRCh38-checked working copy
# of the input (never indexes the user's file in place).
hx_prepare_vcf "$VCF"

# Create output and reference directories
mkdir -p "$OUTDIR" "$REF_DIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║            CLINVAR PATHOGENIC VARIANT SCANNER                             ║"
echo "║            Searching for pathogenic variants in your genome               ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "VCF File: $VCF"
echo "Results: $OUTDIR/"
echo "ClinVar cache: $REF_DIR/"
echo ""

#===============================================================================
# STEP 1: Download ClinVar (if not present or older than CLINVAR_MAX_AGE_DAYS)
# STEP 2: Extract pathogenic variants from ClinVar
# Both run under an exclusive flock on the cache directory so concurrent runs
# (e.g. two frontend sessions) never download or filter twice.
#===============================================================================

CLINVAR_VCF="$REF_DIR/clinvar_GRCh38.vcf.gz"
CLINVAR_PATH_VCF="$REF_DIR/clinvar_pathogenic_only.vcf.gz"

exec 9>"$REF_DIR/.clinvar.lock"
if command -v flock &> /dev/null; then
    flock -w 1800 9 || error_exit "another ClinVar download is in progress (waited 30 min on $REF_DIR/.clinvar.lock)"
else
    warn "flock not available (macOS?): concurrent ClinVar downloads are not serialised"
fi

if clinvar_cache_stale "$CLINVAR_VCF"; then
    if [[ -s "$CLINVAR_VCF" ]]; then
        echo "ℹ️  Cached ClinVar is older than $CLINVAR_MAX_AGE_DAYS days - refreshing..."
    fi
    if ! download_clinvar; then
        if [[ -s "$CLINVAR_VCF" ]]; then
            warn "Using the existing (stale) ClinVar copy in $REF_DIR"
        else
            error_exit "Failed to download ClinVar VCF from $CLINVAR_URL (see wget message above)"
        fi
    fi
else
    echo "ℹ️  Using existing ClinVar database"
fi

CLINVAR_DATE=$(clinvar_release_date "$CLINVAR_VCF")
echo "   ClinVar release (##fileDate): $CLINVAR_DATE"
echo ""

echo "🔬 Extracting Pathogenic/Likely_pathogenic variants from ClinVar..."

# Rebuild the derived file when missing, unindexed, or older than the source
if [[ ! -s "$CLINVAR_PATH_VCF" ]] || [[ "$CLINVAR_PATH_VCF" -ot "$CLINVAR_VCF" ]] \
    || { [[ ! -s "${CLINVAR_PATH_VCF}.csi" ]] && [[ ! -s "${CLINVAR_PATH_VCF}.tbi" ]]; }; then
    if ! bcftools view -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' \
        "$CLINVAR_VCF" -Oz -o "${CLINVAR_PATH_VCF}.tmp"; then
        rm -f "${CLINVAR_PATH_VCF}.tmp"
        error_exit "Failed to filter ClinVar for pathogenic variants"
    fi

    if ! bcftools index -f "${CLINVAR_PATH_VCF}.tmp"; then
        rm -f "${CLINVAR_PATH_VCF}.tmp" "${CLINVAR_PATH_VCF}.tmp.csi"
        error_exit "Failed to index filtered ClinVar file"
    fi

    mv -f "${CLINVAR_PATH_VCF}.tmp.csi" "${CLINVAR_PATH_VCF}.csi"
    mv -f "${CLINVAR_PATH_VCF}.tmp" "$CLINVAR_PATH_VCF"

    PATHOGENIC_COUNT=$(bcftools view -H "$CLINVAR_PATH_VCF" | wc -l)
    echo "   ClinVar contains $PATHOGENIC_COUNT Pathogenic/Likely_pathogenic variants"
else
    echo "   Using existing pathogenic variants file"
fi
echo ""

# Release the cache lock: the comparison below only reads the cache.
if command -v flock &> /dev/null; then
    flock -u 9
fi
exec 9>&-

#===============================================================================
# STEP 3: Find common variants (bcftools isec)
#===============================================================================

echo "🔍 Comparing your genome with ClinVar Pathogenic..."
echo ""

mkdir -p "$OUTDIR/isec_results"

# Number of records scanned (from the index; reported in the coverage line)
TOTAL_RECORDS=$(bcftools index -n "$HX_VCF" 2>/dev/null || echo 0)
[[ $TOTAL_RECORDS =~ ^[0-9]+$ ]] || TOTAL_RECORDS=0

# Find intersection - variants that are BOTH in your VCF AND in ClinVar
# Pathogenic (same position, REF and ALT). $HX_VCF is bgzip+indexed with bare
# contig names, matching the NCBI ClinVar VCF.
if ! bcftools isec -p "$OUTDIR/isec_results" -n=2 -w1 \
    "$HX_VCF" "$CLINVAR_PATH_VCF"; then
    error_exit "Failed to compare VCF with ClinVar"
fi

#===============================================================================
# STEP 4: Analyze results
#===============================================================================

ISEC_FILE="$OUTDIR/isec_results/0000.vcf"
RESULTS_FILE="$OUTDIR/isec_results/pathogenic_carried.vcf"
REPORT="$OUTDIR/PATHOGENIC_VARIANTS_FOUND.txt"

# isec matches on position and alleles only, so an explicit homozygous-
# reference record (0/0) at a pathogenic site would be listed as "found".
# Keep only records in which the sample actually carries an ALT allele.
if [[ -f "$ISEC_FILE" ]]; then
    if ! bcftools view -i 'GT="alt"' "$ISEC_FILE" -o "$RESULTS_FILE"; then
        error_exit "Failed to filter the ClinVar intersection by genotype"
    fi
fi

if [[ -f "$RESULTS_FILE" ]]; then
    # grep -c prints 0 AND exits 1 on no match; `|| echo 0` produced "0\n0"
    FOUND_COUNT=$(grep -cv "^#" "$RESULTS_FILE" || true)
    [[ $FOUND_COUNT =~ ^[0-9]+$ ]] || error_exit "Could not count results in $RESULTS_FILE"

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
                ClinVar release (##fileDate): $CLINVAR_DATE
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
                ClinVar release (##fileDate): $CLINVAR_DATE
================================================================================

⚠️ FOUND $FOUND_COUNT VARIANTS IN CLINVAR PATHOGENIC

DETAILS BELOW:
================================================================================

EOF

        echo "   Details of found variants:"
        echo "   ─────────────────────────────────────────"

        # Process up to 50 variants (top-level code: no `local`; and
        # `((count++))` from 0 returns 1, which aborts under set -e)
        count=0
        while IFS=$'\t' read -r CHROM POS _ REF ALT _ _ _ _ SAMPLE _; do
            [[ "$CHROM" == "#"* ]] && continue
            count=$((count + 1))
            [[ $count -gt 50 ]] && break

            GT="${SAMPLE%%:*}"

            # Get info from ClinVar (first record at the position; no `head`
            # in the pipeline, which would SIGPIPE under pipefail)
            CLINVAR_INFO=$(bcftools query -r "${CHROM}:${POS}-${POS}" \
                -i "POS=${POS} && REF=\"${REF}\" && ALT=\"${ALT}\"" \
                -f '%ID\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\n' \
                "$CLINVAR_PATH_VCF" 2>/dev/null || true)
            CLINVAR_INFO=${CLINVAR_INFO%%$'\n'*}

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
# This scan intersects whole files rather than querying fixed markers: the
# input record set counts as one queried "marker", found when the VCF carries
# at least one record (a header-only file yields 0/1 = "no usable data").
# shellcheck disable=SC2034  # read by hx_report_coverage (scripts/lib/vcf_common.sh)
HX_QUERIED=1
# shellcheck disable=SC2034  # read by hx_report_coverage
HX_FOUND=$([[ ${TOTAL_RECORDS:-0} -gt 0 ]] && echo 1 || echo 0)
hx_report_coverage
