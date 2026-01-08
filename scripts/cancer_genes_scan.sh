#!/bin/bash
#===============================================================================
# COMPREHENSIVE CANCER PREDISPOSITION GENE ANALYSIS
# BRCA1/2, TP53, Lynch syndrome, and many others
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION - Genomic coordinates (GRCh38)
#===============================================================================

SCRIPT_NAME="$(basename "$0")"
OUTDIR="cancer_genes_analysis"

# Gene regions (GRCh38 coordinates)
declare -A GENE_REGIONS=(
    ["BRCA1"]="17:43044295-43125483"
    ["BRCA2"]="13:32315474-32400266"
    ["TP53"]="17:7668402-7687538"
    ["MLH1"]="3:36993350-37050918"
    ["MSH2"]="2:47403067-47709917"
    ["MSH6"]="2:47695530-47810101"
    ["PMS2"]="7:5970925-6009106"
    ["PTEN"]="10:87863625-87971930"
    ["APC"]="5:112707498-112846239"
    ["CDH1"]="16:68771117-68869444"
    ["CHEK2"]="22:28687743-28742422"
    ["PALB2"]="16:23603160-23641310"
    ["ATM"]="11:108222484-108369102"
)

# Known pathogenic variants
declare -A BRCA1_MUTATIONS=(
    ["17:43045677"]="185delAG|Ashkenazi"
    ["17:43057051"]="5382insC|Slavic"
    ["17:43071077"]="M1775R|BRCT"
)

declare -A BRCA2_MUTATIONS=(
    ["13:32354860"]="6174delT|Ashkenazi"
    ["13:32336282"]="4075delGT|Pathogenic"
    ["13:32370557"]="886delGT|Icelandic"
)

declare -A TP53_MUTATIONS=(
    ["17:7674220"]="R175H|Hotspot"
    ["17:7675076"]="R248Q|Hotspot"
    ["17:7675994"]="R273H|Hotspot"
    ["17:7676154"]="R282W|Pathogenic"
)

declare -A LYNCH_MUTATIONS=(
    ["3:37042337"]="MLH1_c.350C>T"
    ["2:47641560"]="MSH2_c.1216C>T"
    ["2:47703303"]="MSH6_c.3261dup"
)

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
    for cmd in bcftools; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        error_exit "Missing required dependencies: ${missing[*]}\nPlease run: ./install.sh"
    fi
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

check_snp() {
    local pos="$1"
    bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1
}

check_region() {
    local region="$1"
    bcftools view -r "$region" -H "$VCF" 2>/dev/null | wc -l
}

get_variants() {
    local region="$1"
    bcftools query -r "$region" -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%GT]\n' "$VCF" 2>/dev/null
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
echo "║              COMPREHENSIVE CANCER PREDISPOSITION GENE ANALYSIS            ║"
echo "║                   BRCA1/2, TP53, Lynch, and others                        ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

#===============================================================================
# BRCA1 and BRCA2
#===============================================================================

echo "🎀 BRCA1 and BRCA2 - Breast and Ovarian Cancer"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

BRCA1_VARIANTS=$(check_region "${GENE_REGIONS[BRCA1]}")
BRCA2_VARIANTS=$(check_region "${GENE_REGIONS[BRCA2]}")

echo "   BRCA1: $BRCA1_VARIANTS variants"
echo "   BRCA2: $BRCA2_VARIANTS variants"
echo ""

BRCA1_FOUND=0
BRCA2_FOUND=0

# Check BRCA1 known mutations
for pos in "${!BRCA1_MUTATIONS[@]}"; do
    IFS='|' read -r name pop <<< "${BRCA1_MUTATIONS[$pos]}"

    GT=$(check_snp "$pos")
    if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 BRCA1 $name: $GT ($pop)"
        ((BRCA1_FOUND++)) || true
    fi
done

# Check BRCA2 known mutations
for pos in "${!BRCA2_MUTATIONS[@]}"; do
    IFS='|' read -r name pop <<< "${BRCA2_MUTATIONS[$pos]}"

    GT=$(check_snp "$pos")
    if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 BRCA2 $name: $GT ($pop)"
        ((BRCA2_FOUND++)) || true
    fi
done

if [[ $BRCA1_FOUND -eq 0 ]] && [[ $BRCA2_FOUND -eq 0 ]]; then
    success "  ✅ No known pathogenic BRCA1/2 mutations found"
fi

# Save all variants
get_variants "${GENE_REGIONS[BRCA1]}" > "$OUTDIR/BRCA1_all.txt"
get_variants "${GENE_REGIONS[BRCA2]}" > "$OUTDIR/BRCA2_all.txt"

BRCA1_NOVEL=$(grep -c '\.$' "$OUTDIR/BRCA1_all.txt" 2>/dev/null || echo "0")
BRCA2_NOVEL=$(grep -c '\.$' "$OUTDIR/BRCA2_all.txt" 2>/dev/null || echo "0")

echo ""
echo "   Variants without rsID (potentially rare):"
echo "   BRCA1: $BRCA1_NOVEL | BRCA2: $BRCA2_NOVEL"
echo ""

#===============================================================================
# TP53 (Li-Fraumeni)
#===============================================================================

echo "🧬 TP53 - Li-Fraumeni Syndrome"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

TP53_VARIANTS=$(check_region "${GENE_REGIONS[TP53]}")
echo "   TP53: $TP53_VARIANTS variants"
echo ""

TP53_FOUND=0

for pos in "${!TP53_MUTATIONS[@]}"; do
    IFS='|' read -r name type <<< "${TP53_MUTATIONS[$pos]}"

    GT=$(check_snp "$pos")
    if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 TP53 $name: $GT ($type)"
        ((TP53_FOUND++)) || true
    fi
done

if [[ $TP53_FOUND -eq 0 ]]; then
    success "  ✅ No pathogenic TP53 mutations found"
fi
echo ""

#===============================================================================
# Lynch Syndrome
#===============================================================================

echo "🔵 Lynch Syndrome - MMR genes"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

MLH1=$(check_region "${GENE_REGIONS[MLH1]}")
MSH2=$(check_region "${GENE_REGIONS[MSH2]}")
MSH6=$(check_region "${GENE_REGIONS[MSH6]}")
PMS2=$(check_region "${GENE_REGIONS[PMS2]}")

echo "   MLH1: $MLH1 | MSH2: $MSH2 | MSH6: $MSH6 | PMS2: $PMS2"
echo ""

LYNCH_FOUND=0

for pos in "${!LYNCH_MUTATIONS[@]}"; do
    name="${LYNCH_MUTATIONS[$pos]}"

    GT=$(check_snp "$pos")
    if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 $name: $GT"
        ((LYNCH_FOUND++)) || true
    fi
done

if [[ $LYNCH_FOUND -eq 0 ]]; then
    success "  ✅ No pathogenic Lynch syndrome mutations found"
fi
echo ""

#===============================================================================
# Other cancer genes
#===============================================================================

echo "🧬 OTHER CANCER PREDISPOSITION GENES"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "   GENE       | VARIANTS | CANCERS"
echo "   ─────────────────────────────────────────────────────────"

declare -A OTHER_GENES=(
    ["PTEN"]="Breast, thyroid"
    ["APC"]="Colorectal"
    ["CDH1"]="Gastric"
    ["CHEK2"]="Breast"
    ["PALB2"]="Breast, pancreatic"
    ["ATM"]="Breast"
)

for gene in PTEN APC CDH1 CHEK2 PALB2 ATM; do
    region="${GENE_REGIONS[$gene]}"
    cancers="${OTHER_GENES[$gene]}"
    count=$(check_region "$region")
    printf "   %-10s | %8d | %s\n" "$gene" "$count" "$cancers"
done

echo ""

# CHEK2 1100delC - special variant
GT=$(check_snp "22:28695868")
if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
    warn "  ⚠️  CHEK2 1100delC: $GT (2x breast cancer risk)"
fi

# MUTYH - recessive
GT1=$(check_snp "1:45332879")
GT2=$(check_snp "1:45340175")
MUTYH_COUNT=0
[[ "${GT1:-}" == *"1"* ]] && ((MUTYH_COUNT++)) || true
[[ "${GT2:-}" == *"1"* ]] && ((MUTYH_COUNT++)) || true

if [[ $MUTYH_COUNT -gt 0 ]]; then
    echo "  ℹ️  MUTYH: $MUTYH_COUNT variant(s) (recessive - need 2 for disease)"
fi

echo ""

#===============================================================================
# SUMMARY
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                              SUMMARY                                      ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

TOTAL=$((BRCA1_FOUND + BRCA2_FOUND + TP53_FOUND + LYNCH_FOUND))

echo "   BRCA1 pathogenic:  $BRCA1_FOUND"
echo "   BRCA2 pathogenic:  $BRCA2_FOUND"
echo "   TP53 pathogenic:   $TP53_FOUND"
echo "   Lynch pathogenic:  $LYNCH_FOUND"
echo "   ────────────────────"
echo "   TOTAL:             $TOTAL"
echo ""

if [[ $TOTAL -eq 0 ]]; then
    success "✅ NO known pathogenic mutations found in cancer genes!"
    echo ""
    echo "   Remember: ~90% of cancers are sporadic, not hereditary."
    echo "   Regular screening is important regardless of genetics."
else
    echo -e "${RED}⚠️  FOUND $TOTAL pathogenic variants!${NC}"
    echo "   → Clinical geneticist consultation recommended"
fi

echo ""
echo "📁 Files: $OUTDIR/BRCA1_all.txt, BRCA2_all.txt"
echo ""
