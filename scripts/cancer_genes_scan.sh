#!/bin/bash
#===============================================================================
# COMPREHENSIVE CANCER PREDISPOSITION GENE ANALYSIS
# BRCA1/2, TP53, Lynch syndrome, and many others
#===============================================================================

set -euo pipefail

# Shared VCF helpers (hx_prepare_vcf, hx_gt_indel, hx_query, hx_count, ...)
# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"

#===============================================================================
# CONFIGURATION - Genomic coordinates (GRCh38, bare contig names; the shared
# library renames chr-prefixed input before any query)
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
    ["CDH1"]="16:68737292-68835537"  # TODO(build): unverified GRCh38 (RefSeq gene span; consistent with the ~-33.9 kb chr16 GRCh37->GRCh38 shift of rs17822931/rs9939609 in lib/GRCh38_COORDINATES.md; old value 68771117-68869444 was the GRCh37 span)
    ["CHEK2"]="22:28687743-28742422"
    ["PALB2"]="16:23603160-23641310"
    ["ATM"]="11:108222484-108369102"
)

# Listed BRCA1/BRCA2 mutation positions. Entries whose name contains
# del/ins/dup are matched against indel records only, the others against SNP
# records only (see check_snp), so a common SNP at an indel position is never
# reported as the named mutation.
# TODO(build): unverified GRCh38 - none of the six positions below could be
# checked against ClinVar/dbSNP offline. They are derived from the standard
# exon layout with the GRCh37->GRCh38 offsets implied by the gene regions
# above (BRCA1 +1,847,983; BRCA2 -574,137): 185delAG (c.68_69delAG, exon 2)
# ~17:43124027, 5382insC (c.5266dupC, exon 20) ~17:43057062, M1775R (c.5324T>G,
# exon 21) ~17:43051071, 6174delT (c.5946delT, exon 11) ~13:32340300, 4075delGT
# (c.3847_3848delGT, exon 11) ~13:32338202, 886delGT (c.658_659delGT, exon 8)
# ~13:32329469. Because a position-only match cannot confirm the named
# mutation, hits at these positions are printed as informational notes and
# are NOT counted as pathogenic findings (use clinvar_scan.sh for a verified
# call).
declare -A BRCA1_MUTATIONS=(
    ["17:43124027"]="185delAG|Ashkenazi"
    ["17:43057062"]="5382insC|Slavic"
    ["17:43051071"]="M1775R|BRCT"
)

declare -A BRCA2_MUTATIONS=(
    ["13:32340300"]="6174delT|Ashkenazi"
    ["13:32338202"]="4075delGT|Pathogenic"
    ["13:32329469"]="886delGT|Icelandic"
)

# TP53 hotspots: "NAME|TYPE|REF|ALT|RSID"; matched with hx_gt_allele so only a
# record with these alleles (and rsID, when the VCF carries one) scores. TP53 is
# on the minus strand, so the genomic alleles are the complement of the cDNA
# change. Positions derived from the GRCh37 dbSNP positions with the +96682
# offset implied by the gene region above.
# TODO(build): unverified GRCh38 - verify against ClinVar/dbSNP before relying
# on them. The former table put "R282W" at 17:7676154, which is rs1042522
# (Pro72Arg), a benign polymorphism most people carry; it is removed.
declare -A TP53_MUTATIONS=(
    ["17:7675088"]="R175H|Hotspot|C|T|rs28934578"
    ["17:7674220"]="R248Q|Hotspot|C|T|rs11540652"
    ["17:7674221"]="R248W|Hotspot|G|A|rs121912651"
    ["17:7673802"]="R273H|Hotspot|C|T|rs28934576"
    ["17:7673776"]="R282W|Pathogenic|G|A|rs28934574"
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
    if ! command -v bcftools &> /dev/null; then
        error_exit "Missing required dependencies: bcftools\nPlease run: ./install.sh"
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

# Genotype of the first called record at POS (any type: most entries in the
# mutation tables are small indels such as 185delAG or 1100delC, so the
# SNP-only hx_gt is not used here). Phase-normalised; empty when there is no
# record or the call is missing.
check_snp() { # POS [NAME]
    local rec
    case "${2:-}" in
        *del*|*ins*|*dup*)
            # Named indel: only a record with a length-changing ALT counts.
            rec=$(hx_gt_indel "$1" 1)
            printf '%s\n' "${rec##* }"
            ;;
        *)
            # Point mutation: SNP records only (hx_gt drops indels/MNPs).
            hx_gt "$1"
            ;;
    esac
}

# Number of records in a gene region. The region is also tallied as a
# coverage marker (present when it holds any called record, hom-ref ALT "."
# included): the 17 point markers below are pathogenic mutations that a
# healthy genome never carries, so without the region probes every clean run
# would report 0 found and the frontend would treat it as "no usable data".
check_region() {
    local region="$1"
    hx_region_present "$region"
    bcftools view -r "$region" -H "$HX_VCF" | wc -l
}

get_variants() {
    hx_query "$1" '%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%GT]\n'
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

# Build (or reuse) a bgzip+indexed, bare-contig, GRCh38-checked working copy;
# every query below uses $HX_VCF through the shared helpers.
hx_prepare_vcf "$VCF"

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

# No verified BRCA1/2 point probes exist (see the TODO(build) note at the
# tables), so these stay 0; hits at the listed positions are only noted.
BRCA1_FOUND=0
BRCA2_FOUND=0
BRCA_NOTED=0

# Check the listed BRCA1 mutation positions (informational only)
for pos in "${!BRCA1_MUTATIONS[@]}"; do
    IFS='|' read -r name pop <<< "${BRCA1_MUTATIONS[$pos]}"

    GT=$(check_snp "$pos" "$name")
    if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  ℹ️  BRCA1: variant at the listed $name position $pos ($GT, $pop) - unverified GRCh38 position, not confirmed as $name"
        ((BRCA_NOTED++)) || true
    fi
done

# Check the listed BRCA2 mutation positions (informational only)
for pos in "${!BRCA2_MUTATIONS[@]}"; do
    IFS='|' read -r name pop <<< "${BRCA2_MUTATIONS[$pos]}"

    GT=$(check_snp "$pos" "$name")
    if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  ℹ️  BRCA2: variant at the listed $name position $pos ($GT, $pop) - unverified GRCh38 position, not confirmed as $name"
        ((BRCA_NOTED++)) || true
    fi
done

if [[ $BRCA_NOTED -eq 0 ]]; then
    success "  ✅ No variants at the listed BRCA1/2 mutation positions (positions unverified; not a clinical BRCA call)"
else
    echo "   These positions are derived, not verified: confirm with clinvar_scan.sh or a clinical test."
fi

# Save all variants
get_variants "${GENE_REGIONS[BRCA1]}" > "$OUTDIR/BRCA1_all.txt"
get_variants "${GENE_REGIONS[BRCA2]}" > "$OUTDIR/BRCA2_all.txt"

# Records whose ID column (5th tab-separated field) is "." - hx_count always
# prints exactly one number (the old `grep -c ... || echo 0` printed "0\n0").
NO_RSID_PATTERN=$'^[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t\\.\t'
BRCA1_NOVEL=$(hx_count "$OUTDIR/BRCA1_all.txt" "$NO_RSID_PATTERN")
BRCA2_NOVEL=$(hx_count "$OUTDIR/BRCA2_all.txt" "$NO_RSID_PATTERN")

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
    IFS='|' read -r name type ref alt rsid <<< "${TP53_MUTATIONS[$pos]}"

    RES=$(hx_gt_allele "$pos" "$ref" "$alt" "$rsid")
    DOS=${RES%%$'\t'*}
    GT=${RES#*$'\t'}
    if [[ -n "$DOS" ]] && [[ "$DOS" != "0" ]]; then
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

    GT=$(check_snp "$pos" "$name")
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
GT=$(check_snp "22:28695868" 1100delC)  # TODO(build): unverified GRCh38
if [[ -n "$GT" ]] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
    warn "  ⚠️  CHEK2 1100delC: $GT (2x breast cancer risk)"
fi

# MUTYH - recessive
GT1=$(check_snp "1:45332879")  # TODO(build): unverified GRCh38 (position-only SNP probe)
GT2=$(check_snp "1:45340175")  # TODO(build): unverified GRCh38 (position-only SNP probe)
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
# Marker coverage (found/queried); the frontend treats 0 found as "no data".
hx_report_coverage
