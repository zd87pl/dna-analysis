#!/bin/bash
#===============================================================================
# CARRIER SCREENING ANALYSIS
# Screens for recessive disease carrier status
#
# DISCLAIMER: FOR EDUCATIONAL PURPOSES ONLY - NOT FOR CLINICAL USE
# Proper carrier screening requires clinical-grade testing
#===============================================================================

set -euo pipefail

# Shared VCF helpers (hx_prepare_vcf, hx_gt_allele, hx_gt_indel, ...)
# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"

# Validate input
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <vcf_file>"
    echo "Screens for common recessive disease carrier variants"
    exit 1
fi

VCF_FILE="$1"
hx_require_arg "$VCF_FILE" "$0 <vcf_file>"

if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: VCF file not found: $VCF_FILE" >&2
    exit 1
fi

# Check for bcftools
if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools is required but not installed" >&2
    exit 1
fi

# Build (or reuse) a bgzip+indexed, bare-contig, GRCh38-checked working copy;
# a bcftools failure now aborts with its message instead of reading as "no data".
hx_prepare_vcf "$VCF_FILE"

echo "==============================================================================="
echo "                       CARRIER SCREENING ANALYSIS"
echo "                    Recessive Disease Carrier Status"
echo "==============================================================================="
echo ""
echo "⚠️  DISCLAIMER: This is for EDUCATIONAL purposes only."
echo "    Clinical carrier screening requires validated laboratory testing."
echo "    Being a carrier typically does NOT affect your own health."
echo ""
echo "==============================================================================="

# Define carrier screening variants
# Format: rsID|Gene|Chr|Pos|Ref|Alt|Disease|Inheritance|Population|Severity
# Build: GRCh38, bare contig names, REF/ALT as on the plus strand of a VCF
# (canonical table: scripts/lib/GRCh38_COORDINATES.md). Rows with Alt "-"
# are deletions of the Ref bases starting at Pos; they are matched by the
# deleted length (see indel_dosage) because a VCF anchors them one base
# upstream. Rows without a canonical entry: TODO(build): unverified GRCh38.
declare -a CARRIER_VARIANTS=(
    # Cystic Fibrosis (CFTR) - Most common in Europeans
    "rs113993960|CFTR|7|117559590|CTT|-|Cystic Fibrosis|AR|European (1:25)|Severe"
    "rs75039782|CFTR|7|117587799|G|A|Cystic Fibrosis (G551D)|AR|European|Severe"
    "rs78655421|CFTR|7|117592219|G|A|Cystic Fibrosis (R553X)|AR|European|Severe"
    "rs74503330|CFTR|7|117540230|G|A|Cystic Fibrosis (W1282X)|AR|Ashkenazi (1:24)|Severe"

    # Sickle Cell Disease / Beta-Thalassemia (HBB)
    "rs334|HBB|11|5227002|T|A|Sickle Cell Disease|AR|African (1:12)|Severe"
    # HbC c.19G>A is the base next to HbS c.20A>T (rs334, 11:5227002) on the
    # minus-strand HBB; IVS-I-1 c.92+1G>A is 73 bp further 5'. Both derived
    # from the rs334 anchor.  TODO(build): unverified GRCh38 (dbSNP not consulted)
    "rs33930165|HBB|11|5227003|C|T|Hemoglobin C (HbC)|AR|West African|Variable"
    "rs33971440|HBB|11|5226929|C|T|Beta-Thalassemia (IVS-I-1)|AR|Mediterranean|Severe"

    # Alpha-Thalassemia (HBA1/HBA2)
    "rs41474145|HBA2|16|173245|G|A|Alpha-Thalassemia|AR|Southeast Asian,African|Variable"

    # Tay-Sachs Disease (HEXA) - Common in Ashkenazi Jews
    # rs76173977 is c.1274_1277dupTATC ("1278insTATC"), a 4-bp INSERTION: REF "-"
    # marks an insertion row, matched by a 4-bp length gain at POS-1/POS whatever
    # the anchor/strand representation (HEXA is on the minus strand).
    # TODO(build): unverified GRCh38 position
    "rs76173977|HEXA|15|72346580|-|TATC|Tay-Sachs Disease|AR|Ashkenazi (1:30)|Severe"
    "rs121907945|HEXA|15|72340467|G|A|Tay-Sachs Disease|AR|Ashkenazi|Severe"

    # Canavan Disease (ASPA) - Ashkenazi Jewish
    "rs28940279|ASPA|17|3483644|G|A|Canavan Disease|AR|Ashkenazi (1:40)|Severe"

    # Familial Dysautonomia (IKBKAP) - Ashkenazi Jewish
    "rs104942556|ELP1|9|108869675|T|C|Familial Dysautonomia|AR|Ashkenazi (1:30)|Severe"

    # Gaucher Disease (GBA) - Ashkenazi Jewish
    "rs76763715|GBA|1|155235843|T|C|Gaucher Disease (N370S)|AR|Ashkenazi (1:15)|Variable"
    "rs421016|GBA|1|155235252|A|G|Gaucher Disease (L444P)|AR|Ashkenazi|Moderate"

    # Niemann-Pick Disease (SMPD1)
    "rs120074117|SMPD1|11|6392706|G|A|Niemann-Pick Type A/B|AR|Ashkenazi|Severe"

    # Phenylketonuria (PAH)
    "rs5030858|PAH|12|102840632|G|A|Phenylketonuria (PKU)|AR|European (1:50)|Treatable"
    "rs5030849|PAH|12|102838369|C|T|Phenylketonuria (PKU)|AR|European|Treatable"

    # Spinal Muscular Atrophy (SMN1) - Note: usually detected by copy number
    "rs121909192|SMN1|5|70951946|C|T|Spinal Muscular Atrophy|AR|Pan-ethnic (1:50)|Severe"

    # Hereditary Hemochromatosis (HFE)
    "rs1800562|HFE|6|26092913|G|A|Hemochromatosis (C282Y)|AR|Northern European|Treatable"
    "rs1799945|HFE|6|26090951|C|G|Hemochromatosis (H63D)|AR|European|Mild"

    # Wilson Disease (ATP7B)
    "rs76151636|ATP7B|13|51937571|G|A|Wilson Disease|AR|Pan-ethnic (1:90)|Treatable"

    # Fanconi Anemia (FANCC) - Ashkenazi Jewish
    "rs104886456|FANCC|9|97862203|G|T|Fanconi Anemia|AR|Ashkenazi|Severe"

    # Bloom Syndrome (BLM) - Ashkenazi Jewish
    # blmAsh is a delins (c.2207_2212delATCTGAinsTAGATTC), matched by exact
    # REF/ALT; a caller that pads it with an anchor base (POS-1, GATCTGA>GTAGATTC)
    # is not matched.  TODO(build): unverified GRCh38 position
    "rs113993962|BLM|15|90631934|ATCTGA|TAGATTC|Bloom Syndrome|AR|Ashkenazi (1:100)|Severe"

    # Maple Syrup Urine Disease
    "rs121964998|BCKDHA|19|41557099|G|A|Maple Syrup Urine Disease|AR|Mennonite|Severe"

    # Medium-chain acyl-CoA dehydrogenase deficiency
    "rs77931234|ACADM|1|75724347|A|G|MCAD Deficiency|AR|Northern European|Treatable"

    # Glycogen Storage Disease Type 1a
    "rs80356483|G6PC|17|42910097|G|A|GSD Type 1a|AR|Ashkenazi|Treatable"

    # Connexin 26 Deafness (GJB2)
    "rs80338939|GJB2|13|20189473|C|-|Connexin 26 Deafness|AR|European (1:30)|Non-lethal"
    "rs72474224|GJB2|13|20189479|G|A|Connexin 26 Deafness|AR|European|Non-lethal"

    # Usher Syndrome (multiple genes)
    "rs111033256|MYO7A|11|76870294|G|A|Usher Syndrome Type 1B|AR|Pan-ethnic|Moderate"

    # Fragile X (FMR1) - X-linked, simplified
    # Note: Full Fragile X requires repeat expansion testing

    # Duchenne/Becker Muscular Dystrophy (DMD) - X-linked
    # Note: Usually requires deletion/duplication testing
)

echo ""
echo "📊 SCREENING FOR CARRIER VARIANTS..."
echo ""

# Counters
TOTAL_CHECKED=0
CARRIERS_FOUND=0
HOMOZYGOUS_FOUND=0

# Results storage by disease category
declare -a CARRIER_RESULTS=()
declare -a HOMOZYGOUS_RESULTS=()

# Dosage (0/1/2) of a deletion of SEQ starting at CHR:POS (table rows with
# Alt "-") or of an insertion of SEQ (rows with Ref "-"). A VCF anchors an
# indel one base upstream (POS-1, REF = anchor + deleted bases, ALT = anchor),
# so accept a called record at POS-1 or POS whose REF/ALT length difference
# equals the indel length. Prints nothing when no such record is called.
# MODE "del": ALT is REF shortened by len(SEQ); MODE "ins": ALT is REF
# lengthened by len(SEQ) (representation-independent: only the length change
# at the anchor is checked, so strand/anchor conventions do not matter).
indel_dosage() { # CHR POS SEQ MODE
    local chr=$1 pos=$2 seq=$3 mode=${4:-del} region rec rpos ref alt gt a i idx n want
    local -a alts gts
    region="${chr}:$((pos - 1))-${pos}"
    want=${#seq}
    [[ "$mode" == ins ]] && want=$(( -want ))
    # Tally the marker (queried/found) and bail out when no indel of at least
    # the expected length is called at the anchor.
    rec=$(hx_gt_indel "$region" "${#seq}")
    [[ -n "$rec" ]] || return 0
    # Scan every called record at the anchor for one ALT whose length change
    # matches exactly (a neighbouring, longer indel may be returned first).
    while IFS=$'\t' read -r rpos ref alt gt _; do
        [[ -n "$rpos" ]] || continue
        [[ "$rpos" -eq $((pos - 1)) || "$rpos" -eq "$pos" ]] || continue
        [[ "$gt" == *.* ]] && continue
        IFS=',' read -r -a alts <<< "$alt"
        idx=0; i=0
        for a in "${alts[@]}"; do
            i=$((i + 1))
            if [[ $(( ${#ref} - ${#a} )) -eq $want ]]; then
                idx=$i
                break
            fi
        done
        [[ $idx -gt 0 ]] || continue
        gt=${gt//|//}
        IFS='/' read -r -a gts <<< "$gt"
        n=0
        for a in "${gts[@]}"; do
            [[ "$a" == "$idx" ]] && n=$((n + 1))
        done
        printf '%s\n' "$n"
        return 0
    done < <(hx_query "$region" '%POS\t%REF\t%ALT\t[%GT\t]\n' -e 'GT="mis"')
    return 0
}

# Check each variant. Counters use X=$((X+1)): `((X++))` returns 1 when X is 0
# and aborted the script under set -e before the first lookup.
for variant_data in "${CARRIER_VARIANTS[@]}"; do
    IFS='|' read -r RSID GENE CHR POS REF ALT DISEASE _ POPULATION SEVERITY <<< "$variant_data"

    TOTAL_CHECKED=$((TOTAL_CHECKED + 1))

    if [[ "$ALT" == "-" ]]; then
        DOSAGE=$(indel_dosage "$CHR" "$POS" "$REF" del)
    elif [[ "$REF" == "-" ]]; then
        DOSAGE=$(indel_dosage "$CHR" "$POS" "$ALT" ins)
    else
        # The record must match POS, REF, ALT (and rsID when the VCF has one);
        # a different variant at the same coordinate is skipped with a
        # warning on stderr instead of being reported as this allele.
        RESULT=$(hx_gt_allele "${CHR}:${POS}-${POS}" "$REF" "$ALT" "$RSID")
        DOSAGE=${RESULT%%$'\t'*}
    fi
    [[ -n "$DOSAGE" ]] || continue

    # Check genotype
    if [[ "$DOSAGE" == "1" ]]; then
        CARRIERS_FOUND=$((CARRIERS_FOUND + 1))
        CARRIER_RESULTS+=("$RSID|$GENE|$DISEASE|$POPULATION|$SEVERITY")
    elif [[ "$DOSAGE" == "2" ]]; then
        HOMOZYGOUS_FOUND=$((HOMOZYGOUS_FOUND + 1))
        HOMOZYGOUS_RESULTS+=("$RSID|$GENE|$DISEASE|$POPULATION|$SEVERITY")
    fi
done

echo "==============================================================================="
echo "                              RESULTS SUMMARY"
echo "==============================================================================="
echo ""
echo "  Variants Screened:     $TOTAL_CHECKED"
echo "  Carrier Status Found:  $CARRIERS_FOUND"
echo "  Homozygous Variants:   $HOMOZYGOUS_FOUND"
echo ""

# Display homozygous findings first (more significant)
if [[ ${#HOMOZYGOUS_RESULTS[@]} -gt 0 ]]; then
    echo "==============================================================================="
    echo "  🔴 HOMOZYGOUS VARIANTS DETECTED - REVIEW RECOMMENDED"
    echo "==============================================================================="
    echo ""
    echo "  These variants are present in both copies of the gene."
    echo "  This may indicate affected status for recessive conditions."
    echo ""

    for result in "${HOMOZYGOUS_RESULTS[@]}"; do
        IFS='|' read -r RSID GENE DISEASE POPULATION SEVERITY <<< "$result"
        echo "  🔴 $DISEASE"
        echo "     Gene: $GENE ($RSID)"
        echo "     Carrier Frequency: $POPULATION"
        echo "     Condition Severity: $SEVERITY"
        echo ""
    done

    echo "  ⚠️  IMPORTANT: Homozygous pathogenic variants may indicate"
    echo "     disease status. Please consult a genetic counselor."
    echo ""
fi

# Display carrier findings
if [[ ${#CARRIER_RESULTS[@]} -gt 0 ]]; then
    echo "==============================================================================="
    echo "  🟡 CARRIER STATUS DETECTED"
    echo "==============================================================================="
    echo ""
    echo "  You carry one copy of variants associated with recessive conditions."
    echo "  Carriers are typically healthy but can pass the variant to children."
    echo ""

    for result in "${CARRIER_RESULTS[@]}"; do
        IFS='|' read -r RSID GENE DISEASE POPULATION SEVERITY <<< "$result"
        echo "  🟡 $DISEASE"
        echo "     Gene: $GENE ($RSID)"
        echo "     Carrier Frequency: $POPULATION"
        echo "     Condition Severity: $SEVERITY"
        echo ""
    done
fi

if [[ ${#CARRIER_RESULTS[@]} -eq 0 && ${#HOMOZYGOUS_RESULTS[@]} -eq 0 ]]; then
    echo "==============================================================================="
    echo "  ✅ NO CARRIER VARIANTS DETECTED"
    echo "==============================================================================="
    echo ""
    echo "  No carrier status detected for the variants analyzed."
    echo "  Note: This screen covers common variants only."
    echo ""
fi

# Reproductive risk information
if [[ ${#CARRIER_RESULTS[@]} -gt 0 ]]; then
    echo "==============================================================================="
    echo "                      REPRODUCTIVE CONSIDERATIONS"
    echo "==============================================================================="
    echo ""
    echo "  📋 If you are a carrier and your partner is also a carrier for the"
    echo "     SAME condition, each pregnancy has a:"
    echo ""
    echo "     • 25% chance of affected child"
    echo "     • 50% chance of carrier child"
    echo "     • 25% chance of non-carrier child"
    echo ""
    echo "  💡 RECOMMENDATION: If family planning is relevant, consider:"
    echo "     • Partner carrier testing for conditions you carry"
    echo "     • Genetic counseling consultation"
    echo "     • Clinical-grade carrier screening"
    echo ""
fi

# Calculate carrier score (lower is better)
CARRIER_SCORE=100
if [[ $HOMOZYGOUS_FOUND -gt 0 ]]; then
    CARRIER_SCORE=$((100 - (HOMOZYGOUS_FOUND * 30)))
    if [[ $CARRIER_SCORE -lt 0 ]]; then CARRIER_SCORE=0; fi
elif [[ $CARRIERS_FOUND -gt 3 ]]; then
    CARRIER_SCORE=60
elif [[ $CARRIERS_FOUND -gt 0 ]]; then
    CARRIER_SCORE=$((90 - (CARRIERS_FOUND * 10)))
fi

echo "==============================================================================="
echo "                           SCREENING SCORE"
echo "==============================================================================="
echo ""
echo "  CARRIER SCORE: $CARRIER_SCORE / 100"
echo ""
echo "  Score interpretation:"
echo "    90-100: Few or no carrier variants detected"
echo "    70-89:  Some carrier variants (common, usually not concerning)"
echo "    50-69:  Multiple carrier variants - consider partner screening"
echo "    <50:    Homozygous variants or many carriers - genetic counseling recommended"
echo ""

echo "==============================================================================="
echo "                             LIMITATIONS"
echo "==============================================================================="
echo ""
echo "  This screening has important limitations:"
echo ""
echo "  1. COVERAGE: Only screens for common variants in selected genes"
echo "  2. DETECTION: Cannot detect:"
echo "     - Copy number variations (deletions/duplications)"
echo "     - Repeat expansions (Fragile X, Huntington's)"
echo "     - Novel or rare mutations"
echo "  3. ETHNICITY: Carrier frequencies vary by population"
echo "  4. ACCURACY: VCF quality affects results"
echo ""
echo "  For comprehensive carrier screening, clinical testing is recommended."
echo ""
echo "==============================================================================="
echo "  ⚠️  NOT FOR CLINICAL USE - EDUCATIONAL PURPOSES ONLY"
echo "==============================================================================="
# Marker coverage (found/queried); the frontend treats 0 found as "no data".
hx_report_coverage
