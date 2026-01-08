#!/bin/bash
#===============================================================================
# CARRIER SCREENING ANALYSIS
# Screens for recessive disease carrier status
#
# DISCLAIMER: FOR EDUCATIONAL PURPOSES ONLY - NOT FOR CLINICAL USE
# Proper carrier screening requires clinical-grade testing
#===============================================================================

set -euo pipefail

# Validate input
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <vcf_file>"
    echo "Screens for common recessive disease carrier variants"
    exit 1
fi

VCF_FILE="$1"

if [[ ! -f "$VCF_FILE" ]]; then
    echo "Error: VCF file not found: $VCF_FILE"
    exit 1
fi

# Check for bcftools
if ! command -v bcftools &> /dev/null; then
    echo "Error: bcftools is required but not installed"
    exit 1
fi

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
declare -a CARRIER_VARIANTS=(
    # Cystic Fibrosis (CFTR) - Most common in Europeans
    "rs113993960|CFTR|7|117559590|CTT|-|Cystic Fibrosis|AR|European (1:25)|Severe"
    "rs75039782|CFTR|7|117587799|G|A|Cystic Fibrosis (G551D)|AR|European|Severe"
    "rs78655421|CFTR|7|117592219|G|A|Cystic Fibrosis (R553X)|AR|European|Severe"
    "rs74503330|CFTR|7|117540230|G|A|Cystic Fibrosis (W1282X)|AR|Ashkenazi (1:24)|Severe"

    # Sickle Cell Disease / Beta-Thalassemia (HBB)
    "rs334|HBB|11|5227002|T|A|Sickle Cell Disease|AR|African (1:12)|Severe"
    "rs33930165|HBB|11|5226773|C|T|Beta-Thalassemia|AR|Mediterranean,Asian|Severe"
    "rs33971440|HBB|11|5226862|C|A|Beta-Thalassemia|AR|Southeast Asian|Severe"

    # Alpha-Thalassemia (HBA1/HBA2)
    "rs41474145|HBA2|16|173245|G|A|Alpha-Thalassemia|AR|Southeast Asian,African|Variable"

    # Tay-Sachs Disease (HEXA) - Common in Ashkenazi Jews
    "rs76173977|HEXA|15|72346580|TATC|-|Tay-Sachs Disease|AR|Ashkenazi (1:30)|Severe"
    "rs121907945|HEXA|15|72340467|G|A|Tay-Sachs Disease|AR|Ashkenazi|Severe"

    # Canavan Disease (ASPA) - Ashkenazi Jewish
    "rs28940279|ASPA|17|3483644|G|A|Canavan Disease|AR|Ashkenazi (1:40)|Severe"

    # Familial Dysautonomia (IKBKAP) - Ashkenazi Jewish
    "rs104942556|ELP1|9|108869675|T|C|Familial Dysautonomia|AR|Ashkenazi (1:30)|Severe"

    # Gaucher Disease (GBA) - Ashkenazi Jewish
    "rs76763715|GBA|1|155204867|T|C|Gaucher Disease|AR|Ashkenazi (1:15)|Variable"
    "rs421016|GBA|1|155205634|G|A|Gaucher Disease (N370S)|AR|Ashkenazi|Moderate"

    # Niemann-Pick Disease (SMPD1)
    "rs120074117|SMPD1|11|6392706|G|A|Niemann-Pick Type A/B|AR|Ashkenazi|Severe"

    # Phenylketonuria (PAH)
    "rs5030858|PAH|12|102840632|G|A|Phenylketonuria (PKU)|AR|European (1:50)|Treatable"
    "rs5030849|PAH|12|102838369|C|T|Phenylketonuria (PKU)|AR|European|Treatable"

    # Spinal Muscular Atrophy (SMN1) - Note: usually detected by copy number
    "rs121909192|SMN1|5|70951946|C|T|Spinal Muscular Atrophy|AR|Pan-ethnic (1:50)|Severe"

    # Hereditary Hemochromatosis (HFE)
    "rs1800562|HFE|6|26091179|G|A|Hemochromatosis (C282Y)|AR|Northern European|Treatable"
    "rs1799945|HFE|6|26090951|C|G|Hemochromatosis (H63D)|AR|European|Mild"

    # Wilson Disease (ATP7B)
    "rs76151636|ATP7B|13|51937571|G|A|Wilson Disease|AR|Pan-ethnic (1:90)|Treatable"

    # Fanconi Anemia (FANCC) - Ashkenazi Jewish
    "rs104886456|FANCC|9|97862203|G|T|Fanconi Anemia|AR|Ashkenazi|Severe"

    # Bloom Syndrome (BLM) - Ashkenazi Jewish
    "rs113993962|BLM|15|90631934|ATCTGA|-|Bloom Syndrome|AR|Ashkenazi (1:100)|Severe"

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

# Check each variant
for variant_data in "${CARRIER_VARIANTS[@]}"; do
    IFS='|' read -r RSID GENE CHR POS REF ALT DISEASE INHERITANCE POPULATION SEVERITY <<< "$variant_data"

    ((TOTAL_CHECKED++))

    # Query the VCF for this position
    RESULT=$(bcftools query -r "${CHR}:${POS}-${POS}" -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")

    if [[ -n "$RESULT" ]]; then
        GT=$(echo "$RESULT" | cut -f5)

        # Check genotype
        if [[ "$GT" == "0/1" || "$GT" == "1/0" || "$GT" == "0|1" || "$GT" == "1|0" ]]; then
            ((CARRIERS_FOUND++))
            CARRIER_RESULTS+=("$RSID|$GENE|$DISEASE|$POPULATION|$SEVERITY")
        elif [[ "$GT" == "1/1" || "$GT" == "1|1" ]]; then
            ((HOMOZYGOUS_FOUND++))
            HOMOZYGOUS_RESULTS+=("$RSID|$GENE|$DISEASE|$POPULATION|$SEVERITY")
        fi
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
