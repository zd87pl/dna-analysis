#!/bin/bash
#===============================================================================
# PHARMACOGENOMICS (PGx) ANALYSIS
# Analyzes genetic variants affecting drug metabolism
#
# DISCLAIMER: FOR EDUCATIONAL PURPOSES ONLY - NOT FOR CLINICAL USE
# Always consult healthcare providers for medication decisions
#===============================================================================

set -euo pipefail

# Validate input
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <vcf_file>"
    echo "Analyzes pharmacogenomic variants affecting drug metabolism"
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
echo "                    PHARMACOGENOMICS (PGx) ANALYSIS"
echo "                    Drug Metabolism Genetic Variants"
echo "==============================================================================="
echo ""
echo "⚠️  DISCLAIMER: This is for EDUCATIONAL purposes only."
echo "    Do NOT use this information to modify any medication regimen."
echo "    Always consult your physician or pharmacist."
echo ""
echo "==============================================================================="

# Define PGx variants with clinical annotations
# Format: rsID|Gene|Chr|Pos|Ref|Alt|Drug_Class|Effect|Clinical_Note
declare -a PGX_VARIANTS=(
    # CYP2D6 - Metabolizes ~25% of all drugs
    "rs3892097|CYP2D6|22|42126611|G|A|Antidepressants,Opioids,Beta-blockers|Poor Metabolizer|Reduced enzyme activity - may need dose adjustment"
    "rs1065852|CYP2D6|22|42130692|G|A|Codeine,Tramadol,Tamoxifen|Reduced Function|Intermediate metabolizer phenotype"
    "rs16947|CYP2D6|22|42127941|G|A|Multiple drugs|Normal/Increased|Common variant, normal function"

    # CYP2C19 - Clopidogrel, PPIs, antidepressants
    "rs4244285|CYP2C19|10|94781859|G|A|Clopidogrel,PPIs,SSRIs|Poor Metabolizer|*2 allele - loss of function"
    "rs4986893|CYP2C19|10|94780653|G|A|Clopidogrel,Voriconazole|Poor Metabolizer|*3 allele - loss of function"
    "rs12248560|CYP2C19|10|94761900|C|T|Clopidogrel,PPIs|Rapid Metabolizer|*17 allele - increased activity"

    # CYP2C9 - Warfarin, NSAIDs
    "rs1799853|CYP2C9|10|94942290|C|T|Warfarin,NSAIDs,Sulfonylureas|Reduced Function|*2 allele - 30% reduced activity"
    "rs1057910|CYP2C9|10|94981296|A|C|Warfarin,Phenytoin|Reduced Function|*3 allele - 80% reduced activity"

    # CYP3A4/CYP3A5 - Largest drug-metabolizing enzyme family
    "rs35599367|CYP3A4|7|99768693|C|T|Statins,Immunosuppressants|Reduced Function|*22 allele - reduced expression"
    "rs776746|CYP3A5|7|99672916|T|C|Tacrolimus,Cyclosporine|Non-expresser|*3 allele - splicing defect"

    # VKORC1 - Warfarin sensitivity
    "rs9923231|VKORC1|16|31096368|G|A|Warfarin|Increased Sensitivity|Lower warfarin dose needed"

    # SLCO1B1 - Statin transport
    "rs4149056|SLCO1B1|12|21176804|T|C|Statins|Reduced Transport|Increased myopathy risk with simvastatin"

    # TPMT - Thiopurines (azathioprine, 6-MP)
    "rs1800462|TPMT|6|18130918|C|G|Azathioprine,6-MP|Reduced Function|*2 allele - toxicity risk"
    "rs1800460|TPMT|6|18130687|C|T|Azathioprine,6-MP|Reduced Function|*3B allele - toxicity risk"
    "rs1142345|TPMT|6|18130348|T|C|Azathioprine,6-MP|Reduced Function|*3C allele - toxicity risk"

    # DPYD - Fluoropyrimidines (5-FU, capecitabine)
    "rs3918290|DPYD|1|97915614|C|T|5-FU,Capecitabine|Deficient|*2A - severe toxicity risk"
    "rs55886062|DPYD|1|98039419|A|C|5-FU,Capecitabine|Deficient|*13 - severe toxicity risk"

    # UGT1A1 - Irinotecan metabolism
    "rs8175347|UGT1A1|2|233757013|TA6|TA7|Irinotecan|Reduced Function|*28 allele - neutropenia risk"

    # CYP1A2 - Caffeine, theophylline, some antipsychotics
    "rs762551|CYP1A2|15|74749576|A|C|Caffeine,Clozapine|Fast Metabolizer|*1F - induced by smoking"

    # NAT2 - Isoniazid, hydralazine, sulfonamides
    "rs1801280|NAT2|8|18257854|T|C|Isoniazid,Hydralazine|Slow Acetylator|*5 allele"
    "rs1799930|NAT2|8|18258103|G|A|Isoniazid,Sulfonamides|Slow Acetylator|*6 allele"
    "rs1799931|NAT2|8|18258370|G|A|Isoniazid,Caffeine|Slow Acetylator|*7 allele"

    # OPRM1 - Opioid response
    "rs1799971|OPRM1|6|154039662|A|G|Opioids|Altered Response|A118G - may need higher doses"

    # COMT - Pain sensitivity, some psych meds
    "rs4680|COMT|22|19963748|G|A|Pain medications,ADHD drugs|Altered Metabolism|Val158Met - affects pain perception"

    # HLA-B*57:01 - Abacavir hypersensitivity
    "rs2395029|HCP5|6|31431780|T|G|Abacavir|Hypersensitivity Risk|HLA-B*57:01 tag - screen before prescribing"

    # HLA-B*15:02 - Carbamazepine (SJS/TEN risk in Asians)
    "rs144012689|HLA-B|6|31356811|C|T|Carbamazepine|SJS/TEN Risk|HLA-B*15:02 tag SNP"

    # Factor V Leiden - Contraceptive risk
    "rs6025|F5|1|169549811|G|A|Oral Contraceptives|Clotting Risk|Factor V Leiden - VTE risk with estrogen"
)

echo ""
echo "📊 ANALYZING PHARMACOGENOMIC VARIANTS..."
echo ""

# Counters
TOTAL_CHECKED=0
VARIANTS_FOUND=0
ACTIONABLE=0

# Results storage
declare -a RESULTS=()

# Check each variant
for variant_data in "${PGX_VARIANTS[@]}"; do
    IFS='|' read -r RSID GENE CHR POS REF ALT DRUGS EFFECT NOTE <<< "$variant_data"

    ((TOTAL_CHECKED++))

    # Query the VCF for this position
    RESULT=$(bcftools query -r "${CHR}:${POS}-${POS}" -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")

    if [[ -n "$RESULT" ]]; then
        GT=$(echo "$RESULT" | cut -f5)
        VCF_ALT=$(echo "$RESULT" | cut -f4)

        # Check if variant allele is present
        if [[ "$GT" == "0/1" || "$GT" == "1/0" || "$GT" == "0|1" || "$GT" == "1|0" ]]; then
            ((VARIANTS_FOUND++))
            STATUS="HETEROZYGOUS"
            RESULTS+=("HET|$RSID|$GENE|$DRUGS|$EFFECT|$NOTE")

            if [[ "$EFFECT" == *"Poor"* || "$EFFECT" == *"Deficient"* || "$EFFECT" == *"Risk"* ]]; then
                ((ACTIONABLE++))
            fi
        elif [[ "$GT" == "1/1" || "$GT" == "1|1" ]]; then
            ((VARIANTS_FOUND++))
            STATUS="HOMOZYGOUS"
            RESULTS+=("HOM|$RSID|$GENE|$DRUGS|$EFFECT|$NOTE")

            if [[ "$EFFECT" == *"Poor"* || "$EFFECT" == *"Deficient"* || "$EFFECT" == *"Risk"* ]]; then
                ((ACTIONABLE++))
            fi
        fi
    fi
done

echo "==============================================================================="
echo "                              RESULTS SUMMARY"
echo "==============================================================================="
echo ""
echo "  Variants Checked:    $TOTAL_CHECKED"
echo "  Variants Found:      $VARIANTS_FOUND"
echo "  Potentially Actionable: $ACTIONABLE"
echo ""

if [[ ${#RESULTS[@]} -gt 0 ]]; then
    echo "==============================================================================="
    echo "                         DETAILED FINDINGS"
    echo "==============================================================================="
    echo ""

    # Group by drug class
    echo "📋 VARIANT DETAILS:"
    echo ""

    for result in "${RESULTS[@]}"; do
        IFS='|' read -r ZYGOSITY RSID GENE DRUGS EFFECT NOTE <<< "$result"

        if [[ "$ZYGOSITY" == "HOM" ]]; then
            echo "🔴 $RSID ($GENE) - HOMOZYGOUS"
        else
            echo "🟡 $RSID ($GENE) - HETEROZYGOUS"
        fi
        echo "   Drug Classes: $DRUGS"
        echo "   Effect: $EFFECT"
        echo "   Note: $NOTE"
        echo ""
    done

    echo "==============================================================================="
    echo "                      DRUG-SPECIFIC IMPLICATIONS"
    echo "==============================================================================="
    echo ""

    # Check for specific high-impact combinations
    for result in "${RESULTS[@]}"; do
        IFS='|' read -r ZYGOSITY RSID GENE DRUGS EFFECT NOTE <<< "$result"

        if [[ "$GENE" == "CYP2C19" && "$EFFECT" == *"Poor"* ]]; then
            echo "⚠️  CLOPIDOGREL (Plavix): Reduced activation - may need alternative antiplatelet"
        fi

        if [[ "$GENE" == "CYP2D6" && "$EFFECT" == *"Poor"* ]]; then
            echo "⚠️  CODEINE/TRAMADOL: Reduced conversion to active metabolite"
            echo "⚠️  TAMOXIFEN: May have reduced efficacy"
        fi

        if [[ "$GENE" == "VKORC1" || "$GENE" == "CYP2C9" ]]; then
            echo "⚠️  WARFARIN: Dose adjustment likely needed - consult pharmacist"
        fi

        if [[ "$GENE" == "SLCO1B1" && "$ZYGOSITY" == "HOM" ]]; then
            echo "⚠️  SIMVASTATIN: Higher myopathy risk - consider alternative statin"
        fi

        if [[ "$GENE" == "DPYD" ]]; then
            echo "🔴 5-FU/CAPECITABINE: DPYD deficiency detected - CRITICAL toxicity risk"
        fi

        if [[ "$GENE" == "TPMT" ]]; then
            echo "⚠️  AZATHIOPRINE/6-MP: Reduced TPMT - dose reduction required"
        fi

        if [[ "$RSID" == "rs2395029" && "$ZYGOSITY" != "" ]]; then
            echo "🔴 ABACAVIR: HLA-B*57:01 positive - hypersensitivity risk - AVOID"
        fi

        if [[ "$RSID" == "rs6025" ]]; then
            echo "⚠️  ORAL CONTRACEPTIVES: Factor V Leiden - increased VTE risk"
        fi
    done
    echo ""
else
    echo "✅ No significant pharmacogenomic variants detected in analyzed positions."
    echo ""
fi

# Calculate score
SCORE=0
if [[ $VARIANTS_FOUND -eq 0 ]]; then
    SCORE=100
elif [[ $ACTIONABLE -eq 0 ]]; then
    SCORE=90
elif [[ $ACTIONABLE -le 2 ]]; then
    SCORE=70
elif [[ $ACTIONABLE -le 5 ]]; then
    SCORE=50
else
    SCORE=30
fi

echo "==============================================================================="
echo "                           METABOLIZER PROFILE"
echo "==============================================================================="
echo ""
echo "  PGx SCORE: $SCORE / 100"
echo ""
echo "  Score interpretation:"
echo "    90-100: Standard drug metabolism expected"
echo "    70-89:  Minor variations - some drugs may need attention"
echo "    50-69:  Moderate variations - discuss with pharmacist"
echo "    <50:    Multiple variations - pharmacogenomic consultation recommended"
echo ""

echo "==============================================================================="
echo "                              IMPORTANT NOTES"
echo "==============================================================================="
echo ""
echo "  1. This analysis covers common PGx variants only"
echo "  2. Many genes (especially CYP2D6) have complex variations not fully"
echo "     captured by single SNP analysis"
echo "  3. Drug response is also affected by:"
echo "     - Other medications (drug-drug interactions)"
echo "     - Age, weight, kidney/liver function"
echo "     - Diet and lifestyle factors"
echo "  4. Always consult a healthcare provider before making medication changes"
echo ""
echo "  📚 Resources:"
echo "     - PharmGKB: https://www.pharmgkb.org"
echo "     - CPIC Guidelines: https://cpicpgx.org"
echo ""
echo "==============================================================================="
echo "  ⚠️  NOT FOR CLINICAL USE - EDUCATIONAL PURPOSES ONLY"
echo "==============================================================================="
