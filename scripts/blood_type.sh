#!/bin/bash
#===============================================================================
# BLOOD TYPE PREDICTION
# Predicts ABO and Rh blood type from genetic variants
#
# DISCLAIMER: FOR EDUCATIONAL PURPOSES ONLY
# Always verify blood type through clinical laboratory testing
#===============================================================================

set -euo pipefail

# Validate input
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <vcf_file>"
    echo "Predicts ABO and Rh blood type from genetic variants"
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
echo "                       BLOOD TYPE PREDICTION"
echo "                      ABO and Rh Factor Analysis"
echo "==============================================================================="
echo ""
echo "⚠️  DISCLAIMER: This is a genetic PREDICTION only."
echo "    Always confirm blood type with laboratory testing before"
echo "    transfusion or medical procedures."
echo ""
echo "==============================================================================="

# ABO Blood Type Genetics:
# - ABO gene on chromosome 9
# - Three alleles: A, B, O
# - A and B are codominant, O is recessive
#
# Key variants:
# rs8176719 (261delG) - Deletion causes O allele
# rs8176746 (796C>A) - Distinguishes B from A
# rs8176747 (803G>C) - Distinguishes B from A
# rs7853989 - Additional A/B distinction

# Rh Factor Genetics:
# - RHD gene on chromosome 1
# - Rh+ (D antigen present) is dominant
# - Rh- (D antigen absent) usually due to RHD deletion
#
# Key variants:
# rs590787 - RHD presence/absence proxy

echo ""
echo "📊 ANALYZING BLOOD TYPE VARIANTS..."
echo ""

# Initialize allele tracking
A_ALLELE=0
B_ALLELE=0
O_ALLELE=0
RH_POSITIVE=0
RH_UNKNOWN=1

# Check ABO variants

# rs8176719 - 261delG (O allele marker)
# Reference G, deletion causes O allele
O_RESULT=$(bcftools query -r "9:133257521-133257521" -f '[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")
if [[ -n "$O_RESULT" ]]; then
    if [[ "$O_RESULT" == "1/1" || "$O_RESULT" == "1|1" ]]; then
        O_ALLELE=2  # Homozygous O
    elif [[ "$O_RESULT" == "0/1" || "$O_RESULT" == "1/0" || "$O_RESULT" == "0|1" || "$O_RESULT" == "1|0" ]]; then
        O_ALLELE=1  # Heterozygous, one O allele
    fi
fi

# rs8176746 (796C>A) - B allele marker
# A allele = C, B allele = A at this position
B_MARKER1=$(bcftools query -r "9:133256264-133256264" -f '%REF\t%ALT\t[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")
if [[ -n "$B_MARKER1" ]]; then
    GT=$(echo "$B_MARKER1" | cut -f3)
    if [[ "$GT" == "1/1" || "$GT" == "1|1" ]]; then
        B_ALLELE=2
    elif [[ "$GT" == "0/1" || "$GT" == "1/0" || "$GT" == "0|1" || "$GT" == "1|0" ]]; then
        B_ALLELE=1
    fi
fi

# rs8176747 (803G>C) - Another B allele marker
B_MARKER2=$(bcftools query -r "9:133256257-133256257" -f '%REF\t%ALT\t[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")

# rs7853989 - A1 vs A2 distinction and A/B
A_MARKER=$(bcftools query -r "9:133256042-133256042" -f '%REF\t%ALT\t[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")
if [[ -n "$A_MARKER" ]]; then
    GT=$(echo "$A_MARKER" | cut -f3)
    if [[ "$GT" == "0/0" || "$GT" == "0|0" ]]; then
        # Reference allele associated with A
        if [[ $B_ALLELE -eq 0 ]]; then
            A_ALLELE=$((2 - O_ALLELE))
        fi
    fi
fi

# Check Rh factor

# rs590787 - Proxy for RHD (not perfect but commonly used)
RH_RESULT=$(bcftools query -r "1:25362249-25362249" -f '[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")
if [[ -n "$RH_RESULT" ]]; then
    RH_UNKNOWN=0
    if [[ "$RH_RESULT" == "0/0" || "$RH_RESULT" == "0|0" ]]; then
        RH_POSITIVE=0  # Likely Rh-
    else
        RH_POSITIVE=1  # Likely Rh+
    fi
fi

# Additional RHD check - rs676785
RH_RESULT2=$(bcftools query -r "1:25409008-25409008" -f '[%GT]\n' "$VCF_FILE" 2>/dev/null || echo "")
if [[ -n "$RH_RESULT2" && $RH_UNKNOWN -eq 1 ]]; then
    RH_UNKNOWN=0
    if [[ "$RH_RESULT2" == "1/1" || "$RH_RESULT2" == "1|1" ]]; then
        RH_POSITIVE=1
    fi
fi

# Determine ABO type
ABO_TYPE="Unknown"
ABO_GENOTYPE="Unknown"

if [[ $O_ALLELE -eq 2 ]]; then
    ABO_TYPE="O"
    ABO_GENOTYPE="OO"
elif [[ $B_ALLELE -ge 1 && $A_ALLELE -eq 0 ]]; then
    if [[ $O_ALLELE -eq 1 ]]; then
        ABO_TYPE="B"
        ABO_GENOTYPE="BO"
    elif [[ $B_ALLELE -eq 2 ]]; then
        ABO_TYPE="B"
        ABO_GENOTYPE="BB"
    else
        ABO_TYPE="B"
        ABO_GENOTYPE="B?"
    fi
elif [[ $A_ALLELE -ge 1 && $B_ALLELE -eq 0 ]]; then
    if [[ $O_ALLELE -eq 1 ]]; then
        ABO_TYPE="A"
        ABO_GENOTYPE="AO"
    elif [[ $A_ALLELE -eq 2 ]]; then
        ABO_TYPE="A"
        ABO_GENOTYPE="AA"
    else
        ABO_TYPE="A"
        ABO_GENOTYPE="A?"
    fi
elif [[ $A_ALLELE -ge 1 && $B_ALLELE -ge 1 ]]; then
    ABO_TYPE="AB"
    ABO_GENOTYPE="AB"
elif [[ $O_ALLELE -eq 0 && $A_ALLELE -eq 0 && $B_ALLELE -eq 0 ]]; then
    # No O markers, check if we have A by default (reference)
    # Most common type is A or O
    ABO_TYPE="A (probable)"
    ABO_GENOTYPE="Likely A"
fi

# Determine Rh factor
if [[ $RH_UNKNOWN -eq 1 ]]; then
    RH_TYPE="Unknown"
elif [[ $RH_POSITIVE -eq 1 ]]; then
    RH_TYPE="Positive (+)"
else
    RH_TYPE="Negative (-)"
fi

# Combine for full blood type
if [[ "$ABO_TYPE" != "Unknown" && "$RH_TYPE" != "Unknown" ]]; then
    if [[ $RH_POSITIVE -eq 1 ]]; then
        FULL_TYPE="${ABO_TYPE}+"
    else
        FULL_TYPE="${ABO_TYPE}-"
    fi
else
    FULL_TYPE="${ABO_TYPE} / Rh ${RH_TYPE}"
fi

echo "==============================================================================="
echo "                           PREDICTION RESULTS"
echo "==============================================================================="
echo ""
echo "  ┌─────────────────────────────────────────────────────────────────────┐"
echo "  │                                                                     │"
printf "  │              PREDICTED BLOOD TYPE:  %-10s                    │\n" "$FULL_TYPE"
echo "  │                                                                     │"
echo "  └─────────────────────────────────────────────────────────────────────┘"
echo ""
echo "  Detailed Analysis:"
echo "  ───────────────────"
echo "  ABO Type:      $ABO_TYPE"
echo "  ABO Genotype:  $ABO_GENOTYPE"
echo "  Rh Factor:     $RH_TYPE"
echo ""

# Blood type information
echo "==============================================================================="
echo "                        BLOOD TYPE INFORMATION"
echo "==============================================================================="
echo ""

case "$ABO_TYPE" in
    "O")
        echo "  🩸 TYPE O - 'Universal Donor'"
        echo ""
        echo "  • Most common blood type globally (~45%)"
        echo "  • Red cells can be given to anyone (O-)"
        echo "  • Can only receive O blood"
        echo "  • Lower risk of heart disease"
        echo "  • Higher risk of ulcers (H. pylori)"
        ;;
    "A"*)
        echo "  🩸 TYPE A"
        echo ""
        echo "  • Second most common (~40%)"
        echo "  • Can donate to A and AB"
        echo "  • Can receive from A and O"
        echo "  • Slightly higher risk of heart disease"
        echo "  • Better adapted to vegetarian diet (historical)"
        ;;
    "B"*)
        echo "  🩸 TYPE B"
        echo ""
        echo "  • Less common (~11%)"
        echo "  • Can donate to B and AB"
        echo "  • Can receive from B and O"
        echo "  • More common in Asian populations"
        echo "  • May have stronger immune system"
        ;;
    "AB")
        echo "  🩸 TYPE AB - 'Universal Recipient'"
        echo ""
        echo "  • Rarest type (~4%)"
        echo "  • Can receive from all blood types"
        echo "  • Can only donate to AB"
        echo "  • Universal plasma donor"
        echo "  • Higher risk of cognitive issues with age"
        ;;
    *)
        echo "  🩸 BLOOD TYPE COULD NOT BE DETERMINED"
        echo ""
        echo "  Insufficient genetic data for reliable prediction."
        ;;
esac

echo ""

if [[ "$RH_TYPE" == "Negative (-)" ]]; then
    echo "  ⚠️  Rh NEGATIVE CONSIDERATIONS:"
    echo "  • About 15% of population is Rh-"
    echo "  • Important for pregnancy (Rh incompatibility)"
    echo "  • Rh- mothers with Rh+ babies may need RhoGAM"
    echo "  • Can only receive Rh- blood"
    echo ""
fi

# Compatibility chart
echo "==============================================================================="
echo "                       TRANSFUSION COMPATIBILITY"
echo "==============================================================================="
echo ""
echo "  Red Blood Cell Compatibility (who can receive from whom):"
echo ""
echo "  Recipient  │ Can Receive From"
echo "  ───────────┼─────────────────────────────────"
echo "   O-        │ O- only"
echo "   O+        │ O-, O+"
echo "   A-        │ O-, A-"
echo "   A+        │ O-, O+, A-, A+"
echo "   B-        │ O-, B-"
echo "   B+        │ O-, O+, B-, B+"
echo "   AB-       │ O-, A-, B-, AB-"
echo "   AB+       │ All types (universal recipient)"
echo ""

# Population frequencies
echo "==============================================================================="
echo "                      POPULATION FREQUENCIES"
echo "==============================================================================="
echo ""
echo "  Blood Type  │ Caucasian │ African │ Hispanic │ Asian"
echo "  ────────────┼───────────┼─────────┼──────────┼───────"
echo "   O+         │   37%     │   47%   │   53%    │  39%"
echo "   O-         │    8%     │    4%   │    4%    │   1%"
echo "   A+         │   33%     │   24%   │   29%    │  27%"
echo "   A-         │    7%     │    2%   │    2%    │   0.5%"
echo "   B+         │    9%     │   18%   │    9%    │  25%"
echo "   B-         │    2%     │    1%   │    1%    │   0.4%"
echo "   AB+        │    3%     │    4%   │    2%    │   7%"
echo "   AB-        │    1%     │    0.3% │   0.2%   │   0.1%"
echo ""

# Calculate confidence
CONFIDENCE="Low"
if [[ "$ABO_TYPE" != "Unknown" && "$RH_TYPE" != "Unknown" ]]; then
    CONFIDENCE="Moderate"
    if [[ $O_ALLELE -eq 2 || ($A_ALLELE -ge 1 && $B_ALLELE -eq 0) || $B_ALLELE -ge 1 ]]; then
        CONFIDENCE="High"
    fi
fi

echo "==============================================================================="
echo "                          PREDICTION CONFIDENCE"
echo "==============================================================================="
echo ""
echo "  Confidence Level: $CONFIDENCE"
echo ""
echo "  Factors affecting accuracy:"
echo "  • Quality of variant calls in VCF"
echo "  • Coverage of ABO gene region"
echo "  • Presence of rare ABO variants"
echo "  • RHD gene deletion detection (limited in SNP data)"
echo ""

echo "==============================================================================="
echo "                              IMPORTANT"
echo "==============================================================================="
echo ""
echo "  ⚠️  This prediction is based on common genetic variants only."
echo ""
echo "  • ALWAYS verify blood type through laboratory testing"
echo "  • NEVER rely on genetic prediction for transfusion decisions"
echo "  • Rare blood types and subtypes cannot be detected"
echo "  • Some variants may not be captured in your VCF file"
echo ""
echo "==============================================================================="
echo "  ⚠️  NOT FOR CLINICAL USE - EDUCATIONAL PURPOSES ONLY"
echo "==============================================================================="
