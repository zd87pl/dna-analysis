#!/bin/bash
#===============================================================================
# BLOOD TYPE PREDICTION
# Predicts ABO and Rh blood type from genetic variants
#
# DISCLAIMER: FOR EDUCATIONAL PURPOSES ONLY
# Always verify blood type through clinical laboratory testing
#===============================================================================

set -euo pipefail

# Shared VCF helpers (hx_prepare_vcf, hx_gt_allele, hx_gt_indel, ...)
# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"

# Validate input
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <vcf_file>"
    echo "Predicts ABO and Rh blood type from genetic variants"
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
# 1 when at least one ABO marker record was present in the VCF (any call);
# 0 when the ABO call rests on absent records only
ABO_RECORDS=0
# 1 when at least one ABO marker carried a non-reference allele
ABO_ALT_EVIDENCE=0

# Every lookup below verifies the record's POS/REF/ALT against the expected
# marker (hx_gt_allele) so an overlapping indel or a different variant at the
# same coordinate can no longer be scored as the ABO/RHD allele.
# Coordinates: GRCh38, see scripts/lib/GRCh38_COORDINATES.md.

# Dosage of the expected ALT allele (0/1/2) from "DOSAGE<TAB>GT"; empty when
# there is no matching, called record.
marker_dosage() { # REGION REF ALT RSID
    local res
    res=$(hx_gt_allele "$1" "$2" "$3" "$4")
    printf '%s\n' "${res%%$'\t'*}"
}

# Check ABO variants

# rs8176719 - 261delG (O allele marker), a 1-bp indel. The GRCh37/GRCh38
# reference sequence carries the deletion (the O allele), so a real VCF
# represents the marker as a 1-bp INSERTION (9:133257521 T>TC, ALT = non-O,
# 0/0 = O/O, 0/1 = X/O). The older deletion form (9:133257520 TC>T, ALT = O)
# is accepted too. Matched by POS and a 1-bp REF/ALT length difference.
# TODO(build): unverified GRCh38 orientation/POS (dbSNP not reachable offline)
O_KNOWN=0  # 1 once the rs8176719 record has been read
O_RESULT=$(hx_gt_indel "9:133257520-133257521" 1)
if [[ -n "$O_RESULT" ]]; then
    read -r O_POS O_REF O_ALT O_GT <<< "$O_RESULT"
    if [[ "$O_POS" == "133257520" || "$O_POS" == "133257521" ]]; then
        # Index (1-based) of the ALT allele that is REF minus one base (O) or
        # REF plus one base (non-O; then allele 0 is O)
        O_IDX=0; NONO_IDX=0; i=0
        IFS=',' read -r -a O_ALTS <<< "$O_ALT"
        for a in "${O_ALTS[@]}"; do
            i=$((i + 1))
            if [[ ${#a} -eq $(( ${#O_REF} - 1 )) && "$O_REF" == "$a"* ]]; then
                O_IDX=$i
                break
            elif [[ ${#a} -eq $(( ${#O_REF} + 1 )) && "$a" == "$O_REF"* ]]; then
                NONO_IDX=$i
                break
            fi
        done
        if [[ $O_IDX -gt 0 || $NONO_IDX -gt 0 ]]; then
            ABO_RECORDS=1
            O_KNOWN=1
            IFS='/' read -r -a O_GTA <<< "$O_GT"
            for a in "${O_GTA[@]}"; do
                if [[ $O_IDX -gt 0 ]]; then
                    [[ "$a" == "$O_IDX" ]] && O_ALLELE=$((O_ALLELE + 1))
                else
                    [[ "$a" == "0" ]] && O_ALLELE=$((O_ALLELE + 1))
                fi
            done
            # A non-reference call at the marker is direct evidence (deletion
            # form: an O allele; insertion form: a non-O allele). A hom-ref
            # insertion-form call (O/O) is reference-only and keeps the
            # confidence capped at Moderate.
            [[ "$O_GT" != "0/0" ]] && ABO_ALT_EVIDENCE=1
        else
            echo "  ℹ️  rs8176719: record at 9:$O_POS is ${O_REF}>${O_ALT}, neither the 261delG deletion nor its insertion form - skipped"
        fi
    fi
fi

# rs8176746 (c.796C>A on the minus strand; plus-strand VCF G>T) - B allele
# marker: non-B (A/O-like) = G (REF), B = T (ALT) at this position
# TODO(build): unverified GRCh38 position (derived, dbSNP not consulted)
B_DOSAGE=$(marker_dosage "9:133256264-133256264" G T rs8176746)
if [[ -n "$B_DOSAGE" ]]; then
    ABO_RECORDS=1
    B_ALLELE=$B_DOSAGE
    [[ $B_ALLELE -gt 0 ]] && ABO_ALT_EVIDENCE=1
fi

# rs8176747 (803G>C) is a second B marker, but its GRCh38 position/alleles
# are not in scripts/lib/GRCh38_COORDINATES.md, so it is not evaluated.

# rs7853989 - A1 vs A2 distinction and A/B (REF G associated with A)
A_DOSAGE=$(marker_dosage "9:133256042-133256042" G C rs7853989)
if [[ -n "$A_DOSAGE" ]]; then
    ABO_RECORDS=1
    if [[ $A_DOSAGE -eq 0 ]]; then
        # Reference allele associated with A; the A count follows from the O
        # count only when the rs8176719 record was actually read (its absence
        # means reference = O, not A).
        if [[ $B_ALLELE -eq 0 && $O_KNOWN -eq 1 ]]; then
            A_ALLELE=$((2 - O_ALLELE))
        fi
    else
        ABO_ALT_EVIDENCE=1
    fi
fi

# rs7853989 absent from a variants-only VCF = reference (A-associated); the
# insertion-form rs8176719 call already fixes the number of non-O alleles.
A_ASSUMED=0
if [[ $O_KNOWN -eq 1 && $B_ALLELE -eq 0 && $A_ALLELE -eq 0 && -z "$A_DOSAGE" ]]; then
    A_ALLELE=$((2 - O_ALLELE))
    [[ $A_ALLELE -gt 0 ]] && A_ASSUMED=1
fi

# Check Rh factor

# rs590787 - Proxy for RHD (not perfect but commonly used); ALT carrier = Rh+
RH_DOSAGE=$(marker_dosage "1:25362249-25362249" G A rs590787)
if [[ -n "$RH_DOSAGE" ]]; then
    RH_UNKNOWN=0
    if [[ $RH_DOSAGE -eq 0 ]]; then
        RH_POSITIVE=0  # Likely Rh-
    else
        RH_POSITIVE=1  # Likely Rh+
    fi
fi

# Additional RHD check - rs676785 (homozygous ALT = Rh+)
RH_DOSAGE2=$(marker_dosage "1:25409008-25409008" C T rs676785)
if [[ -n "$RH_DOSAGE2" && $RH_UNKNOWN -eq 1 ]]; then
    RH_UNKNOWN=0
    if [[ $RH_DOSAGE2 -eq 2 ]]; then
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
    # No O/A/B allele observed and no rs8176719 record. Absence of that record
    # cannot be read as "non-O": the reference carries the O (261delG) allele
    # and a variants-only VCF omits hom-ref sites, so A and O are both open.
    ABO_TYPE="A or O"
    ABO_GENOTYPE="A or O (rs8176719 not in VCF)"
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
if [[ "$ABO_TYPE" != "Unknown" && "$ABO_TYPE" != "A or O" && "$RH_TYPE" != "Unknown" ]]; then
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
echo "  ✅ Predicted blood type: $FULL_TYPE (ABO genotype $ABO_GENOTYPE, Rh $RH_TYPE)"
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
        ;;
    "AB")
        # Must precede "A"*, which would otherwise shadow this arm
        echo "  🩸 TYPE AB - 'Universal Recipient'"
        echo ""
        echo "  • Rarest type (~4%)"
        echo "  • Can receive from all blood types"
        echo "  • Can only donate to AB"
        echo "  • Universal plasma donor"
        ;;
    "A"*)
        echo "  🩸 TYPE A"
        echo ""
        echo "  • Second most common (~40%)"
        echo "  • Can donate to A and AB"
        echo "  • Can receive from A and O"
        ;;
    "B"*)
        echo "  🩸 TYPE B"
        echo ""
        echo "  • Less common (~11%)"
        echo "  • Can donate to B and AB"
        echo "  • Can receive from B and O"
        echo "  • More common in Asian populations"
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

# Calculate confidence. "High" requires at least one ABO marker with a
# non-reference call; a call that rests only on homozygous-reference or
# absent records is capped at "Moderate".
CONFIDENCE="Low"
CONFIDENCE_NOTE=""
if [[ "$ABO_TYPE" == "A or O" ]]; then
    CONFIDENCE_NOTE="rs8176719 (261delG) record absent: A and O cannot be told apart"
elif [[ "$ABO_TYPE" != "Unknown" && "$RH_TYPE" != "Unknown" ]]; then
    CONFIDENCE="Moderate"
    if [[ $O_ALLELE -eq 2 || ($A_ALLELE -ge 1 && $B_ALLELE -eq 0) || $B_ALLELE -ge 1 ]]; then
        CONFIDENCE="High"
    fi
    if [[ $A_ASSUMED -eq 1 && "$CONFIDENCE" == "High" ]]; then
        CONFIDENCE="Moderate"
        CONFIDENCE_NOTE="rs7853989 record absent (reference assumed)"
    fi
    if [[ $ABO_ALT_EVIDENCE -eq 0 && "$CONFIDENCE" == "High" ]]; then
        CONFIDENCE="Moderate"
        if [[ $ABO_RECORDS -eq 1 ]]; then
            CONFIDENCE_NOTE="ABO call rests on homozygous-reference calls only"
        else
            CONFIDENCE_NOTE="ABO call rests on absent marker records only"
        fi
    fi
fi

echo "==============================================================================="
echo "                          PREDICTION CONFIDENCE"
echo "==============================================================================="
echo ""
echo "  Confidence Level: $CONFIDENCE"
if [[ -n "$CONFIDENCE_NOTE" ]]; then
    echo "  ($CONFIDENCE_NOTE)"
fi
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
# Marker coverage (found/queried); the frontend treats 0 found as "no data".
hx_report_coverage
