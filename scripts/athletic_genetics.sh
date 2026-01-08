#!/bin/bash
#===============================================================================
# EXTENDED SPORTS GENETICS PROFILE
# Comprehensive athletic potential analysis
#===============================================================================

set -euo pipefail

#===============================================================================
# CONFIGURATION - SNP positions (GRCh38)
#===============================================================================

SCRIPT_NAME="$(basename "$0")"

# Genomic coordinates for athletic variants
declare -A SNP_POSITIONS=(
    # Power/Strength
    ["ACTN3"]="11:66560624"
    ["AGT"]="1:230845794"
    ["IL6"]="7:22727026"
    ["MSTN"]="2:190379711"
    # Endurance
    ["ACE"]="17:63488529"
    ["PPARGC1A"]="4:23814039"
    ["VEGFA"]="6:43770613"
    ["NFE2L2"]="2:178095031"
    ["PPARA"]="22:46615880"
    # Recovery
    ["CRP"]="1:159682233"
    ["TNF"]="6:31575254"
    ["SOD2"]="6:160113872"
    ["COL1A1"]="17:50201632"
    # Injury
    ["COL5A1"]="9:137684151"
    ["GDF5"]="20:34025756"
    ["MMP3"]="11:102731092"
    # Psychology
    ["COMT"]="22:19951271"
    ["DRD4"]="11:637269"
    ["BDNF"]="11:27679916"
)

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

error_exit() {
    echo -e "${RED}Error: $1${NC}" >&2
    exit 1
}

# Check for required dependencies
check_dependencies() {
    local missing=()
    for cmd in bcftools bc; do
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

# Query SNP genotype
get_genotype() {
    local pos="$1"
    bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1
}

# Query SNP with REF/ALT
get_snp_info() {
    local pos="$1"
    bcftools query -r "$pos" -f '%REF %ALT [%GT]\n' "$VCF" 2>/dev/null | head -1
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

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║           EXTENDED SPORTS GENETICS PROFILE                                ║"
echo "║           Comprehensive athletic potential analysis                       ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

#===============================================================================
# CATEGORY 1: POWER AND STRENGTH
#===============================================================================

echo "💪 POWER AND STRENGTH"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

POWER_SCORE=0

# ACTN3 - key strength gene
result=$(get_snp_info "${SNP_POSITIONS[ACTN3]}")
GT=$(echo "$result" | awk '{print $3}')

if [[ -n "$GT" ]]; then
    case $GT in
        "0/0"|"0|0")
            echo "  ACTN3 (rs1815739): R/R → 🏆 ELITE power genetics!"
            echo "    • Full alpha-actinin-3 function in fast-twitch fibers"
            echo "    • Advantage in: sprint, jumps, throws, weightlifting"
            echo "    • ~18% of population has this genotype"
            POWER_SCORE=$((POWER_SCORE + 3))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  ACTN3 (rs1815739): R/X → Good power/endurance mix"
            echo "    • Good for mixed sports (soccer, basketball)"
            POWER_SCORE=$((POWER_SCORE + 2))
            ;;
        "1/1"|"1|1")
            echo "  ACTN3 (rs1815739): X/X → Endurance genetics"
            echo "    • Better performance in endurance sports"
            echo "    • ~18% of population has this genotype"
            POWER_SCORE=$((POWER_SCORE + 0))
            ;;
    esac
fi

# AGT - angiotensinogen (muscle strength)
result=$(get_genotype "${SNP_POSITIONS[AGT]}")
if [[ -n "$result" ]]; then
    case $result in
        "1/1"|"1|1")
            echo "  AGT (rs699): T/T → Higher muscle strength"
            POWER_SCORE=$((POWER_SCORE + 1))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  AGT (rs699): M/T → Average predisposition"
            ;;
    esac
fi

# IL6 - strength training response
result=$(get_genotype "${SNP_POSITIONS[IL6]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  IL6 (rs1800795): G/G → Better strength training response"
            POWER_SCORE=$((POWER_SCORE + 1))
            ;;
    esac
fi

# MSTN - myostatin (muscle growth brake)
result=$(get_genotype "${SNP_POSITIONS[MSTN]}")
if [[ -n "$result" ]] && [[ "$result" != "0/0" ]] && [[ "$result" != "0|0" ]]; then
    echo "  MSTN: Variant → Potentially better muscle hypertrophy"
    POWER_SCORE=$((POWER_SCORE + 1))
fi

echo ""
echo "  POWER SCORE: $POWER_SCORE / 6"
echo ""

#===============================================================================
# CATEGORY 2: ENDURANCE
#===============================================================================

echo "🏃 ENDURANCE"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

ENDURANCE_SCORE=0

# ACE - angiotensin-converting enzyme
result=$(get_genotype "${SNP_POSITIONS[ACE]}")
if [[ -n "$result" ]]; then
    case $result in
        "1/1"|"1|1")
            echo "  ACE (I/D): I/I → Better endurance, aerobic capacity"
            echo "    • Advantage in: long runs, cycling, swimming"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  ACE (I/D): I/D → Balanced profile"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
        "0/0"|"0|0")
            echo "  ACE (I/D): D/D → Better strength and power"
            ;;
    esac
fi

# PPARGC1A (PGC-1alpha) - mitochondrial biogenesis
result=$(get_genotype "${SNP_POSITIONS[PPARGC1A]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  PPARGC1A (rs8192678): G/G → Better endurance adaptation"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  PPARGC1A (rs8192678): G/A → Average response"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
    esac
fi

# VEGFA - angiogenesis
result=$(get_genotype "${SNP_POSITIONS[VEGFA]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  VEGFA (rs2010963): G/G → Better angiogenesis"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
    esac
fi

# NRF2 - antioxidant response
result=$(get_genotype "${SNP_POSITIONS[NFE2L2]}")
if [[ -n "$result" ]] && [[ "$result" != "0/0" ]] && [[ "$result" != "0|0" ]]; then
    echo "  NFE2L2: Variant → Better antioxidant protection"
    ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
fi

# PPARA - fat metabolism as fuel
result=$(get_genotype "${SNP_POSITIONS[PPARA]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  PPARA (rs4253778): G/G → Better fat burning"
            ENDURANCE_SCORE=$((ENDURANCE_SCORE + 1))
            ;;
    esac
fi

echo ""
echo "  ENDURANCE SCORE: $ENDURANCE_SCORE / 7"
echo ""

#===============================================================================
# CATEGORY 3: RECOVERY
#===============================================================================

echo "🔄 RECOVERY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

RECOVERY_SCORE=0

# IL6 - post-training inflammation
result=$(get_genotype "${SNP_POSITIONS[IL6]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  IL6 (rs1800795): G/G → Faster recovery, less inflammation"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
        "1/1"|"1|1")
            echo "  IL6 (rs1800795): C/C → Slower recovery, more DOMS"
            ;;
    esac
fi

# CRP - C-reactive protein
result=$(get_genotype "${SNP_POSITIONS[CRP]}")
if [[ -n "$result" ]] && [[ "$result" == "0/0" || "$result" == "0|0" ]]; then
    echo "  CRP (rs1205): Lower inflammation level"
    RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
fi

# TNF - tumor necrosis factor
result=$(get_genotype "${SNP_POSITIONS[TNF]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  TNF (rs1800629): G/G → Lower inflammatory response"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
    esac
fi

# SOD2 - superoxide dismutase
result=$(get_genotype "${SNP_POSITIONS[SOD2]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  SOD2 (rs4880): Heterozygote → Optimal antioxidant protection"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 2))
            ;;
        "0/0"|"1/1"|"0|0"|"1|1")
            echo "  SOD2 (rs4880): Homozygote → Less optimal protection"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
    esac
fi

# COL1A1 - collagen (tendon recovery)
result=$(get_genotype "${SNP_POSITIONS[COL1A1]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  COL1A1 (rs1800012): G/G → Better collagen, fewer injuries"
            RECOVERY_SCORE=$((RECOVERY_SCORE + 1))
            ;;
    esac
fi

echo ""
echo "  RECOVERY SCORE: $RECOVERY_SCORE / 7"
echo ""

#===============================================================================
# CATEGORY 4: INJURY RISK
#===============================================================================

echo "🩹 INJURY RISK"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

INJURY_RISK=0

# COL5A1 - ACL and tendon injuries
result=$(get_genotype "${SNP_POSITIONS[COL5A1]}")
if [[ -n "$result" ]]; then
    case $result in
        "1/1"|"1|1")
            echo "  COL5A1 (rs12722): T/T → Elevated ACL injury risk"
            INJURY_RISK=$((INJURY_RISK + 2))
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  COL5A1 (rs12722): C/T → Moderate injury risk"
            INJURY_RISK=$((INJURY_RISK + 1))
            ;;
        "0/0"|"0|0")
            echo "  COL5A1 (rs12722): C/C → Lower tendon injury risk"
            ;;
    esac
fi

# GDF5 - cartilage development
result=$(get_genotype "${SNP_POSITIONS[GDF5]}")
if [[ -n "$result" ]]; then
    case $result in
        "1/1"|"1|1")
            echo "  GDF5 (rs143383): A/A → Elevated osteoarthritis risk"
            INJURY_RISK=$((INJURY_RISK + 1))
            ;;
    esac
fi

# MMP3 - matrix degradation
result=$(get_genotype "${SNP_POSITIONS[MMP3]}")
if [[ -n "$result" ]] && [[ "$result" != "0/0" ]] && [[ "$result" != "0|0" ]]; then
    echo "  MMP3: Variant → Higher Achilles tendon injury risk"
    INJURY_RISK=$((INJURY_RISK + 1))
fi

echo ""
if [[ $INJURY_RISK -le 1 ]]; then
    echo -e "  INJURY RISK: ${GREEN}LOW ✅${NC}"
elif [[ $INJURY_RISK -le 2 ]]; then
    echo -e "  INJURY RISK: ${YELLOW}MODERATE 🟡${NC}"
else
    echo -e "  INJURY RISK: ${RED}ELEVATED 🔴${NC}"
    echo "  → Recommendation: Longer warm-ups, work on flexibility"
fi
echo ""

#===============================================================================
# CATEGORY 5: SPORTS PSYCHOLOGY
#===============================================================================

echo "🧠 SPORTS PSYCHOLOGY"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# COMT - competition stress resilience
result=$(get_genotype "${SNP_POSITIONS[COMT]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  COMT (rs4680): Val/Val → 'Warrior'"
            echo "    • Better performance under pressure"
            echo "    • Lower baseline dopamine = better during stress"
            echo "    • Good for: competitions, matches, important events"
            ;;
        "0/1"|"0|1"|"1/0"|"1|0")
            echo "  COMT (rs4680): Val/Met → Balanced profile"
            echo "    • Good cognitive flexibility"
            ;;
        "1/1"|"1|1")
            echo "  COMT (rs4680): Met/Met → 'Worrier'"
            echo "    • Better focus and precision in calm conditions"
            echo "    • May need stress management techniques before competitions"
            echo "    • Good for: technical training, precision sports"
            ;;
    esac
fi

# DRD4 - sensation seeking / risk
result=$(get_genotype "${SNP_POSITIONS[DRD4]}")
if [[ -n "$result" ]] && [[ "$result" != "0/0" ]] && [[ "$result" != "0|0" ]]; then
    echo "  DRD4: 7R Variant → Higher sensation seeking"
    echo "    • Good for: extreme sports, high risk"
fi

# BDNF - motor learning
result=$(get_genotype "${SNP_POSITIONS[BDNF]}")
if [[ -n "$result" ]]; then
    case $result in
        "0/0"|"0|0")
            echo "  BDNF (rs6265): Val/Val → Faster motor learning"
            ;;
        "0/1"|"1/1"|"0|1"|"1|1")
            echo "  BDNF (rs6265): Met carrier → Better long-term motor memory"
            ;;
    esac
fi

echo ""

#===============================================================================
# PROFILE SUMMARY
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    SPORTS PROFILE SUMMARY                                 ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Calculate total scores
POWER_PCT=$(echo "scale=0; $POWER_SCORE * 100 / 6" | bc)
ENDURANCE_PCT=$(echo "scale=0; $ENDURANCE_SCORE * 100 / 7" | bc)
RECOVERY_PCT=$(echo "scale=0; $RECOVERY_SCORE * 100 / 7" | bc)

echo "SCORES:"
echo ""
printf "  Power:      "
for ((i=1; i<=10; i++)); do
    if [[ $i -le $((POWER_PCT / 10)) ]]; then
        printf "█"
    else
        printf "░"
    fi
done
echo " ${POWER_PCT}%"

printf "  Endurance:  "
for ((i=1; i<=10; i++)); do
    if [[ $i -le $((ENDURANCE_PCT / 10)) ]]; then
        printf "█"
    else
        printf "░"
    fi
done
echo " ${ENDURANCE_PCT}%"

printf "  Recovery:   "
for ((i=1; i<=10; i++)); do
    if [[ $i -le $((RECOVERY_PCT / 10)) ]]; then
        printf "█"
    else
        printf "░"
    fi
done
echo " ${RECOVERY_PCT}%"

echo ""

# Sports recommendations
echo "RECOMMENDED DISCIPLINES:"
echo ""

if [[ $POWER_SCORE -ge 4 ]] && [[ $ENDURANCE_SCORE -le 3 ]]; then
    echo "  🥇 POWER AND SPEED SPORTS:"
    echo "     • Sprint (100m, 200m)"
    echo "     • Weightlifting"
    echo "     • Jumps"
    echo "     • Throws"
    echo "     • CrossFit (strength components)"
    echo "     • American football"
elif [[ $ENDURANCE_SCORE -ge 4 ]] && [[ $POWER_SCORE -le 3 ]]; then
    echo "  🥇 ENDURANCE SPORTS:"
    echo "     • Long distance running (5km+)"
    echo "     • Road cycling"
    echo "     • Long distance swimming"
    echo "     • Triathlon"
    echo "     • Mountain running"
elif [[ $POWER_SCORE -ge 3 ]] && [[ $ENDURANCE_SCORE -ge 3 ]]; then
    echo "  🥇 MIXED SPORTS (POWER + ENDURANCE):"
    echo "     • Soccer"
    echo "     • Basketball"
    echo "     • Tennis"
    echo "     • MMA / Judo / Wrestling"
    echo "     • CrossFit"
    echo "     • Middle distance running (800m-1500m)"
fi

echo ""
echo "TRAINING TIPS:"
echo ""

if [[ $RECOVERY_SCORE -le 3 ]]; then
    echo "  ⚠️ Recovery: Slower - plan more rest days"
    echo "     • 48-72h between intense workouts"
    echo "     • Recovery techniques: sleep, massage, sauna"
fi

if [[ $INJURY_RISK -ge 2 ]]; then
    echo "  ⚠️ Injuries: Elevated risk - be careful!"
    echo "     • Longer warm-ups (15-20 min)"
    echo "     • Work on flexibility and stability"
    echo "     • Avoid sudden load increases"
fi

echo ""
