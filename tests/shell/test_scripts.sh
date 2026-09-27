#!/usr/bin/env bash
# =============================================================================
# tests/shell/test_scripts.sh - regression tests for individual analysis
# scripts on small synthetic GRCh38 VCFs (no bats). Exit status is non-zero on
# any failure.
#
# Usage: bash tests/shell/test_scripts.sh
# Requires: bash >= 4, bcftools. Fixtures are built under a temporary
# directory from data/sample_genome.vcf and removed on exit; nothing is
# written into the repository (scripts run from the temporary directory with
# HELIXIGHT_WORK_DIR set).
# =============================================================================
set -uo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SCRIPTS="$ROOT/scripts"
SAMPLE="$ROOT/data/sample_genome.vcf"

for tool in bcftools; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "SKIP: $tool not installed" >&2
        exit 0
    fi
done
[[ -f $SAMPLE ]] || { echo "FAIL: sample not found: $SAMPLE"; exit 1; }

TMP=$(mktemp -d "${TMPDIR:-/tmp}/hx_scripts_test.XXXXXX")
trap 'rm -rf "$TMP"' EXIT
export HELIXIGHT_WORK_DIR="$TMP/work"
cd "$TMP" || exit 1

PASS=0
FAIL=0
pass() { PASS=$((PASS + 1)); echo "PASS: $1"; }
fail() { FAIL=$((FAIL + 1)); echo "FAIL: $1${2:+ -- $2}"; }
assert_eq() { # NAME EXPECTED ACTUAL
    if [[ $2 == "$3" ]]; then pass "$1"; else fail "$1" "expected [$2] got [$3]"; fi
}
assert_contains() { # NAME NEEDLE HAYSTACK
    if [[ $3 == *"$2"* ]]; then pass "$1"; else fail "$1" "[$2] not found in output"; fi
}
assert_not_contains() { # NAME NEEDLE HAYSTACK
    if [[ $3 != *"$2"* ]]; then pass "$1"; else fail "$1" "[$2] unexpectedly found in output"; fi
}

# mk FILE RECORD... : a VCF with the sample's header and the given records
# (tab-separated, written with '|' for readability).
mk() {
    local out=$1; shift
    {
        grep '^##' "$SAMPLE"
        grep '^#CHROM' "$SAMPLE"
        local r
        for r in "$@"; do printf '%s\n' "$r" | tr '|' '\t'; done
    } >"$out"
}
# sample_with ID GT [POS REF ALT]: the sample with one record's GT (and
# optionally POS/REF/ALT) rewritten.
sample_with() {
    local id=$1 gt=$2 pos=${3:-} ref=${4:-} alt=${5:-}
    awk -v id="$id" -v gt="$gt" -v pos="$pos" -v ref="$ref" -v alt="$alt" 'BEGIN{OFS="\t"}
        $3==id { $10=gt; if (pos!="") { $2=pos; $4=ref; $5=alt } } { print }' "$SAMPLE"
}
# run SCRIPT VCF: stdout -> $OUT, stderr -> $ERR, exit status -> $RC
run() {
    local d="$TMP/run.$1" vcf=$2
    [[ $vcf == /* ]] || vcf="$TMP/$vcf"
    rm -rf "$d"; mkdir -p "$d"
    OUT=$(cd "$d" && bash "$SCRIPTS/$1" "$vcf" 2>"$TMP/stderr.txt")
    RC=$?
    ERR=$(cat "$TMP/stderr.txt")
}

# ---------------------------------------------------------------------------
# cancer_genes_scan.sh: point-mutation probes match alleles / variant type
# ---------------------------------------------------------------------------
mk benign.vcf '17|7676154|rs1042522|G|C|100|PASS|.|GT|1/1'
run cancer_genes_scan.sh benign.vcf
assert_eq "cancer: rc 0 on a benign-only VCF" 0 "$RC"
assert_not_contains "cancer: rs1042522 (Pro72Arg) is not reported as a TP53 hotspot" "🔴" "$OUT"
assert_contains "cancer: benign-only VCF -> TOTAL 0" "TOTAL:             0" "$OUT"

mk hotspot.vcf '17|7675088|rs28934578|C|T|100|PASS|.|GT|0/1'
run cancer_genes_scan.sh hotspot.vcf
assert_contains "cancer: TP53 R175H hotspot (C>T, rs28934578) is reported" "🔴 TP53 R175H: 0/1" "$OUT"

mk snp_at_indel.vcf '17|43124027|.|A|G|100|PASS|.|GT|0/1'
run cancer_genes_scan.sh snp_at_indel.vcf
assert_not_contains "cancer: a SNP at the 185delAG position is not reported as the deletion" "185delAG" "$OUT"

# Indels at the (unverified, derived) BRCA positions are informational only
mk brca_indels.vcf '13|32340300|.|GT|G|100|PASS|.|GT|0/1' '17|43057062|.|T|TG|100|PASS|.|GT|0/1'
run cancer_genes_scan.sh brca_indels.vcf
assert_eq "cancer: rc 0 with indels at listed BRCA positions" 0 "$RC"
assert_contains "cancer: BRCA2 6174delT position hit is an informational note" "ℹ️  BRCA2: variant at the listed 6174delT position" "$OUT"
assert_contains "cancer: BRCA1 5382insC position hit is an informational note" "ℹ️  BRCA1: variant at the listed 5382insC position" "$OUT"
assert_not_contains "cancer: listed-position hits are not red findings" "🔴" "$OUT"
assert_contains "cancer: listed-position hits do not count as pathogenic" "TOTAL:             0" "$OUT"
assert_not_contains "cancer: no 'FOUND N pathogenic' for listed-position hits" "pathogenic variants!" "$OUT"

# A region full of explicit hom-ref records (ALT ".") counts as covered
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {$5="."; $10="0/0"; print}' "$SAMPLE" >homref_all.vcf
run cancer_genes_scan.sh homref_all.vcf
assert_eq "cancer: rc 0 on an all-hom-ref VCF" 0 "$RC"
assert_contains "cancer: hom-ref record inside BRCA2 counts as data (not 0/30)" "HX_DATA_COVERAGE: 1/30 markers" "$OUT"

mk cdh1.vcf '16|68737400|.|A|G|100|PASS|.|GT|0/1' '16|68860000|.|A|G|100|PASS|.|GT|0/1'
run cancer_genes_scan.sh cdh1.vcf
assert_contains "cancer: CDH1 counted in its GRCh38 window (1 of 2 records)" "CDH1       |        1 |" "$OUT"

# ---------------------------------------------------------------------------
# carrier_screening.sh: insertion (HEXA) and delins (BLM) rows are scored
# ---------------------------------------------------------------------------
{
    cat "$SAMPLE"
    printf '15\t72346579\trs76173977\tT\tTTATC\t100\tPASS\t.\tGT\t0/1\n'
    printf '15\t90631934\trs113993962\tATCTGA\tTAGATTC\t100\tPASS\t.\tGT\t0/1\n'
} >ashkenazi.vcf
run carrier_screening.sh ashkenazi.vcf
assert_eq "carrier: rc 0" 0 "$RC"
assert_contains "carrier: Tay-Sachs TATC insertion scored" "Gene: HEXA (rs76173977)" "$OUT"
assert_contains "carrier: Bloom delins scored" "Gene: BLM (rs113993962)" "$OUT"

# HBB rows next to rs334: HbC (c.19G>A, 11:5227003) and IVS-I-1 (11:5226929)
mk hbb.vcf '11|5226929|rs33971440|C|T|100|PASS|.|GT|0/1' '11|5227003|rs33930165|C|T|100|PASS|.|GT|0/1'
run carrier_screening.sh hbb.vcf
assert_eq "carrier: rc 0 on HBB fixture" 0 "$RC"
assert_contains "carrier: HbC rs33930165 at 11:5227003 scored" "Gene: HBB (rs33930165)" "$OUT"
assert_contains "carrier: IVS-I-1 rs33971440 at 11:5226929 scored" "Gene: HBB (rs33971440)" "$OUT"

# ---------------------------------------------------------------------------
# mega_analysis.sh: indel F508del and PAH at the GRCh38 position
# ---------------------------------------------------------------------------
sample_with rs113993960 0/1 >f508.vcf
run mega_analysis.sh f508.vcf
assert_eq "mega: rc 0" 0 "$RC"
assert_contains "mega: CFTR F508del carrier detected (indel record)" "CFTR F508del: NOSICIEL" "$OUT"

sample_with rs5030858 0/1 >pku.vcf
run mega_analysis.sh pku.vcf
assert_contains "mega: PAH rs5030858 carrier detected at 12:102840632" "PAH rs5030858: NOSICIEL" "$OUT"

# Tay-Sachs: the TATC insertion at the carrier_screening.sh position is a
# carrier; a SNP at the old GRCh37 locus (15:72638892) is not.
run mega_analysis.sh ashkenazi.vcf
assert_contains "mega: HEXA TATC insertion (15:72346579) -> carrier" "HEXA: Wariant Tay-Sachsa wykryty" "$OUT"
mk hexa_grch37_snp.vcf '15|72638892|rsX|G|A|100|PASS|.|GT|0/1'
run mega_analysis.sh hexa_grch37_snp.vcf
assert_contains "mega: SNP at the GRCh37 HEXA locus -> no Tay-Sachs call" "HEXA: Nie wykryto" "$OUT"
mk hexa_snp_at_pos.vcf '15|72346579|.|T|G|100|PASS|.|GT|0/1'
run mega_analysis.sh hexa_snp_at_pos.vcf
assert_contains "mega: SNP at the HEXA insertion position -> no Tay-Sachs call" "HEXA: Nie wykryto" "$OUT"

# ---------------------------------------------------------------------------
# quick_insights.sh: MC1R counts rs1805007 alleles at the canonical position
# ---------------------------------------------------------------------------
run quick_insights.sh "$SAMPLE"
assert_eq "quick: rc 0" 0 "$RC"
assert_contains "quick: hom-ref rs1805007 -> no MC1R variant" "Warianty MC1R: 0" "$OUT"
sample_with rs1805007 0/1 >red.vcf
run quick_insights.sh red.vcf
assert_contains "quick: het rs1805007 -> 1 MC1R variant" "Warianty MC1R: 1" "$OUT"

# ---------------------------------------------------------------------------
# athletic_genetics.sh: explicit hom-ref records (ALT ".") are interpreted
# ---------------------------------------------------------------------------
sample_with rs1815739 0/0 66560624 C . >homref.vcf
run athletic_genetics.sh homref.vcf
assert_eq "athletic: rc 0 with ALT-. record" 0 "$RC"
assert_contains "athletic: ACTN3 R/R read from an ALT-. hom-ref record" "ACTN3 (rs1815739): R/R" "$OUT"

# Unverified GRCh37-range probes (COL5A1 9:137684151, NFE2L2 2:178095031) are
# skipped by all three fitness scripts; PPARA rs4253778 G/G is the endurance
# genotype in both athletic_genetics.sh and triathlon_genetics.sh.
mk oldprobes.vcf '2|178095031|rsZ|C|T|100|PASS|.|GT|0/1' '9|137684151|rsY|C|T|100|PASS|.|GT|1/1' '22|46150500|rs4253778|G|C|100|PASS|.|GT|0/0'
run athletic_genetics.sh oldprobes.vcf
assert_eq "athletic: rc 0 on old-probe fixture" 0 "$RC"
assert_not_contains "athletic: variant at 9:137684151 is not an ACL injury call" "Elevated ACL injury risk" "$OUT"
assert_not_contains "athletic: variant at 2:178095031 is not an antioxidant call" "Better antioxidant protection" "$OUT"
assert_contains "athletic: COL5A1 skipped" "COL5A1 (rs12722): position unverified - skipped" "$OUT"
assert_contains "athletic: PPARA G/G -> better fat burning, endurance 1/6" "ENDURANCE SCORE: 1 / 6" "$OUT"
assert_contains "athletic: coverage counts 18 markers" "HX_DATA_COVERAGE: 1/18 markers" "$OUT"
run triathlon_genetics.sh oldprobes.vcf
assert_eq "triathlon: rc 0 on old-probe fixture" 0 "$RC"
assert_not_contains "triathlon: variant at 9:137684151 is not an injury call" "wyższe ryzyko kontuzji" "$OUT"
assert_contains "triathlon: PPARA G/G is the favourable genotype" "✅ PPARA rs4253778: G/G" "$OUT"
assert_contains "triathlon: coverage counts 30 markers" "HX_DATA_COVERAGE: 1/30 markers" "$OUT"
mk ppara_cc.vcf '22|46150500|rs4253778|G|C|100|PASS|.|GT|1/1'
run triathlon_genetics.sh ppara_cc.vcf
assert_contains "triathlon: PPARA C/C is the unfavourable genotype" "🔴 PPARA rs4253778: C/C" "$OUT"
run athletic_genetics.sh ppara_cc.vcf
assert_contains "athletic: PPARA C/C -> no endurance point" "ENDURANCE SCORE: 0 / 6" "$OUT"
run personalized_protocol.sh oldprobes.vcf
assert_eq "protocol: rc 0 on old-probe fixture" 0 "$RC"
assert_not_contains "protocol: variant at 9:137684151 is not a collagen risk call" "Warianty ryzyka kontuzji wykryte" "$OUT"
assert_contains "protocol: coverage counts 26 markers" "HX_DATA_COVERAGE: 0/26 markers" "$OUT"

# ---------------------------------------------------------------------------
# pharmacogenomics.sh: normal-function alleles are informational, not red
# ---------------------------------------------------------------------------
run pharmacogenomics.sh "$SAMPLE"
assert_eq "pgx: rc 0" 0 "$RC"
assert_not_contains "pgx: CYP2D6*2 rs16947 (normal function) is not a red finding" "🔴 rs16947" "$OUT"
assert_contains "pgx: rs16947 listed as informational" "ℹ️  rs16947 (CYP2D6): Common variant, normal function - no dosing impact" "$OUT"
assert_contains "pgx: Variants Found excludes the normal-function allele" "Variants Found:      7" "$OUT"
assert_contains "pgx: coverage still counts rs16947 as data" "HX_DATA_COVERAGE: 20/26 markers" "$OUT"

# ---------------------------------------------------------------------------
# blood_type.sh: absent rs8176719 is not read as a definitive A
# ---------------------------------------------------------------------------
grep -v rs8176719 "$SAMPLE" >noabo.vcf
run blood_type.sh noabo.vcf
assert_eq "blood: rc 0" 0 "$RC"
assert_contains "blood: no rs8176719 record -> A or O, not AA" "A or O (rs8176719 not in VCF)" "$OUT"
assert_contains "blood: no rs8176719 record -> Confidence Low" "Confidence Level: Low" "$OUT"
run blood_type.sh "$SAMPLE"
assert_contains "blood: sample (T>TC 0/1) -> A+ AO" "A+ (ABO genotype AO" "$OUT"
assert_contains "blood: sample -> all 5 markers read" "HX_DATA_COVERAGE: 5/5 markers" "$OUT"
assert_not_contains "blood: sample rs8176746 record (G>T) matches the probe" "expected" "$ERR"

# rs8176746 is minus-strand c.796C>A = plus-strand G>T: a B/B genome
mk bb.vcf '1|25362249|rs590787|G|A|100|PASS|.|GT|0/1' '9|133256042|rs7853989|G|C|100|PASS|.|GT|1/1' '9|133256264|rs8176746|G|T|100|PASS|.|GT|1/1'
run blood_type.sh bb.vcf
assert_eq "blood: rc 0 on B/B fixture" 0 "$RC"
assert_contains "blood: plus-strand G>T rs8176746 1/1 -> B+ BB" "B+ (ABO genotype BB" "$OUT"
assert_not_contains "blood: no REF/ALT mismatch warning for G>T" "expected" "$ERR"

# Variants-only A/O genome: rs8176719 insertion 0/1 and no rs7853989 record
mk ao_vo.vcf '1|25362249|rs590787|G|A|100|PASS|.|GT|0/1' '9|133257521|rs8176719|T|TC|100|PASS|.|GT|0/1'
run blood_type.sh ao_vo.vcf
assert_eq "blood: rc 0 on variants-only A/O fixture" 0 "$RC"
assert_contains "blood: variants-only A/O -> A+ AO" "A+ (ABO genotype AO" "$OUT"
assert_contains "blood: variants-only A/O -> Confidence Moderate" "Confidence Level: Moderate" "$OUT"
assert_contains "blood: variants-only A/O -> reference assumed note" "rs7853989 record absent" "$OUT"
mk aa_vo.vcf '1|25362249|rs590787|G|A|100|PASS|.|GT|0/1' '9|133257521|rs8176719|T|TC|100|PASS|.|GT|1/1'
run blood_type.sh aa_vo.vcf
assert_contains "blood: variants-only A/A -> ABO genotype AA" "ABO genotype AA" "$OUT"

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]]
