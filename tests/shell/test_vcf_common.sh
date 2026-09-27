#!/usr/bin/env bash
# =============================================================================
# tests/shell/test_vcf_common.sh - plain bash test runner for
# scripts/lib/vcf_common.sh (no bats). Exit status is non-zero on any failure.
#
# Usage: bash tests/shell/test_vcf_common.sh
# Requires: bash >= 3.2, bcftools, bgzip, tabix, gzip.
# Fixtures are built under a temporary directory from data/sample_genome.vcf.gz
# and removed on exit; nothing is written into the repository.
# =============================================================================
set -uo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
LIB="$ROOT/scripts/lib/vcf_common.sh"
SAMPLE="$ROOT/data/sample_genome.vcf.gz"

for tool in bcftools bgzip tabix gzip; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "SKIP: $tool not installed" >&2
        exit 0
    fi
done
[[ -f $LIB ]] || { echo "FAIL: library not found: $LIB"; exit 1; }
[[ -f $SAMPLE ]] || { echo "FAIL: sample not found: $SAMPLE"; exit 1; }

TMP=$(mktemp -d "${TMPDIR:-/tmp}/hx_lib_test.XXXXXX")
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
    if [[ $3 == *"$2"* ]]; then pass "$1"; else fail "$1" "[$2] not found in [$3]"; fi
}
assert_empty() { # NAME ACTUAL
    if [[ -z $2 ]]; then pass "$1"; else fail "$1" "expected empty, got [$2]"; fi
}

# run_lib VCF CODE: run CODE in a fresh bash under `set -euo pipefail` with the
# library sourced and HX_VCF prepared from VCF. stdout -> $OUT, stderr -> $ERR,
# exit status -> $RC.
run_lib() {
    local vcf=$1 code=$2
    OUT=$(LIB="$LIB" VCF="$vcf" bash -euo pipefail -c \
        'source "$LIB"; hx_prepare_vcf "$VCF"; '"$code" 2>"$TMP/stderr.txt")
    RC=$?
    ERR=$(cat "$TMP/stderr.txt")
}
# run_raw CODE: same but without preparing a VCF.
run_raw() {
    local code=$1
    OUT=$(LIB="$LIB" bash -euo pipefail -c 'source "$LIB"; '"$code" 2>"$TMP/stderr.txt")
    RC=$?
    ERR=$(cat "$TMP/stderr.txt")
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
zcat "$SAMPLE" >plain.vcf
gzip -c plain.vcf >plaingz.vcf.gz
bgzip -c plain.vcf >noidx.vcf.gz
grep '^#' plain.vcf >headeronly.vcf
sed 's/^##contig=<ID=1,length=248956422>/##contig=<ID=1,length=249250621>/' plain.vcf >grch37.vcf
grep -v '^##contig=' plain.vcf >nocontig.vcf
bgzip -c nocontig.vcf >nocontig.vcf.gz
gzip -c nocontig.vcf >nocontiggz.vcf.gz
{ grep '^#' plain.vcf; grep -v '^#' plain.vcf | tac; } >unsorted.vcf
printf 'hello, not a vcf\n' >garbage.txt
# ##contig lines without length= (the build cannot be verified either way)
sed -E 's/^##contig=<ID=([^,>]+),length=[0-9]+>/##contig=<ID=\1>/' plain.vcf >nolen.vcf
# Records using FORMAT (DP, GQ) and INFO (FOO) tags the header does not declare:
# bcftools query only warns, annotate/sort refuse them (build path only).
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {$9="GT:DP:GQ"; $10=$10":30:99"; print}' plain.vcf >fmtextra.vcf
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {$8=$8";FOO=1"; print}' plain.vcf >undefinfo.vcf
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {$8=$8";DP=7"; $9="GT:DP"; $10=$10":7"; print}' plain.vcf >sametag.vcf

# chr-prefixed copy of the sample (bgzip + tabix, so the fast path must reject it)
{
    for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' "$c" "$c"; done
    printf 'MT\tchrM\n'
} >to_chr.txt
bcftools annotate --rename-chrs to_chr.txt -Oz -o chr.vcf.gz plain.vcf 2>/dev/null
tabix -f -p vcf chr.vcf.gz
# bcftools view (not annotate) so the undeclared tags survive into the copy
bcftools view fmtextra.vcf 2>/dev/null | sed -E '/^#/!s/^([0-9XYM]+)\t/chr\1\t/; s/^##contig=<ID=([0-9XY]+|MT)/##contig=<ID=chr\1/' | sed 's/^##contig=<ID=chrMT/##contig=<ID=chrM/; /^#/!s/^chrMT\t/chrM\t/' | bgzip -c >fmtextra_chr.vcf.gz
tabix -f -p vcf fmtextra_chr.vcf.gz

# Edge-case VCF: phased 1|1, missing ./., multi-allelic 1/2, half call ./1,
# ALT-swapped record, unphased 1/0, hom-ref record with ALT ".", an indel,
# mixed SNP+indel records (SNP het, indel het, missing), an MNP, a symbolic
# ALT and an MT record.
{
    grep '^##' plain.vcf | grep -v '^##contig=<ID=MT'
    printf '##contig=<ID=MT,length=16569>\n'
    grep '^#CHROM' plain.vcf
    printf '1\t1000\trsPHASED\tA\tG\t100\tPASS\tGENE=EDGE\tGT\t1|1\n'
    printf '1\t2000\trsMISSING\tC\tT\t100\tPASS\tGENE=EDGE\tGT\t./.\n'
    printf '1\t3000\trsMULTI\tG\tA,T\t100\tPASS\tGENE=EDGE\tGT\t1/2\n'
    printf '1\t4000\trsHALF\tT\tC\t100\tPASS\tGENE=EDGE\tGT\t./1\n'
    printf '1\t5000\trsSWAP\tG\tA\t100\tPASS\tGENE=EDGE\tGT\t0/1\n'
    printf '1\t6000\trsHOMREF\tA\t.\t100\tPASS\tGENE=EDGE\tGT\t0/0\n'
    printf '1\t7000\trsINDEL\tAT\tA\t100\tPASS\tGENE=EDGE\tGT\t0/1\n'
    printf '1\t8000\trsUNPHASED\tC\tG\t100\tPASS\tGENE=EDGE\tGT\t1/0\n'
    printf '1\t9000\trsMIXED\tC\tT,CTT\t100\tPASS\tGENE=EDGE\tGT\t0/1\n'
    printf '1\t9100\trsMIXEDINDEL\tC\tT,CTT\t100\tPASS\tGENE=EDGE\tGT\t0/2\n'
    printf '1\t9200\trsMIXEDMISS\tC\tT,CTT\t100\tPASS\tGENE=EDGE\tGT\t./.\n'
    printf '1\t9300\trsMIXEDFIRST\tC\tCTT,T\t100\tPASS\tGENE=EDGE\tGT\t0/1\n'
    printf '1\t9600\trsMIXEDFIRSTSNP\tC\tCTT,T\t100\tPASS\tGENE=EDGE\tGT\t0/2\n'
    printf '1\t9400\trsMNP\tCA\tTG\t100\tPASS\tGENE=EDGE\tGT\t0/1\n'
    printf '1\t9500\trsSYM\tC\t<DEL>\t100\tPASS\tGENE=EDGE\tGT\t0/1\n'
    printf 'MT\t100\trsMT\tA\tG\t100\tPASS\tGENE=EDGE\tGT\t1/1\n'
} >edge.vcf
bcftools annotate --rename-chrs to_chr.txt -o edge_chr.vcf edge.vcf 2>/dev/null

# ---------------------------------------------------------------------------
# Sourcing is silent and safe under set -euo pipefail
# ---------------------------------------------------------------------------
run_raw 'true'
assert_eq "source: rc 0 under set -euo pipefail" 0 "$RC"
assert_empty "source: prints nothing to stdout" "$OUT"
run_raw 'source "$LIB"; echo again'
assert_eq "source: safe to source twice" "again" "$OUT"

# ---------------------------------------------------------------------------
# hx_die / hx_require_arg
# ---------------------------------------------------------------------------
run_raw 'hx_die "boom"'
assert_eq "hx_die: default exit code 2" 2 "$RC"
assert_eq "hx_die: message on stderr" "ERROR: boom" "$ERR"
run_raw 'hx_die "boom" 5'
assert_eq "hx_die: custom exit code" 5 "$RC"
run_raw 'hx_require_arg "" "myscript <file.vcf.gz>"'
assert_eq "hx_require_arg: empty -> exit 2" 2 "$RC"
assert_contains "hx_require_arg: usage on stderr" "myscript <file.vcf.gz>" "$ERR"
run_raw 'hx_require_arg "x" "usage"; echo ok'
assert_eq "hx_require_arg: non-empty passes" "ok" "$OUT"

# ---------------------------------------------------------------------------
# hx_validate_vcf
# ---------------------------------------------------------------------------
run_raw "hx_validate_vcf '$SAMPLE' && echo ok"
assert_eq "hx_validate_vcf: BGZF sample" "ok" "$OUT"
run_raw "hx_validate_vcf '$TMP/plain.vcf' && echo ok"
assert_eq "hx_validate_vcf: plain vcf" "ok" "$OUT"
run_raw "hx_validate_vcf '$TMP/plaingz.vcf.gz' && echo ok"
assert_eq "hx_validate_vcf: plain gzip" "ok" "$OUT"
run_raw "hx_validate_vcf '$TMP/garbage.txt' || echo bad"
assert_eq "hx_validate_vcf: garbage -> 1" "bad" "$OUT"
run_raw "hx_validate_vcf '$TMP/does_not_exist.vcf' || echo bad"
assert_eq "hx_validate_vcf: missing -> 1" "bad" "$OUT"

# ---------------------------------------------------------------------------
# hx_prepare_vcf
# ---------------------------------------------------------------------------
run_lib "$SAMPLE" 'echo "$HX_VCF"'
assert_eq "prepare: BGZF+index+bare contigs used as is" "$SAMPLE" "$OUT"
if [[ -z $(ls "$HELIXIGHT_WORK_DIR"/*.vcf.gz 2>/dev/null) ]]; then
    pass "prepare: no working copy for the ready sample"
else
    fail "prepare: unexpected working copy for the ready sample"
fi

run_lib plain.vcf 'echo "$HX_VCF"; hx_gt 11:66560624'
assert_eq "prepare: plain .vcf -> genotype readable" "0/1" "${OUT##*$'\n'}"
PLAIN_WORK=${OUT%%$'\n'*}
assert_contains "prepare: plain .vcf working copy under HELIXIGHT_WORK_DIR" "$HELIXIGHT_WORK_DIR/" "$PLAIN_WORK"
[[ -s $PLAIN_WORK.tbi ]] && pass "prepare: working copy is tabix-indexed" || fail "prepare: missing $PLAIN_WORK.tbi"
assert_eq "prepare: working copy is BGZF" "1" "$(bcftools index -s "$PLAIN_WORK" >/dev/null 2>&1 && echo 1 || echo 0)"

run_lib plain.vcf 'echo "$HX_VCF"'
assert_eq "prepare: working copy reused on second call" "$PLAIN_WORK" "$OUT"
assert_eq "prepare: exactly one working copy per input" 1 "$(ls "$HELIXIGHT_WORK_DIR"/*.vcf.gz | wc -l)"

run_lib plaingz.vcf.gz 'hx_gt 11:66560624'
assert_eq "prepare: plain gzip -> genotype readable" "0/1" "$OUT"

run_lib noidx.vcf.gz 'hx_gt 11:66560624'
assert_eq "prepare: BGZF without index -> genotype readable" "0/1" "$OUT"
[[ ! -e noidx.vcf.gz.tbi && ! -e noidx.vcf.gz.csi ]] \
    && pass "prepare: never indexes the user's file in place" \
    || fail "prepare: created an index next to the user's file"

run_lib chr.vcf.gz 'echo "$HX_VCF"; hx_gt 11:66560624; bcftools view -h "$HX_VCF" | grep -c "^##contig=<ID=chr"'
assert_contains "prepare: chr-prefixed input -> working copy" "$HELIXIGHT_WORK_DIR/" "${OUT%%$'\n'*}"
assert_eq "prepare: chr-prefixed input -> bare region query works, no chr contigs left" $'0/1\n0' "${OUT#*$'\n'}"

run_lib edge_chr.vcf 'hx_gt 1:1000; hx_gt MT:100'
assert_eq "prepare: chr1/chrM renamed to 1/MT" $'1/1\n1/1' "$OUT"

run_lib unsorted.vcf 'bcftools view -H "$HX_VCF" | cut -f1,2'
assert_eq "prepare: unsorted input is sorted" "$(bcftools view -H "$SAMPLE" | cut -f1,2)" "$OUT"

# No ##contig lines and no index: bcftools 1.19 cannot sort/annotate such a
# file, so the working copy is reheadered with synthesised contig lines.
for f in nocontig.vcf nocontig.vcf.gz nocontiggz.vcf.gz; do
    run_lib "$f" 'hx_gt 11:66560624'
    assert_eq "prepare: unindexed $f without ##contig lines -> rc 0, genotype readable" "0|0/1" "$RC|$OUT"
    assert_contains "prepare: $f -> no-contig warning on stderr" "no ##contig lines" "$ERR"
done
run_lib nocontig.vcf 'hx_gt 11:66560624'
assert_contains "prepare: no-contig warning also on the cached second call" "no ##contig lines" "$ERR"

run_lib nolen.vcf 'hx_gt 11:66560624'
assert_eq "prepare: ##contig without length -> rc 0, genotype readable" "0|0/1" "$RC|$OUT"
assert_contains "prepare: ##contig without length -> cannot-verify warning" "contig 1 has no length" "$ERR"
run_lib plain.vcf 'true'
assert_empty "prepare: full ##contig lines -> no warning" "$ERR"

# Undeclared INFO/FORMAT tags: the build path declares them in the copy's header
for f in fmtextra.vcf undefinfo.vcf sametag.vcf fmtextra_chr.vcf.gz; do
    run_lib "$f" 'hx_gt 11:66560624; hx_report_coverage'
    assert_eq "prepare: $f (undeclared tags) -> rc 0, genotype readable" $'0|0/1\nHX_DATA_COVERAGE: 1/1 markers' "$RC|$OUT"
done
run_lib fmtextra.vcf 'bcftools view -h "$HX_VCF" | grep -c "^##FORMAT=<ID=\(DP\|GQ\),"'
assert_eq "prepare: undeclared FORMAT tags declared in the working copy" "2" "$OUT"
run_lib sametag.vcf 'bcftools view -h "$HX_VCF" | grep -c "^##\(INFO\|FORMAT\)=<ID=DP,"'
assert_eq "prepare: same tag undeclared as INFO and FORMAT -> both declared" "2" "$OUT"

# Default work dir (HELIXIGHT_WORK_DIR unset): per user, mode 700, private files
run_raw "env -u HELIXIGHT_WORK_DIR TMPDIR='$TMP/tmpdefault' bash -euo pipefail -c 'source \"\$LIB\"; hx_prepare_vcf plain.vcf; echo \"\$HX_VCF\"'"
DEF_WORK=$OUT
assert_eq "prepare: default work dir is per user" "$TMP/tmpdefault/helixight-$(id -u)/" "${DEF_WORK%/*}/"
if [[ -d $TMP/tmpdefault/helixight-$(id -u) ]]; then
    assert_eq "prepare: default work dir mode 700" "700" "$(stat -c '%a' "$TMP/tmpdefault/helixight-$(id -u)")"
    assert_eq "prepare: working copy not group/world readable" "600" "$(stat -c '%a' "$DEF_WORK")"
    assert_eq "prepare: index not group/world readable" "600" "$(stat -c '%a' "$DEF_WORK.tbi")"
else
    fail "prepare: default work dir not created"
fi
run_raw 'umask 022; source "$LIB"; hx_prepare_vcf plain.vcf; umask'
assert_eq "prepare: caller umask restored" "0022" "$OUT"

run_lib grch37.vcf 'echo unreachable'
assert_eq "prepare: GRCh37 header -> exit 3" 3 "$RC"
assert_contains "prepare: GRCh37 message" "GRCh37" "$ERR"
assert_empty "prepare: GRCh37 prints nothing to stdout" "$OUT"

run_lib headeronly.vcf 'hx_gt 11:66560624; hx_report_coverage'
assert_eq "prepare: header-only VCF -> rc 0" 0 "$RC"
assert_eq "prepare: header-only -> no genotype, coverage 0/1" "HX_DATA_COVERAGE: 0/1 markers" "$OUT"

run_lib garbage.txt 'echo unreachable'
assert_eq "prepare: garbage -> exit 2" 2 "$RC"
assert_contains "prepare: garbage -> ERROR on stderr" "ERROR:" "$ERR"

run_lib does_not_exist.vcf 'echo unreachable'
assert_eq "prepare: missing file -> exit 2" 2 "$RC"

run_raw 'hx_prepare_vcf ""'
assert_eq "prepare: empty path -> exit 2" 2 "$RC"

# ---------------------------------------------------------------------------
# hx_gt
# ---------------------------------------------------------------------------
run_lib edge.vcf 'hx_gt 1:1000'
assert_eq "hx_gt: phased 1|1 -> 1/1" "1/1" "$OUT"
run_lib edge.vcf 'hx_gt 1:2000'
assert_empty "hx_gt: missing ./. -> empty" "$OUT"
run_lib edge.vcf 'hx_gt 1:3000'
assert_eq "hx_gt: multi-allelic 1/2 kept" "1/2" "$OUT"
run_lib edge.vcf 'hx_gt 1:4000'
assert_empty "hx_gt: half call ./1 -> empty" "$OUT"
run_lib edge.vcf 'hx_gt 1:8000'
assert_eq "hx_gt: 1/0 -> 0/1" "0/1" "$OUT"
run_lib edge.vcf 'hx_gt 1:6000'
assert_eq "hx_gt: hom-ref record with ALT . -> 0/0" "0/0" "$OUT"
run_lib edge.vcf 'hx_gt 1:7000'
assert_empty "hx_gt: indel record -> empty" "$OUT"
run_lib edge.vcf 'hx_gt 1:9000'
assert_eq "hx_gt: SNP+indel record, SNP het -> 0/1" "0/1" "$OUT"
run_lib edge.vcf 'hx_gt 1:9100'
assert_empty "hx_gt: SNP+indel record, indel het (0/2) -> empty" "$OUT"
run_lib edge.vcf 'hx_gt 1:9300'
assert_empty "hx_gt: indel+SNP record, indel het (0/1) -> empty" "$OUT"
run_lib edge.vcf 'hx_gt 1:9600'
assert_eq "hx_gt: indel+SNP record, SNP het (0/2) -> 0/2 (record's own ALT index)" "0/2" "$OUT"
run_lib edge.vcf 'hx_gt 1:9200'
assert_empty "hx_gt: SNP+indel record, missing call -> empty" "$OUT"
run_lib edge.vcf 'hx_gt 1:9400'
assert_empty "hx_gt: MNP record -> empty" "$OUT"
run_lib edge.vcf 'hx_gt 1:9500'
assert_empty "hx_gt: symbolic ALT -> empty" "$OUT"
run_lib edge.vcf 'a=$(hx_gt 1:9000); b=$(hx_gt 1:9100); hx_report_coverage'
assert_eq "hx_gt: mixed records tallied (SNP het found, indel het not)" "HX_DATA_COVERAGE: 1/2 markers" "$OUT"
run_lib edge.vcf 'hx_gt 1:999999'
assert_eq "hx_gt: no record -> empty, rc 0" "0|" "$RC|$OUT"
run_lib "$SAMPLE" 'hx_gt 15:74749576; hx_gt 12:111803962'
assert_eq "hx_gt: sample 1/1 and 0/0" $'1/1\n0/0' "$OUT"

# ---------------------------------------------------------------------------
# hx_gt_indel
# ---------------------------------------------------------------------------
run_lib "$SAMPLE" 'hx_gt_indel 17:63488529 250'
read -r pos ref alt gt <<<"$OUT"
assert_eq "hx_gt_indel: ACE record found" "63488529 288 T 0/1" "$pos ${#ref} $alt $gt"
run_lib "$SAMPLE" 'hx_gt_indel 17:63488529 300'
assert_empty "hx_gt_indel: MINLEN larger than the indel -> empty" "$OUT"
run_lib "$SAMPLE" 'hx_gt_indel 11:66560624 1'
assert_empty "hx_gt_indel: SNP record -> empty" "$OUT"
run_lib "$SAMPLE" 'hx_gt_indel 9:133257520-133257521 1'
assert_eq "hx_gt_indel: ABO 1-bp insertion (non-O allele; reference = 261delG)" "133257521 T TC 0/1" "$OUT"
run_lib edge.vcf 'hx_gt_indel 1:7000 1'
assert_eq "hx_gt_indel: edge indel" "7000 AT A 0/1" "$OUT"
run_lib edge.vcf 'hx_gt_indel 1:2000 1'
assert_empty "hx_gt_indel: missing call -> empty" "$OUT"

# ---------------------------------------------------------------------------
# hx_region_present
# ---------------------------------------------------------------------------
run_lib edge.vcf 'hx_region_present 1:5500-6500; hx_region_present 1:1500-2500; hx_region_present 1:999000-999999; hx_report_coverage'
assert_eq "hx_region_present: hom-ref ALT . counts, missing call and no record do not" "HX_DATA_COVERAGE: 1/3 markers" "$OUT"
run_lib edge.vcf 'hx_region_present 1:6500-7500'
assert_eq "hx_region_present: prints nothing, rc 0" "0|" "$RC|$OUT"

# ---------------------------------------------------------------------------
# hx_gt_allele
# ---------------------------------------------------------------------------
run_lib "$SAMPLE" 'hx_gt_allele 11:66560624 C T rs1815739'
assert_eq "hx_gt_allele: het with matching rsID" $'1\t0/1' "$OUT"
assert_empty "hx_gt_allele: no warning on match" "$ERR"
run_lib "$SAMPLE" 'hx_gt_allele 15:74749576-74749576 A C'
assert_eq "hx_gt_allele: hom alt -> dosage 2" $'2\t1/1' "$OUT"
run_lib "$SAMPLE" 'hx_gt_allele 12:111803962 G A'
assert_eq "hx_gt_allele: hom ref -> dosage 0" $'0\t0/0' "$OUT"
run_lib "$SAMPLE" 'hx_gt_allele 11:66560624 C T rs9999999'
assert_empty "hx_gt_allele: rsID mismatch -> empty" "$OUT"
assert_contains "hx_gt_allele: rsID mismatch -> warning" "expected C>T - skipped" "$ERR"
run_lib "$SAMPLE" 'hx_gt_allele 11:66560624 T C'
assert_empty "hx_gt_allele: swapped REF/ALT -> empty" "$OUT"
assert_contains "hx_gt_allele: swapped REF/ALT -> warning names the record" "has REF/ALT C/T, expected T>C" "$ERR"
run_lib edge.vcf 'hx_gt_allele 1:5000 A G'
assert_eq "hx_gt_allele: ALT-swapped edge record -> empty, rc 0" "0|" "$RC|$OUT"
assert_contains "hx_gt_allele: ALT-swapped edge record -> warning" "warning: 1:5000 has REF/ALT G/A, expected A>G - skipped" "$ERR"
run_lib edge.vcf 'hx_gt_allele 1:3000 G T'
assert_eq "hx_gt_allele: multi-allelic, second ALT" $'1\t1/2' "$OUT"
run_lib edge.vcf 'hx_gt_allele 1:3000 G A'
assert_eq "hx_gt_allele: multi-allelic, first ALT" $'1\t1/2' "$OUT"
run_lib edge.vcf 'hx_gt_allele 1:9000 C T'
assert_eq "hx_gt_allele: SNP+indel record, SNP allele" $'1\t0/1' "$OUT"
run_lib edge.vcf 'hx_gt_allele 1:9000 c t'
assert_eq "hx_gt_allele: lower-case REF/ALT arguments" $'1\t0/1' "$OUT"
run_lib edge.vcf 'hx_gt_allele 1:3000 G C'
assert_empty "hx_gt_allele: multi-allelic, absent ALT -> empty" "$OUT"
assert_contains "hx_gt_allele: multi-allelic, absent ALT -> warning" "expected G>C" "$ERR"
run_lib edge.vcf 'hx_gt_allele 1:6000 A G'
assert_eq "hx_gt_allele: hom-ref record with ALT . -> dosage 0" $'0\t0/0' "$OUT"
run_lib edge.vcf 'hx_gt_allele 1:2000 C T'
assert_eq "hx_gt_allele: missing ./. -> empty, no warning" "|" "$OUT|$ERR"
run_lib edge.vcf 'hx_gt_allele 1:4000 T C'
assert_eq "hx_gt_allele: half call ./1 -> empty, no warning" "|" "$OUT|$ERR"
run_lib edge.vcf 'hx_gt_allele 1:1000 A G'
assert_eq "hx_gt_allele: phased 1|1 -> 2 and normalised GT" $'2\t1/1' "$OUT"
run_lib edge.vcf 'hx_gt_allele 1:8000 C G'
assert_eq "hx_gt_allele: 1/0 -> dosage 1, GT 0/1" $'1\t0/1' "$OUT"
run_lib edge.vcf 'hx_gt_allele 1:999999 A G'
assert_eq "hx_gt_allele: no record -> empty, no warning, rc 0" "0||" "$RC|$OUT|$ERR"
run_lib "$SAMPLE" 'hx_gt_allele 9:133257521 G A rs8176719'
assert_empty "hx_gt_allele: overlapping deletion at a SNP position is not scored" "$OUT"
assert_contains "hx_gt_allele: overlapping deletion -> warning" "has REF/ALT T/TC, expected G>A" "$ERR"
run_lib edge.vcf 'hx_gt_allele 1:1000'
assert_eq "hx_gt_allele: missing arguments -> exit 2" 2 "$RC"

# ---------------------------------------------------------------------------
# hx_alt_dosage
# ---------------------------------------------------------------------------
run_raw 'for g in 0/0 0/1 1/1 "1|1" "0|1" 1/2 "./1" "./." "" 1 "1/0"; do hx_alt_dosage "$g"; done'
assert_eq "hx_alt_dosage" $'0\n1\n2\n2\n1\n2\nNA\nNA\nNA\n1\n1' "$OUT"

# ---------------------------------------------------------------------------
# hx_count
# ---------------------------------------------------------------------------
printf 'BRCA1 x\nBRCA1 y\nother\n' >count.txt
run_raw "hx_count '$TMP/count.txt' BRCA1"
assert_eq "hx_count: matches" "2" "$OUT"
run_raw "hx_count '$TMP/count.txt' ZZZ"
assert_eq "hx_count: no match prints exactly one 0" "0" "$OUT"
run_raw "hx_count '$TMP/no_such_file' BRCA1"
assert_eq "hx_count: missing file -> 0, rc 0" "0|0" "$RC|$OUT"
run_raw "hx_count '$TMP/count.txt' '^-x'"
assert_eq "hx_count: pattern starting with a dash" "0" "$OUT"

# ---------------------------------------------------------------------------
# hx_query
# ---------------------------------------------------------------------------
run_lib "$SAMPLE" "hx_query 11:66560624 '%REF %ALT [%GT]\n'"
assert_eq "hx_query: format and region" "C T 0/1" "$OUT"
run_lib "$SAMPLE" "hx_query 11:66560624 '%POS\n' -e 'GT=\"het\"'; echo end"
assert_eq "hx_query: extra bcftools args are passed through" "end" "$OUT"
run_lib "$SAMPLE" "hx_query 11:66560624 '%NO_SUCH_TAG\n'; echo unreachable"
[[ $RC -ne 0 && $OUT != *unreachable* ]] \
    && pass "hx_query: bcftools error aborts a set -e caller" \
    || fail "hx_query: bcftools error did not abort the caller"
run_raw "hx_query 11:66560624 '%POS\n'"
assert_eq "hx_query: without hx_prepare_vcf -> exit 2" 2 "$RC"

# ---------------------------------------------------------------------------
# hx_report_coverage (counters must survive $(...) subshells)
# ---------------------------------------------------------------------------
run_lib "$SAMPLE" 'a=$(hx_gt 11:66560624); b=$(hx_gt 1:999999); c=$(hx_gt_indel 17:63488529 250); d=$(hx_gt_allele 12:111803962 G A); e=$(hx_gt_allele 11:66560624 T C); hx_report_coverage; echo "$HX_FOUND/$HX_QUERIED"'
assert_eq "hx_report_coverage: counts found/queried across subshells" $'HX_DATA_COVERAGE: 3/5 markers\n3/5' "$OUT"
run_lib "$SAMPLE" 'hx_report_coverage'
assert_eq "hx_report_coverage: no queries -> 0/0" "HX_DATA_COVERAGE: 0/0 markers" "$OUT"
run_lib "$SAMPLE" 'x=$(hx_gt 11:66560624); [[ -s $HX_COUNTER_FILE ]] && echo present; hx_report_coverage >/dev/null; [[ -e $HX_COUNTER_FILE ]] || echo gone'
assert_eq "hx_report_coverage: tally file written by queries and removed by the report" $'present\ngone' "$OUT"
run_raw '[[ -e $HX_COUNTER_FILE ]] || echo none'
assert_eq "source: no tally file created before the first query" "none" "$OUT"

# ---------------------------------------------------------------------------
# Sample fixture sanity (what the scripts rely on)
# ---------------------------------------------------------------------------
n_records=$(bcftools view -H "$SAMPLE" | wc -l)
n_ids=$(bcftools query -f '%ID\n' "$SAMPLE" | sort -u | wc -l)
n_pos=$(bcftools query -f '%CHROM:%POS\n' "$SAMPLE" | sort -u | wc -l)
assert_eq "sample: unique IDs and CHROM:POS" "$n_records/$n_records" "$n_ids/$n_pos"
[[ $n_records -ge 50 ]] && pass "sample: has >= 50 records ($n_records)" || fail "sample: only $n_records records"
gts=$(bcftools query -f '[%GT]\n' "$SAMPLE" | sort -u | tr '\n' ' ')
assert_eq "sample: mix of 0/0, 0/1 and 1/1" "0/0 0/1 1/1 " "$gts"
if [[ -f $ROOT/data/sample_genome.vcf ]]; then
    cmp -s <(zcat "$SAMPLE") "$ROOT/data/sample_genome.vcf" \
        && pass "sample: plain .vcf identical to the .vcf.gz" \
        || fail "sample: plain .vcf differs from the .vcf.gz"
fi

echo
echo "Results: $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]]
