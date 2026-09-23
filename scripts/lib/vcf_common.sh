#!/usr/bin/env bash
# shellcheck shell=bash
# =============================================================================
# scripts/lib/vcf_common.sh - shared VCF helpers for the Helixight scripts
# =============================================================================
#
# Sourcing idiom (from any script in scripts/):
#   source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"
#
# The library is safe to source under `set -euo pipefail`: it defines only
# functions and a few globals, never exits on load, never uses `head` or
# other early-exit consumers in pipelines (SIGPIPE under pipefail), and never
# prints anything to stdout except the documented outputs below.
# Requires: bash >= 3.2 (no bash-4-only expansions such as ${var^^}, so the
# stock macOS /bin/bash works; some scripts still need bash >= 4 for
# `declare -A`), bcftools (>= 1.10). Optional: md5sum (or md5) for the
# work-copy cache key.
#
# API
# ---
# hx_die MSG [CODE]
#     Print "ERROR: MSG" to stderr and exit CODE (default 2).
# hx_require_arg VALUE USAGE
#     hx_die with "USAGE" when VALUE is empty (replaces `${1:-default.vcf.gz}`).
# hx_warn MSG
#     Print "warning: MSG" to stderr. Never exits.
# hx_validate_vcf PATH
#     `bcftools view -h PATH` succeeds (plain, gzip or BGZF). Returns 1 on
#     failure. Replaces `zcat | head -1 | grep -q` checks (which SIGPIPE under
#     pipefail for any file larger than ~64 KB decompressed).
# hx_prepare_vcf PATH
#     Sets and exports HX_VCF: a bgzip-compressed, indexed, sorted VCF with
#     bare contig names (chr1 -> 1, chrM/chrMT/M -> MT). When PATH is already
#     BGZF + indexed with bare contigs it is used as is; otherwise a working
#     copy is built under ${HELIXIGHT_WORK_DIR:-${TMPDIR:-/tmp}/helixight-$UID}/
#     <md5 of "path|size|mtime">.vcf.gz and reused on later calls. The copy is
#     a full copy of the genome: the default directory is per user and created
#     mode 700, and every file of the build is written under umask 077. Nothing
#     removes the copy (delete the directory to purge it). Exits 2 on
#     unreadable/invalid input, 3 when the file looks like GRCh37 (contig 1
#     length 249250621); warns on stderr when the build cannot be verified:
#     no ##contig lines in the header (bcftools 1.19 cannot sort or annotate
#     such a file unless it is tabix-indexed, so the working copy is first
#     reheadered with one ##contig=<ID=x> line per contig seen in the records)
#     or a contig 1 line without length=. Records that use INFO/FORMAT tags
#     the header does not declare (bcftools query only warns, annotate/sort
#     refuse them) are handled by declaring the tags as Type=String in the
#     working copy's header. Never indexes or modifies the caller's file in
#     place.
# hx_query REGION FORMAT [EXTRA bcftools query args...]
#     `bcftools query -r REGION -f FORMAT EXTRA... "$HX_VCF"`. stderr is not
#     suppressed; returns bcftools' status. Empty output means "no record".
# hx_gt REGION
#     Prints the genotype of the first SNP record at REGION with the phase
#     normalised ("|" -> "/", "1/0" -> "0/1"); missing and half calls are
#     excluded, and so are records with no SNP allele at all (pure indel/MNP/
#     symbolic: -e 'GT="mis" || (TYPE!~"snp" && TYPE!~"ref")'). A mixed
#     record (e.g. C>T,CTT, TYPE=SNP,INDEL) is kept and prints its GT with the
#     record's own ALT indexes, unless a called allele is the non-SNP one
#     (0/2 for C>T,CTT), in which case nothing is printed so callers never
#     score an indel carrier as a SNP genotype. Explicit hom-ref records (ALT
#     ".", as written by `bcftools convert --tsv2vcf` and all-sites callers)
#     are kept and print 0/0. Prints nothing when there is no usable record.
# hx_gt_indel REGION MINLEN
#     Like hx_gt but prints "POS REF ALT GT" of the first record (any type)
#     whose |len(REF) - len(ALT)| >= MINLEN (e.g. the ACE Alu I/D at
#     17:63488529, MINLEN 250). Prints nothing when there is none.
# hx_gt_allele REGION REF ALT [RSID]
#     Fetches every record at REGION and picks the one whose POS equals the
#     region start, whose REF equals REF and whose comma-split ALT list
#     contains ALT (a hom-ref record with ALT "." and matching REF also
#     matches, dosage 0). If RSID is given, a record whose ID is neither "."
#     nor RSID does not match. Prints "DOSAGE<TAB>GT" (dosage 0/1/2 = number of
#     GT alleles equal to the matching ALT index; GT normalised as in hx_gt).
#     Prints nothing for "no record" or for a missing/half call. When records
#     exist but none matches, prints nothing and writes
#     "warning: REGION has REF/ALT <x>/<y>, expected REF>ALT - skipped" to
#     stderr (this is what stops a wrong variant at the same position from
#     being scored silently).
# hx_region_present REGION
#     Tallies REGION as one coverage marker, found when it holds at least one
#     called record of any type (hom-ref ALT "." included). Prints nothing.
# hx_alt_dosage GT
#     Prints 0/1/2 (number of non-reference alleles) for 0/0, 0/1, 1/1, 1|0,
#     1/2 ...; prints NA for missing/half calls or an empty GT.
# hx_count FILE PATTERN
#     Prints the number of lines of FILE matching PATTERN (grep -c), always
#     exactly one number (grep -c prints 0 AND exits 1 on no match, so
#     `grep -c ... || echo 0` produced "0\n0"). Missing file -> 0.
# hx_report_coverage
#     Prints "HX_DATA_COVERAGE: <found>/<queried> markers" from the counters
#     that hx_gt, hx_gt_indel and hx_gt_allele maintain (HX_QUERIED is
#     incremented on every call, HX_FOUND whenever a genotype is printed).
#     Because callers run these functions inside `$(...)` subshells, the tally
#     is also appended to a per-process file (HX_COUNTER_FILE, under the work
#     dir) which hx_report_coverage reads, removes, and copies back into
#     HX_QUERIED/HX_FOUND. Scripts call it once at the end; the frontend
#     parses it and treats found == 0 as "no usable data" instead of a
#     successful run.
#
# Globals: HX_VCF (exported by hx_prepare_vcf), HX_QUERIED, HX_FOUND,
#          HX_COUNTER_FILE.
# =============================================================================

# Guard against double sourcing (counters must not be reset).
if [[ -n ${HX_VCF_COMMON_LOADED:-} ]]; then
    return 0 2>/dev/null || true
fi
HX_VCF_COMMON_LOADED=1

HX_VCF="${HX_VCF:-}"
HX_QUERIED=0
HX_FOUND=0
# Per-process tally file ($$ is the top-level script's PID even in subshells).
# The file is created lazily on the first tally and removed by
# hx_report_coverage; a stale file left by an earlier process with the same
# PID is removed here so counts never carry over.
HX_COUNTER_FILE="${HX_COUNTER_FILE:-${HELIXIGHT_WORK_DIR:-${TMPDIR:-/tmp}/helixight-${UID:-$(id -u)}}/hx_coverage.$$.txt}"
rm -f "$HX_COUNTER_FILE" 2>/dev/null || true

# Internal: record one query (Q) or one found genotype (F).
_hx_tally() {
    if [[ $1 == F ]]; then
        HX_FOUND=$((HX_FOUND + 1))
    else
        HX_QUERIED=$((HX_QUERIED + 1))
    fi
    {
        mkdir -p "$(dirname "$HX_COUNTER_FILE")" \
            && printf '%s\n' "$1" >>"$HX_COUNTER_FILE"
    } 2>/dev/null || true
}

hx_die() {
    printf 'ERROR: %s\n' "${1:-unknown error}" >&2
    exit "${2:-2}"
}

hx_warn() {
    printf 'warning: %s\n' "$*" >&2
}

hx_require_arg() {
    if [[ -z ${1:-} ]]; then
        hx_die "missing argument. Usage: ${2:-<file.vcf|file.vcf.gz>}"
    fi
}

hx_validate_vcf() {
    [[ -n ${1:-} && -f $1 && -r $1 ]] || return 1
    bcftools view -h "$1" >/dev/null 2>&1
}

# Internal: normalise a genotype string ("1|0" -> "0/1", "2/1" -> "1/2").
_hx_norm_gt() {
    local gt=${1//|//} a b
    if [[ $gt == */* ]]; then
        a=${gt%%/*}
        b=${gt#*/}
        if [[ $a =~ ^[0-9]+$ && $b =~ ^[0-9]+$ && $a -gt $b ]]; then
            gt="$b/$a"
        fi
    fi
    printf '%s' "$gt"
}

# Internal: 0 when GT is a complete call (no missing allele), 1 otherwise.
_hx_gt_is_called() {
    [[ -n $1 && $1 != *.* ]]
}

# Internal: upper-case a string (${var^^} is bash 4 only; macOS ships 3.2).
_hx_upper() {
    printf '%s' "$1" | tr '[:lower:]' '[:upper:]'
}

# Internal: INFO/FORMAT tags used by the records of a VCF but not declared in
# its header, one "INFO TAG" / "FORMAT TAG" line each (from the bcftools
# parser warnings; full pass over the file).
_hx_undeclared_tags() {
    bcftools view -H "$1" 2>&1 >/dev/null \
        | grep -oE "^\[W::vcf_parse[A-Za-z0-9_]*\] (INFO|FORMAT) '[^']+'" \
        | sed -E "s/^\[[^]]*\] //; s/'//g" | sort -u
}

# Internal: md5 of a string, with a fallback for systems without md5sum.
_hx_md5() {
    local sum
    if command -v md5sum >/dev/null 2>&1; then
        sum=$(printf '%s' "$1" | md5sum)
        sum=${sum%% *}
    elif command -v md5 >/dev/null 2>&1; then
        sum=$(printf '%s' "$1" | md5 -q)
    else
        sum=$(printf '%s' "$1" | cksum)
        sum=${sum%% *}
    fi
    printf '%s' "$sum"
}

# Internal: "size|mtime" of a file (GNU stat, then BSD stat, then 0|0).
_hx_stat_key() {
    local s
    if s=$(stat -c '%s|%Y' -- "$1" 2>/dev/null); then
        printf '%s' "$s"
    elif s=$(stat -f '%z|%m' -- "$1" 2>/dev/null); then
        printf '%s' "$s"
    else
        printf '0|0'
    fi
}

# Internal: write the chr -> bare contig rename map once per work dir.
_hx_rename_map() {
    local dir=$1 map=$1/rename_chrs.txt tmp c
    if [[ ! -s $map ]]; then
        tmp="$map.$$.tmp"
        {
            for c in $(seq 1 22) X Y; do
                printf 'chr%s\t%s\n' "$c" "$c"
            done
            printf 'chrM\tMT\nchrMT\tMT\nM\tMT\n'
        } >"$tmp" && mv -f "$tmp" "$map"
    fi
    printf '%s' "$map"
}

hx_prepare_vcf() {
    local src=${1:-} abs dir base work_dir key work tmp log map idx first_contig hdr re len
    local no_contig="" fixhdr taghdr tags kind tag old_umask

    [[ -n $src ]] || hx_die "hx_prepare_vcf: no VCF path given"
    [[ -f $src && -r $src ]] || hx_die "VCF file not found or not readable: $src"
    hx_validate_vcf "$src" || hx_die "not a readable VCF: $src (input must be a valid VCF; run: bcftools sort -Oz -o x.vcf.gz in.vcf && bcftools index -t x.vcf.gz)"

    dir=$(cd "$(dirname "$src")" && pwd) || hx_die "cannot resolve directory of $src"
    base=$(basename "$src")
    abs="$dir/$base"

    # Checked on the SOURCE header: a working copy built below gets synthesised
    # ##contig lines (without length), so the build guard could not tell.
    hdr=$(bcftools view -h "$abs" 2>/dev/null) || hx_die "cannot read header of $src"
    if [[ $hdr != *'##contig=<'* ]]; then
        no_contig=1
        hx_warn "no ##contig lines in $src - cannot verify the genome build (Helixight assumes GRCh38)"
    fi

    # Fast path: BGZF + index present and bare contig names -> use as is.
    idx=""
    if idx=$(bcftools index -s "$abs" 2>/dev/null); then
        first_contig=${idx%%$'\n'*}
        first_contig=${first_contig%%$'\t'*}
        if [[ $first_contig != chr* && $first_contig != M ]]; then
            HX_VCF=$abs
        fi
    fi

    if [[ -z ${HX_VCF:-} || $HX_VCF != "$abs" ]]; then
        # The working copy is a full copy of the genome: the default location
        # is per user and private, and everything built below is written
        # under umask 077 (restored before returning).
        old_umask=$(umask)
        umask 077
        if [[ -n ${HELIXIGHT_WORK_DIR:-} ]]; then
            work_dir=$HELIXIGHT_WORK_DIR
            mkdir -p "$work_dir" || hx_die "cannot create work directory: $work_dir"
        else
            work_dir="${TMPDIR:-/tmp}/helixight-${UID:-$(id -u)}"
            { mkdir -p "$work_dir" && chmod 700 "$work_dir"; } \
                || hx_die "cannot create work directory: $work_dir"
        fi
        key=$(_hx_md5 "$abs|$(_hx_stat_key "$abs")")
        work="$work_dir/$key.vcf.gz"
        if [[ -s $work && ( -s $work.tbi || -s $work.csi ) ]]; then
            HX_VCF=$work
        else
            map=$(_hx_rename_map "$work_dir")
            tmp="$work.$$.tmp"
            log="$work_dir/$key.log"
            : >"$log"
            fixhdr=""
            if [[ -n $no_contig ]]; then
                # Synthesise ##contig lines from the contigs seen in the records
                # (full scan of the input; only done for files without them).
                fixhdr="$work_dir/$key.hdr"
                if ! {
                    printf '%s\n' "$hdr" | grep -v '^#CHROM'
                    bcftools view -H "$abs" 2>>"$log" | cut -f1 | awk '!seen[$1]++ {printf "##contig=<ID=%s>\n", $1}'
                    printf '%s\n' "$hdr" | grep '^#CHROM'
                } >"$fixhdr" 2>>"$log"; then
                    rm -f "$fixhdr"
                    hx_die "could not build a contig header for $src: $(tr '\n' ' ' <"$log" 2>/dev/null)"
                fi
            fi
            taghdr=""
            while :; do
                if (
                    set -o pipefail
                    if [[ -n $fixhdr ]]; then
                        # reheader cannot read gzip input itself (-dcf passes plain
                        # text through), and BGZF is gzip.
                        gzip -dcf "$abs" | bcftools reheader -h "$fixhdr" - 2>>"$log"
                    else
                        cat "$abs"
                    fi \
                        | bcftools annotate --rename-chrs "$map" - 2>>"$log" \
                        | bcftools sort -T "$work_dir/sort.$$.XXXXXX" -Oz -o "$tmp" - 2>>"$log"
                ) && bcftools index -f -t "$tmp" 2>>"$log"; then
                    break
                fi
                rm -f "$tmp" "$tmp.tbi"
                if [[ -z $taghdr ]] && grep -q 'is not defined in the header' "$log" 2>/dev/null; then
                    # annotate/sort refuse records that use INFO/FORMAT tags the
                    # header does not declare (query only warns): declare them
                    # as Type=String in the working copy's header and retry once.
                    tags=$(_hx_undeclared_tags "$abs")
                    if [[ -n $tags ]]; then
                        taghdr="$work_dir/$key.tags.hdr"
                        if {
                            if [[ -n $fixhdr ]]; then
                                grep -v '^#CHROM' "$fixhdr"
                            else
                                printf '%s\n' "$hdr" | grep -v '^#CHROM'
                            fi
                            while read -r kind tag; do
                                [[ -n $kind && -n $tag ]] || continue
                                printf '##%s=<ID=%s,Number=.,Type=String,Description="not declared in the input header">\n' "$kind" "$tag"
                            done <<<"$tags"
                            printf '%s\n' "$hdr" | grep '^#CHROM'
                        } >"$taghdr" 2>>"$log"; then
                            [[ -n $fixhdr ]] && rm -f "$fixhdr"
                            fixhdr=$taghdr
                            : >"$log"
                            continue
                        fi
                        rm -f "$taghdr"
                    fi
                fi
                rm -f "$fixhdr"
                hx_die "could not prepare VCF $src: $(grep -v -e '^Writing to' -e '^Merging' -e '^Cleaning' -e '^Done' "$log" 2>/dev/null | tr '\n' ' ')
input must be a valid VCF; run: bcftools sort -Oz -o x.vcf.gz in.vcf && bcftools index -t x.vcf.gz$(grep -q 'is not defined in the header' "$log" 2>/dev/null && printf '\n(the records use INFO/FORMAT tags the header does not declare: add ##INFO/##FORMAT lines for them with bcftools reheader -h or bcftools annotate -h before sorting)')"
            done
            mv -f "$tmp.tbi" "$work.tbi" && mv -f "$tmp" "$work" \
                || hx_die "could not move prepared VCF into place: $work"
            [[ -n $fixhdr ]] && rm -f "$fixhdr"
            HX_VCF=$work
        fi
        umask "$old_umask"
    fi

    # Build guard: GRCh37 chr1 is 249250621 bp, GRCh38 chr1 is 248956422 bp.
    # (Read from the prepared copy, whose contig names are bare.)
    if [[ -z $no_contig ]]; then
        hdr=$(bcftools view -h "$HX_VCF" 2>/dev/null) || hx_die "cannot read header of $HX_VCF"
        re="##contig=<ID=1,[^>"$'\n'"]*length=([0-9]+)"
        if [[ $hdr =~ $re ]]; then
            len=${BASH_REMATCH[1]}
            if [[ $len == 249250621 ]]; then
                hx_die "VCF appears to be GRCh37 (chr1 length 249250621); Helixight requires GRCh38" 3
            fi
        elif [[ $hdr == *'##contig=<ID=1'[,\>]* ]]; then
            hx_warn "contig 1 has no length in the header of $src - cannot verify the genome build (Helixight assumes GRCh38)"
        fi
    fi

    export HX_VCF
}

hx_query() {
    local region=${1:-} fmt=${2:-}
    [[ -n ${HX_VCF:-} ]] || hx_die "hx_query: HX_VCF is not set (call hx_prepare_vcf first)"
    [[ -n $region && -n $fmt ]] || hx_die "hx_query: usage: hx_query REGION FORMAT [bcftools query args]"
    bcftools query -r "$region" -f "$fmt" "${@:3}" "$HX_VCF"
}

hx_gt() {
    local out rc=0 ref alt gt a
    local -a alts called
    _hx_tally Q
    # TYPE!~"snp" (not TYPE!="snp"): a record with a SNP allele next to an
    # indel allele (C>T,CTT, TYPE=SNP,INDEL) still carries the SNP call.
    out=$(hx_query "${1:-}" '%REF\t%ALT\t[%GT\t]\n' -e 'GT="mis" || (TYPE!~"snp" && TYPE!~"ref")') || rc=$?
    [[ $rc -eq 0 ]] || return "$rc"
    out=${out%%$'\n'*}
    IFS=$'\t' read -r ref alt gt _ <<<"$out"
    _hx_gt_is_called "$gt" || return 0
    gt=$(_hx_norm_gt "$gt")
    if [[ $alt != '.' ]]; then
        # Mixed record: a called allele that is not a SNP allele (length
        # differs from REF) must not be reported as a SNP genotype.
        IFS=',' read -r -a alts <<<"$alt"
        IFS='/' read -r -a called <<<"$gt"
        for a in "${called[@]}"; do
            [[ $a =~ ^[0-9]+$ && $a -gt 0 ]] || continue
            [[ ${#alts[a-1]} -eq ${#ref} ]] || return 0
        done
    fi
    _hx_tally F
    printf '%s\n' "$gt"
}

hx_region_present() {
    local out rc=0
    _hx_tally Q
    out=$(hx_query "${1:-}" '%POS\n' -e 'GT="mis"') || rc=$?
    [[ $rc -eq 0 ]] || return "$rc"
    [[ -n $out ]] && _hx_tally F
    return 0
}

hx_gt_indel() {
    local region=${1:-} minlen=${2:-1} out pos ref alt gt a diff rc=0
    _hx_tally Q
    out=$(hx_query "$region" '%POS\t%REF\t%ALT\t[%GT\t]\n' -e 'GT="mis"') || rc=$?
    [[ $rc -eq 0 ]] || return "$rc"
    [[ -n $out ]] || return 0
    while IFS=$'\t' read -r pos ref alt gt _; do
        [[ -n $pos ]] || continue
        _hx_gt_is_called "$gt" || continue
        local IFS_save=$IFS
        IFS=',' read -r -a alts <<<"$alt"
        IFS=$IFS_save
        for a in "${alts[@]}"; do
            [[ $a == '.' || $a == '<'* || $a == '*' ]] && continue
            diff=$(( ${#ref} - ${#a} ))
            [[ $diff -lt 0 ]] && diff=$(( -diff ))
            if [[ $diff -ge $minlen ]]; then
                _hx_tally F
                printf '%s %s %s %s\n' "$pos" "$ref" "$alt" "$(_hx_norm_gt "$gt")"
                return 0
            fi
        done
    done <<<"$out"
    return 0
}

hx_gt_allele() {
    local region=${1:-} want_ref=${2:-} want_alt=${3:-} rsid=${4:-}
    local start out pos id ref alt gt seen="" i idx dosage allele a rc=0
    local -a alts gt_alleles
    _hx_tally Q
    [[ -n $region && -n $want_ref && -n $want_alt ]] \
        || hx_die "hx_gt_allele: usage: hx_gt_allele REGION REF ALT [RSID]"
    start=${region#*:}
    start=${start%%-*}
    want_ref=$(_hx_upper "$want_ref")
    want_alt=$(_hx_upper "$want_alt")
    out=$(hx_query "$region" '%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n') || rc=$?
    [[ $rc -eq 0 ]] || return "$rc"
    [[ -n $out ]] || return 0
    while IFS=$'\t' read -r pos id ref alt gt _; do
        [[ -n $pos ]] || continue
        seen="${seen:+$seen, }$ref/$alt"
        [[ $pos == "$start" ]] || continue
        [[ $(_hx_upper "$ref") == "$want_ref" ]] || continue
        if [[ -n $rsid && $id != '.' && $id != "$rsid" ]]; then
            continue
        fi
        idx=0
        if [[ $alt != '.' ]]; then
            local IFS_save=$IFS
            IFS=',' read -r -a alts <<<"$(_hx_upper "$alt")"
            IFS=$IFS_save
            i=0
            for a in "${alts[@]}"; do
                i=$((i + 1))
                if [[ $a == "$want_alt" ]]; then
                    idx=$i
                    break
                fi
            done
            [[ $idx -gt 0 ]] || continue
        fi
        # Matching record found: report its dosage or nothing for a missing call.
        _hx_gt_is_called "$gt" || return 0
        gt=$(_hx_norm_gt "$gt")
        dosage=0
        local IFS_save=$IFS
        IFS='/' read -r -a gt_alleles <<<"$gt"
        IFS=$IFS_save
        for allele in "${gt_alleles[@]}"; do
            [[ $idx -gt 0 && $allele == "$idx" ]] && dosage=$((dosage + 1))
        done
        _hx_tally F
        printf '%s\t%s\n' "$dosage" "$gt"
        return 0
    done <<<"$out"
    hx_warn "$region has REF/ALT ${seen:-?}, expected $want_ref>$want_alt - skipped"
    return 0
}

hx_alt_dosage() {
    local gt=${1:-} n=0 a
    local -a alleles
    if ! _hx_gt_is_called "$gt"; then
        printf 'NA\n'
        return 0
    fi
    gt=${gt//|//}
    local IFS_save=$IFS
    IFS='/' read -r -a alleles <<<"$gt"
    IFS=$IFS_save
    for a in "${alleles[@]}"; do
        if [[ ! $a =~ ^[0-9]+$ ]]; then
            printf 'NA\n'
            return 0
        fi
        [[ $a != 0 ]] && n=$((n + 1))
    done
    printf '%s\n' "$n"
}

hx_count() {
    local c
    c=$(grep -c -- "${2:-}" "${1:-/dev/null}" 2>/dev/null || true)
    c=${c%%$'\n'*}
    [[ $c =~ ^[0-9]+$ ]] || c=0
    printf '%s\n' "$c"
}

hx_report_coverage() {
    local q=${HX_QUERIED:-0} f=${HX_FOUND:-0}
    if [[ -r ${HX_COUNTER_FILE:-} ]]; then
        q=$(hx_count "$HX_COUNTER_FILE" '^Q$')
        f=$(hx_count "$HX_COUNTER_FILE" '^F$')
        rm -f "$HX_COUNTER_FILE"
    fi
    HX_QUERIED=$q
    HX_FOUND=$f
    printf 'HX_DATA_COVERAGE: %s/%s markers\n' "$f" "$q"
}
