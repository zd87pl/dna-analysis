#!/bin/bash
#===============================================================================
# EKSTRAKCJA mtDNA dla Haplogrep
#===============================================================================

# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"
set -eo pipefail

VCF="${1:-}"
hx_require_arg "$VCF" "$(basename "$0") <variants.vcf.gz> [output_dir]"
OUTPUT_DIR="${2:-mtdna_analysis}"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║         EKSTRAKCJA mtDNA DLA HAPLOGREP                            ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

# Zindeksowana, posortowana kopia robocza z ujednoliconymi nazwami kontigów
# (chrM/chrMT/M -> MT); zwykły .vcf / .vcf.gz jest konwertowany, nie odrzucany.
hx_prepare_vcf "$VCF"
VCF="$HX_VCF"

#===============================================================================
# 1. Ekstrakcja VCF tylko dla MT
#===============================================================================

echo "[1/3] Ekstrakcja wariantów mtDNA do VCF..."

# Kontig MT wykrywany z indeksu (tylko kontigi z rekordami); po normalizacji
# przez bibliotekę nazwa to "MT". Wyzerowanie chroni przed MT_CHROM z otoczenia.
MT_CHROM=""
MT_CHROM=$(bcftools index -s "$VCF" | cut -f1 | grep -Exm1 'MT|chrM|chrMT|M' || true)

if [ -z "$MT_CHROM" ]; then
    echo "⚠️ Brak wariantów mtDNA w pliku (exome/chip?)"
    echo "⚠️ Brak wariantów mtDNA w pliku (exome/chip?)" >&2
    # The MT contig is the one queried "marker"; none found here.
    # shellcheck disable=SC2034  # read by hx_report_coverage (scripts/lib/vcf_common.sh)
    HX_QUERIED=1
    HX_FOUND=0
    hx_report_coverage
    exit 0
fi

mkdir -p "$OUTPUT_DIR"

COUNT=$(bcftools view -r "$MT_CHROM" -H "$VCF" | wc -l | tr -d ' ')
echo "  Znaleziono $COUNT wariantów jako '$MT_CHROM'"

# Eksportuj do VCF
bcftools view -r "$MT_CHROM" "$VCF" -Oz -o "$OUTPUT_DIR/mtdna_variants.vcf.gz"
bcftools view -r "$MT_CHROM" "$VCF" -o "$OUTPUT_DIR/mtdna_variants.vcf"

echo "  ✓ Zapisano: $OUTPUT_DIR/mtdna_variants.vcf"

#===============================================================================
# 2. Format HSD (Haplogrep Sample Data)
#===============================================================================

echo ""
echo "[2/3] Tworzenie pliku HSD..."

# Format HSD: SampleID	Range	Haplogroup	Polymorphisms
# Polymorphisms to lista wariantów oddzielona tabulatorem

# Zbierz warianty w formacie Haplogrep (pozycja + zmiana)
VARIANTS=$(bcftools view -r "$MT_CHROM" -H "$VCF" | \
    awk '{
        pos = $2
        ref = $4
        alt = $5
        
        # Haplogrep format: pozycjaZmiana (np. 73G, 263G, 7028T)
        # Dla SNP: pozycja + alternatywny allel
        # Dla insercji: pozycja.1Xbaza (np. 309.1C)
        # Dla delecji: pozycjad (np. 522-523d)
        
        if (length(ref) == 1 && length(alt) == 1) {
            # SNP
            printf "%d%s\t", pos, alt
        } else if (length(alt) > length(ref)) {
            # Insercja
            ins_base = substr(alt, length(ref)+1)
            for (i = 1; i <= length(ins_base); i++) {
                printf "%d.%d%s\t", pos, i, substr(ins_base, i, 1)
            }
        } else if (length(ref) > length(alt)) {
            # Delecja
            printf "%d-%dd\t", pos+1, pos+length(ref)-1
        }
    }')

# Zapisz HSD
cat > "$OUTPUT_DIR/mtdna_haplogrep.hsd" << EOF
SampleID	Range	Haplogroup	Polymorphisms
Sample1	1-16569		$VARIANTS
EOF

echo "  ✓ Zapisano: $OUTPUT_DIR/mtdna_haplogrep.hsd"

#===============================================================================
# 3. Format tekstowy (do wklejenia)
#===============================================================================

echo ""
echo "[3/3] Lista wariantów do wklejenia..."

# Prosta lista wariantów
bcftools view -r "$MT_CHROM" -H "$VCF" | \
    awk '{
        pos = $2
        ref = $4
        alt = $5
        
        if (length(ref) == 1 && length(alt) == 1) {
            print pos alt
        } else if (length(alt) > length(ref)) {
            ins = substr(alt, length(ref)+1)
            for (i = 1; i <= length(ins); i++) {
                print pos "." i substr(ins, i, 1)
            }
        } else {
            print pos "-" (pos+length(ref)-1) "d"
        }
    }' > "$OUTPUT_DIR/mtdna_variants_list.txt"

echo "  ✓ Zapisano: $OUTPUT_DIR/mtdna_variants_list.txt"

#===============================================================================
# PODSUMOWANIE
#===============================================================================

echo ""
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║                         GOTOWE!                                   ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "Pliki do użycia z Haplogrep (https://haplogrep.i-med.ac.at/):"
echo ""
echo "📁 OPCJA 1: Upload pliku VCF"
echo "   Plik: $OUTPUT_DIR/mtdna_variants.vcf"
echo ""
echo "📁 OPCJA 2: Upload pliku HSD"
echo "   Plik: $OUTPUT_DIR/mtdna_haplogrep.hsd"
echo ""
echo "📋 OPCJA 3: Wklej listę wariantów (najprostsze!)"
echo "   Na stronie wybierz 'Paste Data' i wklej:"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
tr '\n' ' ' < "$OUTPUT_DIR/mtdna_variants_list.txt"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "✅ Znaleziono $COUNT wariantów mtDNA - pliki dla Haplogrep zapisane w $OUTPUT_DIR/"
echo ""
# The MT contig is the one queried "marker"; found when it carries records.
# shellcheck disable=SC2034  # read by hx_report_coverage (scripts/lib/vcf_common.sh)
HX_QUERIED=1
# shellcheck disable=SC2034  # read by hx_report_coverage
HX_FOUND=$([[ ${COUNT:-0} -gt 0 ]] && echo 1 || echo 0)
hx_report_coverage
