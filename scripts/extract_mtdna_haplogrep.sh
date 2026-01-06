#!/bin/bash
#===============================================================================
# EKSTRAKCJA mtDNA dla Haplogrep
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"
OUTPUT_DIR="${2:-.}"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║         EKSTRAKCJA mtDNA DLA HAPLOGREP                            ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

mkdir -p "$OUTPUT_DIR"

#===============================================================================
# 1. Ekstrakcja VCF tylko dla MT
#===============================================================================

echo "[1/3] Ekstrakcja wariantów mtDNA do VCF..."

# Próbuj różne nazwy chromosomu MT
for MT_NAME in "MT" "chrM" "chrMT" "M"; do
    COUNT=$(bcftools view -r "$MT_NAME" -H "$VCF" 2>/dev/null | wc -l | tr -d ' ')
    if [ "$COUNT" -gt "0" ]; then
        echo "  Znaleziono $COUNT wariantów jako '$MT_NAME'"
        
        # Eksportuj do VCF
        bcftools view -r "$MT_NAME" "$VCF" -Oz -o "$OUTPUT_DIR/mtdna_variants.vcf.gz" 2>/dev/null
        bcftools view -r "$MT_NAME" "$VCF" -o "$OUTPUT_DIR/mtdna_variants.vcf" 2>/dev/null
        
        MT_CHROM="$MT_NAME"
        break
    fi
done

if [ -z "$MT_CHROM" ]; then
    echo "❌ Nie znaleziono wariantów mtDNA!"
    exit 1
fi

echo "  ✓ Zapisano: $OUTPUT_DIR/mtdna_variants.vcf"

#===============================================================================
# 2. Format HSD (Haplogrep Sample Data)
#===============================================================================

echo ""
echo "[2/3] Tworzenie pliku HSD..."

# Format HSD: SampleID	Range	Haplogroup	Polymorphisms
# Polymorphisms to lista wariantów oddzielona tabulatorem

# Zbierz warianty w formacie Haplogrep (pozycja + zmiana)
VARIANTS=$(bcftools view -r "$MT_CHROM" -H "$VCF" 2>/dev/null | \
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
bcftools view -r "$MT_CHROM" -H "$VCF" 2>/dev/null | \
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
cat "$OUTPUT_DIR/mtdna_variants_list.txt" | tr '\n' ' '
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
