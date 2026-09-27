#!/bin/bash
#===============================================================================
# STR EXPANSION ANALYSIS
# Analiza chorób powtórzeń (Huntington, Fragile X, ALS, etc.)
# Wymaga: ExpansionHunter
#===============================================================================

set -o pipefail

BAM="${1:-}"
# Referencja: ta sama, do ktorej wyrownano BAM (Helixight zaklada GRCh38;
# domyslna nazwa pliku jest historyczna). Wspolrzedne regionow STR ponizej
# pochodza z oryginalnego skryptu i nie zostaly zweryfikowane dla GRCh38.
REFERENCE="${2:-human_g1k_v37.fasta}"  # TODO(build): unverified GRCh38
OUTPUT_DIR="${3:-./str_analysis}"

# Skrypt wymaga wyrownanego BAM/CRAM (nie VCF!): sprawdz przed utworzeniem katalogu
[[ -f "$BAM" && "$BAM" =~ \.(bam|cram)$ ]] || {
    echo "str_analysis.sh requires an aligned BAM/CRAM file, got: ${BAM:-<none>}" >&2
    exit 2
}
samtools quickcheck "$BAM" || exit 2
# Contig naming of the BAM ("chr4" vs "4"): the fallback regions below are
# bare and get this prefix so `samtools depth -r` parses them on either style.
HX_CHR=$(samtools view -H "$BAM" | awk '$1=="@SQ" && !d {sub("SN:","",$2); print ($2 ~ /^chr/ ? "chr" : ""); d=1}')

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║         STR EXPANSION ANALYSIS                                    ║"
echo "║         Choroby powtórzeń tandemowych                             ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

mkdir -p "$OUTPUT_DIR"

#===============================================================================
# SPRAWDZENIE EXPANSIONHUNTER
#===============================================================================

if ! command -v ExpansionHunter &> /dev/null; then
    echo "⚠️  ExpansionHunter nie znaleziony!"
    echo ""
    echo "Instalacja na macOS:"
    echo ""
    echo "# Opcja 1: Homebrew (jeśli dostępny)"
    echo "brew install brewsci/bio/expansion-hunter"
    echo ""
    echo "# Opcja 2: Pobranie binarki"
    echo "wget https://github.com/Illumina/ExpansionHunter/releases/download/v5.0.0/ExpansionHunter-v5.0.0-macosx_x86_64.tar.gz"
    echo "tar xzf ExpansionHunter-v5.0.0-macosx_x86_64.tar.gz"
    echo "sudo mv ExpansionHunter-v5.0.0-macosx_x86_64/bin/ExpansionHunter /usr/local/bin/"
    echo ""
    echo "# Dla Apple Silicon (M1/M2/M3/M4) może wymagać Rosetta:"
    echo "softwareupdate --install-rosetta"
    echo ""
    
    # Spróbuj alternatywnej analizy
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Alternatywna analiza: Sprawdzam pokrycie regionów STR..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    # HTT - Huntington (chr4:3076604-3076693, CAG repeat)
    echo ""
    echo "=== HUNTINGTON (HTT) - CAG repeat ==="
    echo "Region: chr4:3076604-3076693"
    echo "Normalny: <36 powtórzeń | Patogenny: >40 powtórzeń"
    samtools depth -r "${HX_CHR}4:3076604-3076693" "$BAM" | \
        awk '{sum+=$3; cnt++} END {print "Średnie pokrycie regionu: " (cnt>0 ? sum/cnt : 0) "x"}' \
        || { echo "samtools depth failed for region ${HX_CHR}4:3076604-3076693 in $BAM (see the samtools message above: contig not in the BAM header / genome build mismatch, or missing .bai index)" >&2; exit 2; }
    echo ""
    
    # FMR1 - Fragile X (chrX:146993569-146993628, CGG repeat)
    echo "=== FRAGILE X (FMR1) - CGG repeat ==="
    echo "Region: chrX:146993569-146993628"
    echo "Normalny: <45 | Premutacja: 55-200 | Pełna mutacja: >200"
    samtools depth -r "${HX_CHR}X:146993569-146993628" "$BAM" | \
        awk '{sum+=$3; cnt++} END {print "Średnie pokrycie regionu: " (cnt>0 ? sum/cnt : 0) "x"}' \
        || { echo "samtools depth failed for region ${HX_CHR}X:146993569-146993628 in $BAM (see the samtools message above: contig not in the BAM header / genome build mismatch, or missing .bai index)" >&2; exit 2; }
    echo ""
    
    # C9orf72 - ALS/FTD (chr9:27573528-27573546, GGGGCC repeat)
    echo "=== ALS/FTD (C9orf72) - GGGGCC repeat ==="
    echo "Region: chr9:27573528-27573546"
    echo "Normalny: <30 | Patogenny: >60 (często setki-tysiące)"
    samtools depth -r "${HX_CHR}9:27573528-27573546" "$BAM" | \
        awk '{sum+=$3; cnt++} END {print "Średnie pokrycie regionu: " (cnt>0 ? sum/cnt : 0) "x"}' \
        || { echo "samtools depth failed for region ${HX_CHR}9:27573528-27573546 in $BAM (see the samtools message above: contig not in the BAM header / genome build mismatch, or missing .bai index)" >&2; exit 2; }
    echo ""
    
    # DMPK - Dystrofia miotoniczna typ 1 (chr19:46273463-46273524, CTG repeat)
    echo "=== DYSTROFIA MIOTONICZNA TYP 1 (DMPK) - CTG repeat ==="
    echo "Region: chr19:46273463-46273524"
    echo "Normalny: 5-37 | Patogenny: >50"
    samtools depth -r "${HX_CHR}19:46273463-46273524" "$BAM" | \
        awk '{sum+=$3; cnt++} END {print "Średnie pokrycie regionu: " (cnt>0 ? sum/cnt : 0) "x"}' \
        || { echo "samtools depth failed for region ${HX_CHR}19:46273463-46273524 in $BAM (see the samtools message above: contig not in the BAM header / genome build mismatch, or missing .bai index)" >&2; exit 2; }
    echo ""
    
    # FXN - Ataksja Friedreicha (chr9:71652202-71652240, GAA repeat)
    echo "=== ATAKSJA FRIEDREICHA (FXN) - GAA repeat ==="
    echo "Region: chr9:71652202-71652240"
    echo "Normalny: <33 | Patogenny: >66"
    samtools depth -r "${HX_CHR}9:71652202-71652240" "$BAM" | \
        awk '{sum+=$3; cnt++} END {print "Średnie pokrycie regionu: " (cnt>0 ? sum/cnt : 0) "x"}' \
        || { echo "samtools depth failed for region ${HX_CHR}9:71652202-71652240 in $BAM (see the samtools message above: contig not in the BAM header / genome build mismatch, or missing .bai index)" >&2; exit 2; }
    echo ""
    
    echo "⚠️  UWAGA: Analiza pokrycia NIE określa liczby powtórzeń!"
    echo "   Niskie pokrycie może sugerować ekspansję (trudne do sekwencjonowania)"
    echo "   Dla dokładnej analizy zainstaluj ExpansionHunter"
    echo ""
    
    exit 0
fi

#===============================================================================
# PEŁNA ANALIZA EXPANSIONHUNTER
#===============================================================================

echo "✅ ExpansionHunter znaleziony!"
echo ""

# Katalog z definicjami powtórzeń
CATALOG_DIR="$OUTPUT_DIR/catalogs"
mkdir -p "$CATALOG_DIR"

# Utwórz katalog najważniejszych chorób
cat > "$CATALOG_DIR/disease_strs.json" << 'EOF'
{
  "LocusId": "HTT",
  "LocusStructure": "(CAG)*CAACAG(CCG)*",
  "ReferenceRegion": "4:3076604-3076693",
  "VariantType": "Repeat",
  "Disease": "Huntington disease"
},
{
  "LocusId": "FMR1",
  "LocusStructure": "(CGG)*",
  "ReferenceRegion": "X:146993569-146993628",
  "VariantType": "Repeat",
  "Disease": "Fragile X syndrome"
},
{
  "LocusId": "C9orf72",
  "LocusStructure": "(GGGGCC)*",
  "ReferenceRegion": "9:27573528-27573546",
  "VariantType": "Repeat",
  "Disease": "ALS/FTD"
},
{
  "LocusId": "DMPK",
  "LocusStructure": "(CTG)*",
  "ReferenceRegion": "19:46273463-46273524",
  "VariantType": "Repeat",
  "Disease": "Myotonic dystrophy type 1"
},
{
  "LocusId": "FXN",
  "LocusStructure": "(GAA)*",
  "ReferenceRegion": "9:71652202-71652240",
  "VariantType": "Repeat",
  "Disease": "Friedreich ataxia"
},
{
  "LocusId": "ATXN1",
  "LocusStructure": "(CAG)*",
  "ReferenceRegion": "6:16327865-16327954",
  "VariantType": "Repeat",
  "Disease": "Spinocerebellar ataxia type 1"
},
{
  "LocusId": "ATXN2",
  "LocusStructure": "(CAG)*",
  "ReferenceRegion": "12:112036754-112036823",
  "VariantType": "Repeat",
  "Disease": "Spinocerebellar ataxia type 2"
},
{
  "LocusId": "ATXN3",
  "LocusStructure": "(CAG)*",
  "ReferenceRegion": "14:92537355-92537387",
  "VariantType": "Repeat",
  "Disease": "Machado-Joseph disease (SCA3)"
},
{
  "LocusId": "AR",
  "LocusStructure": "(CAG)*",
  "ReferenceRegion": "X:66765159-66765227",
  "VariantType": "Repeat",
  "Disease": "Spinal and bulbar muscular atrophy (Kennedy disease)"
}
EOF

echo "Uruchamiam ExpansionHunter..."

# Uruchom ExpansionHunter
ExpansionHunter \
    --reads "$BAM" \
    --reference "$REFERENCE" \
    --variant-catalog "$CATALOG_DIR/disease_strs.json" \
    --output-prefix "$OUTPUT_DIR/str_results" \
    2> "$OUTPUT_DIR/expansionhunter.log"

#===============================================================================
# ANALIZA WYNIKÓW
#===============================================================================

if [ -f "$OUTPUT_DIR/str_results.vcf" ]; then
    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════╗"
    echo "║                    WYNIKI ANALIZY STR                             ║"
    echo "╚═══════════════════════════════════════════════════════════════════╝"
    echo ""
    
    {
        echo "RAPORT ANALIZY STR (Short Tandem Repeats)"
        echo "Data: $(date)"
        echo "========================================"
        echo ""
        
        # Parsuj VCF
        while read -r line; do
            [[ "$line" == "#"* ]] && continue
            
            CHROM=$(echo "$line" | cut -f1)
            POS=$(echo "$line" | cut -f2)
            ID=$(echo "$line" | cut -f3)
            INFO=$(echo "$line" | cut -f8)
            
            # Wyciągnij liczbę powtórzeń
            REPCN=$(echo "$INFO" | grep -oP 'REPCN=\K[^;]+' || echo "N/A")
            
            echo "=== $ID ==="
            echo "Pozycja: $CHROM:$POS"
            echo "Liczba powtórzeń: $REPCN"
            
            # Interpretacja per locus
            case $ID in
                HTT)
                    echo "Choroba: Huntington"
                    echo "Normalny: <36 | Pośredni: 36-39 | Patogenny: ≥40"
                    if [[ "$REPCN" =~ ^[0-9]+$ ]] && [ "$REPCN" -lt 36 ]; then
                        echo "Status: ✅ NORMALNY"
                    elif [[ "$REPCN" =~ ^[0-9]+$ ]] && [ "$REPCN" -ge 40 ]; then
                        echo "Status: 🔴 PATOGENNY - skonsultuj z genetykiem!"
                    else
                        echo "Status: ⚠️ Sprawdź ręcznie"
                    fi
                    ;;
                FMR1)
                    echo "Choroba: Fragile X"
                    echo "Normalny: <45 | Premutacja: 55-200 | Pełna mutacja: >200"
                    ;;
                C9orf72)
                    echo "Choroba: ALS/FTD"
                    echo "Normalny: <30 | Patogenny: >60"
                    ;;
                DMPK)
                    echo "Choroba: Dystrofia miotoniczna typ 1"
                    echo "Normalny: 5-37 | Patogenny: >50"
                    ;;
                FXN)
                    echo "Choroba: Ataksja Friedreicha"
                    echo "Normalny: <33 | Nosiciel: 33-65 | Patogenny: >66"
                    ;;
            esac
            
            echo ""
            
        done < "$OUTPUT_DIR/str_results.vcf"
        
    } | tee "$OUTPUT_DIR/str_report.txt"
    
    echo ""
    echo "Pliki wynikowe:"
    echo "  • $OUTPUT_DIR/str_results.vcf"
    echo "  • $OUTPUT_DIR/str_results.json"
    echo "  • $OUTPUT_DIR/str_report.txt"
    
else
    echo "❌ Błąd: Nie wygenerowano pliku wynikowego"
    echo "Sprawdź log: $OUTPUT_DIR/expansionhunter.log"
fi
