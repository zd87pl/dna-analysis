#!/bin/bash
#===============================================================================
# CLINVAR PATHOGENIC VARIANT SCANNER
# Znajduje wszystkie patogenne warianty w Twoim genomie
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"
OUTDIR="clinvar_analysis"
mkdir -p "$OUTDIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║            CLINVAR PATHOGENIC VARIANT SCANNER                             ║"
echo "║            Szukam patogennych wariantow w Twoim genomie                   ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Plik VCF: $VCF"
echo "Wyniki: $OUTDIR/"
echo ""

#===============================================================================
# KROK 1: Pobierz ClinVar (jesli brak)
#===============================================================================

CLINVAR_VCF="$OUTDIR/clinvar_GRCh38.vcf.gz"

if [ ! -f "$CLINVAR_VCF" ]; then
    echo "📥 Pobieram baze ClinVar (GRCh38)..."
    echo "   To moze potrwac kilka minut (~70MB)..."
    
    wget -q --show-progress -O "$CLINVAR_VCF" \
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
    
    wget -q -O "${CLINVAR_VCF}.tbi" \
        "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
    
    if [ -f "$CLINVAR_VCF" ]; then
        echo "✅ ClinVar pobrany!"
    else
        echo "❌ Blad pobierania ClinVar!"
        exit 1
    fi
else
    echo "✅ ClinVar juz pobrany, uzywam istniejacego pliku"
fi
echo ""

#===============================================================================
# KROK 2: Wyodrebnij tylko patogenne z ClinVar
#===============================================================================

echo "🔬 Wyodrebniam warianty Pathogenic/Likely_pathogenic z ClinVar..."

CLINVAR_PATH_VCF="$OUTDIR/clinvar_pathogenic_only.vcf.gz"

if [ ! -f "$CLINVAR_PATH_VCF" ]; then
    bcftools view -i 'INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic"' \
        "$CLINVAR_VCF" -Oz -o "$CLINVAR_PATH_VCF" 2>/dev/null
    bcftools index -f "$CLINVAR_PATH_VCF" 2>/dev/null
    
    PATHOGENIC_COUNT=$(bcftools view -H "$CLINVAR_PATH_VCF" 2>/dev/null | wc -l)
    echo "   ClinVar zawiera $PATHOGENIC_COUNT wariantow Pathogenic/Likely_pathogenic"
else
    echo "   Uzywam istniejacego pliku pathogenic"
fi
echo ""

#===============================================================================
# KROK 3: Znajdz wspolne warianty (bcftools isec)
#===============================================================================

echo "🔍 Porownuje Twoj genom z ClinVar Pathogenic..."
echo ""

mkdir -p "$OUTDIR/isec_results"

# Upewnij sie ze VCF jest zindeksowany
if [ ! -f "${VCF}.tbi" ] && [ ! -f "${VCF}.csi" ]; then
    echo "   Indeksuje Twoj VCF..."
    bcftools index -f "$VCF" 2>/dev/null
fi

# Znajdz przeciecie - warianty ktore sa ZAROWNO w Twoim VCF JAK I w ClinVar Pathogenic
bcftools isec -p "$OUTDIR/isec_results" -n=2 -w1 \
    "$VCF" "$CLINVAR_PATH_VCF" 2>/dev/null

#===============================================================================
# KROK 4: Analizuj wyniki
#===============================================================================

RESULTS_FILE="$OUTDIR/isec_results/0000.vcf"
REPORT="$OUTDIR/PATHOGENIC_VARIANTS_FOUND.txt"

if [ -f "$RESULTS_FILE" ]; then
    FOUND_COUNT=$(grep -v "^#" "$RESULTS_FILE" 2>/dev/null | wc -l)
    
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    if [ $FOUND_COUNT -eq 0 ]; then
        echo "✅ Nie znaleziono wariantow patogennych w Twoim genomie!"
        echo ""
        echo "   To DOBRA wiadomosc - zadne ze znanych wariantow patogennych"
        echo "   z ClinVar nie zostaly wykryte w Twoim VCF."
        
        cat << EOF > "$REPORT"
================================================================================
                RAPORT CLINVAR - WARIANTY PATOGENNE
                Wygenerowano: $(date)
================================================================================

WYNIK: NIE ZNALEZIONO wariantow patogennych

Twoj genom zostal porownany z baza ClinVar zawierajaca wszystkie znane
warianty sklasyfikowane jako Pathogenic lub Likely_pathogenic.

Nie znaleziono zadnych dopasowań.

UWAGA: To nie oznacza braku JAKIEKOLWIEK ryzyka genetycznego.
- ClinVar nie zawiera wszystkich patogennych wariantow
- Niektore warianty moga nie byc jeszcze sklasyfikowane
- Analiza nie obejmuje CNV, STR, ani wariantow strukturalnych

================================================================================
EOF
    else
        echo "🔴 ZNALEZIONO $FOUND_COUNT wariantow obecnych w ClinVar Pathogenic!"
        echo ""
        
        cat << EOF > "$REPORT"
================================================================================
                RAPORT CLINVAR - WARIANTY PATOGENNE
                Wygenerowano: $(date)
================================================================================

⚠️ ZNALEZIONO $FOUND_COUNT WARIANTOW W CLINVAR PATHOGENIC

SZCZEGOLY PONIZEJ:
================================================================================

EOF
        
        echo "   Szczegoly znalezionych wariantow:"
        echo "   ─────────────────────────────────────────"
        
        grep -v "^#" "$RESULTS_FILE" | head -50 | while read -r line; do
            CHROM=$(echo "$line" | cut -f1)
            POS=$(echo "$line" | cut -f2)
            RSID=$(echo "$line" | cut -f3)
            REF=$(echo "$line" | cut -f4)
            ALT=$(echo "$line" | cut -f5)
            GT=$(echo "$line" | cut -f10 | cut -d: -f1)
            
            # Pobierz info z ClinVar
            CLINVAR_INFO=$(bcftools query -r "${CHROM}:${POS}" \
                -f '%ID\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\n' \
                "$CLINVAR_PATH_VCF" 2>/dev/null | head -1)
            
            CV_ID=$(echo "$CLINVAR_INFO" | cut -f1)
            CLNSIG=$(echo "$CLINVAR_INFO" | cut -f2)
            CLNDN=$(echo "$CLINVAR_INFO" | cut -f3 | sed 's/_/ /g')
            REVIEW=$(echo "$CLINVAR_INFO" | cut -f4)
            
            # Genotyp
            case $GT in
                "0/1"|"0|1"|"1/0"|"1|0") GENO="HETEROZYGOTA" ;;
                "1/1"|"1|1") GENO="HOMOZYGOTA" ;;
                *) GENO="$GT" ;;
            esac
            
            echo ""
            echo "   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "   📍 chr${CHROM}:${POS} ${REF}>${ALT}"
            echo "      ClinVar ID: $CV_ID"
            echo "      Klasyfikacja: $CLNSIG"
            echo "      Choroba: $CLNDN"
            echo "      Twoj genotyp: $GENO"
            echo "      Review status: $REVIEW"
            
            # Zapisz do raportu
            cat << EOF >> "$REPORT"
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WARIANT: chr${CHROM}:${POS} ${REF}>${ALT}
ClinVar ID: $CV_ID
Klasyfikacja: $CLNSIG
Choroba: $CLNDN
Twoj genotyp: $GENO
Review status: $REVIEW
Link: https://www.ncbi.nlm.nih.gov/clinvar/?term=${CV_ID}

EOF
        done
        
        if [ $FOUND_COUNT -gt 50 ]; then
            echo ""
            echo "   ... i $(($FOUND_COUNT - 50)) wiecej (zobacz pelny raport)"
        fi
    fi
    
    echo ""
else
    echo "⚠️ Blad w porownaniu - sprawdz czy VCF jest poprawny"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📁 Wyniki zapisane:"
echo "   $REPORT"
echo "   $OUTDIR/isec_results/"
echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║  ⚠️ UWAGA: Wyniki wymagaja interpretacji przez genetyka klinicznego!      ║"
echo "║  Nie wszystkie warianty 'pathogenic' powoduja chorobe u kazdego.          ║"
echo "║  Penetracja, ekspresja i czynniki srodowiskowe maja znaczenie.            ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
