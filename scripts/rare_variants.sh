#!/bin/bash
#===============================================================================
#
#    ██████╗  █████╗ ██████╗ ███████╗    ██╗   ██╗ █████╗ ██████╗ ██╗ █████╗ ███╗   ██╗████████╗███████╗
#    ██╔══██╗██╔══██╗██╔══██╗██╔════╝    ██║   ██║██╔══██╗██╔══██╗██║██╔══██╗████╗  ██║╚══██╔══╝██╔════╝
#    ██████╔╝███████║██████╔╝█████╗      ██║   ██║███████║██████╔╝██║███████║██╔██╗ ██║   ██║   ███████╗
#    ██╔══██╗██╔══██║██╔══██╗██╔══╝      ╚██╗ ██╔╝██╔══██║██╔══██╗██║██╔══██║██║╚██╗██║   ██║   ╚════██║
#    ██║  ██║██║  ██║██║  ██║███████╗     ╚████╔╝ ██║  ██║██║  ██║██║██║  ██║██║ ╚████║   ██║   ███████║
#    ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝      ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝
#
#    ANALIZA RZADKICH WARIANTOW
#    Jak unikalne jest Twoje DNA?
#
#===============================================================================

# Shared VCF helpers (hx_prepare_vcf, hx_count, hx_report_coverage, ...)
# shellcheck source=lib/vcf_common.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/lib/vcf_common.sh"

hx_require_arg "${1:-}" "$(basename "$0") <variants.vcf.gz>"
VCF="$1"
OUTDIR="rare_variants_analysis"

# Build (or reuse) a bgzip+indexed, bare-contig, GRCh38-checked working copy;
# every bcftools call below reads $HX_VCF.
hx_prepare_vcf "$VCF"
mkdir -p "$OUTDIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    ANALIZA RZADKICH WARIANTOW                             ║"
echo "║              Jak unikalne jest Twoje DNA vs populacja?                    ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "Plik VCF: $VCF"
echo "Wyniki: $OUTDIR/"
echo ""

#===============================================================================
# PODSTAWOWE STATYSTYKI
#===============================================================================

echo "📊 PODSTAWOWE STATYSTYKI TWOJEGO GENOMU"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Calkowita liczba wariantow
TOTAL_VARIANTS=$(bcftools view -H "$HX_VCF" | wc -l)
echo "   Calkowita liczba wariantow: $TOTAL_VARIANTS"

# SNP vs Indel
SNPS=$(bcftools view -v snps -H "$HX_VCF" | wc -l)
INDELS=$(bcftools view -v indels -H "$HX_VCF" | wc -l)
echo "   SNP: $SNPS"
echo "   Indele: $INDELS"

# Heterozygoty vs homozygoty
HETS=$(bcftools view -H "$HX_VCF" | awk -F'\t' '{print $NF}' | grep -c "0/1\|0|1\|1/0\|1|0")
HOMS=$(bcftools view -H "$HX_VCF" | awk -F'\t' '{print $NF}' | grep -c "1/1\|1|1")
echo "   Heterozygoty: $HETS"
echo "   Homozygoty ALT: $HOMS"

# awk instead of bc (no extra dependency; bc prints nothing useful on x/0)
if [ "$HOMS" -gt 0 ]; then
    HET_HOM_RATIO=$(awk -v h="$HETS" -v m="$HOMS" 'BEGIN{printf "%.2f", h/m}')
else
    HET_HOM_RATIO="N/A"
fi
echo "   Stosunek Het/Hom: $HET_HOM_RATIO (norma: 1.5-2.0)"

echo ""

#===============================================================================
# ANALIZA CZESTOSCI ALLELICZNYCH (na podstawie rsID)
#===============================================================================

echo "📊 ROZKLAD CZESTOSCI WARIANTOW (szacunek)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Analizuje warianty pod katem czestosci w populacji..."
echo "(Na podstawie struktury VCF - pelna analiza wymaga bazy gnomAD)"
echo ""

# Policz warianty z rsID vs bez
VARIANTS_WITH_RS=$(bcftools query -f '%ID\n' "$HX_VCF" | grep -c "^rs")
VARIANTS_NO_RS=$((TOTAL_VARIANTS - VARIANTS_WITH_RS))

echo "   Warianty ze znanym rsID (w dbSNP): $VARIANTS_WITH_RS"
echo "   Warianty BEZ rsID (potencjalnie rzadkie/nowe): $VARIANTS_NO_RS"
echo ""

if [ "$TOTAL_VARIANTS" -gt 0 ]; then
    PERCENT_NOVEL=$(awk -v n="$VARIANTS_NO_RS" -v t="$TOTAL_VARIANTS" 'BEGIN{printf "%.1f%%", n*100/t}')
else
    PERCENT_NOVEL="N/A"
fi
echo "   Procent potencjalnie rzadkich/nowych: ${PERCENT_NOVEL}"
echo ""

# Kategoryzacja (szacunkowa)
echo "   📈 SZACUNKOWA KATEGORYZACJA:"
echo ""
echo "   Kategoria                  | Typowy genom | Twoj genom (szac.)"
echo "   ─────────────────────────────────────────────────────────────"
echo "   Czeste (MAF > 5%)          | ~3,500,000   | ~$((VARIANTS_WITH_RS * 80 / 100))"
echo "   Srednie (MAF 1-5%)         | ~500,000     | ~$((VARIANTS_WITH_RS * 15 / 100))"
echo "   Rzadkie (MAF 0.1-1%)       | ~150,000     | ~$((VARIANTS_WITH_RS * 4 / 100))"
echo "   Bardzo rzadkie (MAF <0.1%) | ~50,000      | ~$((VARIANTS_WITH_RS * 1 / 100))"
echo "   Potencjalnie nowe          | ~10,000      | ~$VARIANTS_NO_RS"
echo ""

#===============================================================================
# UNIKALNE WARIANTY - ANALIZA
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    ANALIZA UNIKALNOSCI GENOMU                             ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Singletonty (warianty obecne tylko u Ciebie w typowej kohorcie)
echo "🔬 POTENCJALNE SINGLETONTY (warianty bez rsID):"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Zapisz warianty bez rsID do pliku
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%GT]\n' "$HX_VCF" | \
    grep -v "^#" | \
    awk '$5 == "." {print}' | \
    head -1000 > "$OUTDIR/novel_variants.txt"

NOVEL_COUNT=$(wc -l < "$OUTDIR/novel_variants.txt")
echo "   Znaleziono $NOVEL_COUNT potencjalnie nowych wariantow (bez rsID)"
echo ""

# Rozklad po chromosomach
echo "   Rozklad nowych wariantow po chromosomach:"
echo ""
for chr in {1..22} X Y; do
    # hx_count always prints exactly one number (grep -c || echo 0 gave "0\n0")
    count=$(hx_count "$OUTDIR/novel_variants.txt" "^${chr}\s\|^chr${chr}\s")
    if [ "$count" -gt 0 ]; then
        bar=$(printf '%*s' "$((count / 50))" '' | tr ' ' '█')
        printf "   chr%-2s: %6d %s\n" "$chr" "$count" "$bar"
    fi
done

echo ""

#===============================================================================
# RZADKIE WARIANTY W WAZNYCH GENACH
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║              RZADKIE WARIANTY W KLINICZNYCH GENACH                        ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "Sprawdzam regiony waznych genow pod katem rzadkich wariantow..."
echo ""

# Lista waznych genow i ich lokalizacji (GRCh38)
CLINICAL_GENES="
BRCA1|17:43044295-43125483
BRCA2|13:32315474-32400266
TP53|17:7668402-7687538
PTEN|10:87863625-87971930
APC|5:112707498-112846239
MLH1|3:36993350-37050918
MSH2|2:47403067-47709917
ATM|11:108222484-108369102
CHEK2|22:28687743-28742422
PALB2|16:23603160-23641310
LDLR|19:11089463-11133820
APOB|2:21001429-21044073
PCSK9|1:55039548-55064852
MYBPC3|11:47331336-47352689
MYH7|14:23412738-23435706
KCNQ1|11:2444990-2849109
SCN5A|3:38548062-38691163
RYR2|1:237042491-237833588
LMNA|1:156052364-156109878
PKP2|12:32824653-32911934
TTN|2:178525989-178830802
DMD|X:31097677-33339609
CFTR|7:117287120-117715971
HFE|6:26087281-26098343
"

: > "$OUTDIR/clinical_genes_variants.txt"

echo "   GEN         | WARIANTY | BEZ rsID | KOMENTARZ"
echo "   ───────────────────────────────────────────────────────────────"

while IFS='|' read -r gene region; do
    [ -z "$gene" ] && continue
    
    # Policz warianty w tym genie ($HX_VCF has bare contig names, so the old
    # 'chr' retry is no longer needed)
    total=$(bcftools view -r "$region" -H "$HX_VCF" | wc -l)

    if [ "$total" -gt 0 ]; then
        # Policz bez rsID (grep -c always prints exactly one number; the old
        # `|| echo "0"` appended a second "0")
        no_rs=$(bcftools query -r "$region" -f '%ID\n' "$HX_VCF" | grep -c "^\.$")

        comment=""
        if [ "$no_rs" -gt 5 ]; then
            comment="⚠️ Wiele rzadkich!"
        elif [ "$no_rs" -gt 0 ]; then
            comment="Kilka rzadkich"
        fi
        
        printf "   %-12s | %8d | %8s | %s\n" "$gene" "$total" "$no_rs" "$comment"
        
        # Zapisz szczegoly do pliku
        bcftools query -r "$region" -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "$HX_VCF" | \
            awk -v gene="$gene" '{print gene"\t"$0}' >> "$OUTDIR/clinical_genes_variants.txt"
    fi
    
done <<< "$CLINICAL_GENES"

echo ""

#===============================================================================
# WARIANTY EXONOWE BEZ RSID
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                  POTENCJALNIE FUNKCJONALNE NOWE WARIANTY                  ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "Szukam wariantow bez rsID w regionach kodujacych..."
echo "(To sa najbardziej prawdopodobne do bycia funkcjonalnymi)"
echo ""

# Szukaj specyficznych pozycji znanych genow
echo "   Przykladowe potencjalnie nowe warianty w waznych genach:"
echo ""

# Sprawdz kilka kluczowych genow pod katem nowych wariantow
FOUND_NOVEL=0

for gene_info in "BRCA1|17:43044295-43125483" "BRCA2|13:32315474-32400266" "TP53|17:7668402-7687538"; do
    gene=$(echo "$gene_info" | cut -d'|' -f1)
    region=$(echo "$gene_info" | cut -d'|' -f2)
    
    # Znajdz warianty bez rsID
    novel=$(bcftools query -r "$region" -f '%CHROM:%POS\t%REF>%ALT\t%ID\n' "$HX_VCF" | grep "\.$" | head -3)

    if [ -n "$novel" ]; then
        echo "   $gene:"
        echo "$novel" | while read line; do
            pos=$(echo "$line" | cut -f1)
            change=$(echo "$line" | cut -f2)
            echo "     → $pos $change (BEZ rsID - potencjalnie nowy!)"
            FOUND_NOVEL=$((FOUND_NOVEL + 1))
        done
        echo ""
    fi
done

echo ""

#===============================================================================
# Ti/Tv RATIO - JAKOSC WARIANTOW
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    JAKOSC WARIANTOW (Ti/Tv)                               ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Policz tranzycje i transwersje
# Tranzycje: A<>G, C<>T
# Transwersje: wszystkie inne

TI=$(bcftools query -f '%REF %ALT\n' "$HX_VCF" | \
    awk '(($1=="A" && $2=="G") || ($1=="G" && $2=="A") ||
          ($1=="C" && $2=="T") || ($1=="T" && $2=="C")) {count++}
         END {print count+0}')

TV=$(bcftools query -f '%REF %ALT\n' "$HX_VCF" | \
    awk '!(($1=="A" && $2=="G") || ($1=="G" && $2=="A") ||
           ($1=="C" && $2=="T") || ($1=="T" && $2=="C")) &&
         length($1)==1 && length($2)==1 {count++}
         END {print count+0}')

TITV_RATIO="N/A"
if [ -n "$TI" ] && [ -n "$TV" ] && [ "$TV" -gt 0 ]; then
    # awk instead of bc: no extra dependency
    TITV_RATIO=$(awk -v t="$TI" -v v="$TV" 'BEGIN{printf "%.2f", t/v}')
    echo "   Tranzycje (Ti): $TI"
    echo "   Transwersje (Tv): $TV"
    echo "   Stosunek Ti/Tv: $TITV_RATIO"
    echo ""

    if awk -v r="$TITV_RATIO" 'BEGIN{exit !(r >= 2.0 && r <= 2.2)}'; then
        echo "   ✅ Ti/Tv w normie dla WGS (2.0-2.2)"
    elif awk -v r="$TITV_RATIO" 'BEGIN{exit !(r >= 2.8 && r <= 3.3)}'; then
        echo "   ✅ Ti/Tv w normie dla WES/exome (2.8-3.3)"
    elif awk -v r="$TITV_RATIO" 'BEGIN{exit !(r < 2.0)}'; then
        echo "   ⚠️ Ti/Tv nizsze niz oczekiwane - moze wskazywac na wiecej falszywie pozytywnych"
    fi
else
    echo "   Nie udalo sie obliczyc Ti/Tv"
fi

echo ""

#===============================================================================
# PODSUMOWANIE
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                              PODSUMOWANIE                                 ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

echo "📊 TWOJ GENOM W LICZBACH:"
echo ""
echo "   ┌─────────────────────────────────────────────────────────────┐"
echo "   │ Calkowite warianty:    $(printf '%12d' $TOTAL_VARIANTS)                     │"
echo "   │ SNP:                   $(printf '%12d' $SNPS)                     │"
echo "   │ Indele:                $(printf '%12d' $INDELS)                     │"
echo "   │ Warianty ze znanym ID: $(printf '%12d' $VARIANTS_WITH_RS)                     │"
echo "   │ Potencjalnie nowe:     $(printf '%12d' $VARIANTS_NO_RS)                     │"
echo "   │ Ti/Tv ratio:           $(printf '%12s' $TITV_RATIO)                     │"
echo "   └─────────────────────────────────────────────────────────────┘"
echo ""

echo "💡 INTERPRETACJA:"
echo ""
echo "   • Kazdy czlowiek ma ~4-5 milionow wariantow vs genom referencyjny"
echo "   • ~10,000-50,000 z nich to warianty rzadkie (MAF < 0.1%)"
echo "   • ~50-100 to warianty de novo (nowe mutacje)"
echo "   • Twoj genom zawiera $VARIANTS_NO_RS wariantow bez znanego rsID"
echo "     → To nie znaczy, ze sa chorobowe - wiekszosc jest neutralna!"
echo ""
echo "   ✅ Przeanalizowano $TOTAL_VARIANTS wariantow: $VARIANTS_NO_RS bez rsID (potencjalnie rzadkich/nowych)"
echo ""

echo "📁 PLIKI WYNIKOWE:"
echo "   • $OUTDIR/novel_variants.txt - warianty bez rsID"
echo "   • $OUTDIR/clinical_genes_variants.txt - warianty w genach klinicznych"
echo ""
# This script summarises whole-file statistics rather than fixed markers: the
# input record set counts as one queried "marker", found when the VCF carries
# at least one record (a header-only file yields 0/1 = "no usable data").
# shellcheck disable=SC2034  # read by hx_report_coverage (scripts/lib/vcf_common.sh)
HX_QUERIED=1
# shellcheck disable=SC2034  # read by hx_report_coverage
HX_FOUND=$([[ ${TOTAL_VARIANTS:-0} -gt 0 ]] && echo 1 || echo 0)
hx_report_coverage
