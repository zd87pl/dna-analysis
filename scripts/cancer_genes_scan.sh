#!/bin/bash
#===============================================================================
# KOMPLEKSOWA ANALIZA GENOW PREDYSPOZYCJI DO RAKA
# BRCA1/2, TP53, Lynch syndrome, i wiele innych
#===============================================================================

VCF="${1:-saryd_variants.vcf.gz}"
OUTDIR="cancer_genes_analysis"
mkdir -p "$OUTDIR"

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║              KOMPLEKSOWA ANALIZA GENOW PREDYSPOZYCJI DO RAKA              ║"
echo "║                   BRCA1/2, TP53, Lynch, i inne                            ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

check_snp() {
    local pos="$1"
    bcftools query -r "$pos" -f '[%GT]\n' "$VCF" 2>/dev/null | head -1
}

check_region() {
    local region="$1"
    bcftools view -r "$region" -H "$VCF" 2>/dev/null | wc -l
}

get_variants() {
    local region="$1"
    bcftools query -r "$region" -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t[%GT]\n' "$VCF" 2>/dev/null
}

#===============================================================================
# BRCA1 i BRCA2
#===============================================================================

echo "🎀 BRCA1 i BRCA2 - Rak piersi i jajnika"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

BRCA1_VARIANTS=$(check_region "17:43044295-43125483")
BRCA2_VARIANTS=$(check_region "13:32315474-32400266")

echo "   BRCA1: $BRCA1_VARIANTS wariantow"
echo "   BRCA2: $BRCA2_VARIANTS wariantow"
echo ""

BRCA1_FOUND=0
BRCA2_FOUND=0

# BRCA1 znane mutacje
for variant in "17:43045677|185delAG|Ashkenazi" "17:43057051|5382insC|Slavic" "17:43071077|M1775R|BRCT"; do
    pos=$(echo "$variant" | cut -d'|' -f1)
    name=$(echo "$variant" | cut -d'|' -f2)
    pop=$(echo "$variant" | cut -d'|' -f3)
    
    GT=$(check_snp "$pos")
    if [ -n "$GT" ] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 BRCA1 $name: $GT ($pop)"
        BRCA1_FOUND=$((BRCA1_FOUND + 1))
    fi
done

# BRCA2 znane mutacje
for variant in "13:32354860|6174delT|Ashkenazi" "13:32336282|4075delGT|Pathogenic" "13:32370557|886delGT|Icelandic"; do
    pos=$(echo "$variant" | cut -d'|' -f1)
    name=$(echo "$variant" | cut -d'|' -f2)
    pop=$(echo "$variant" | cut -d'|' -f3)
    
    GT=$(check_snp "$pos")
    if [ -n "$GT" ] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 BRCA2 $name: $GT ($pop)"
        BRCA2_FOUND=$((BRCA2_FOUND + 1))
    fi
done

if [ $BRCA1_FOUND -eq 0 ] && [ $BRCA2_FOUND -eq 0 ]; then
    echo "  ✅ Nie znaleziono znanych patogennych mutacji BRCA1/2"
fi

# Zapisz wszystkie warianty
get_variants "17:43044295-43125483" > "$OUTDIR/BRCA1_all.txt"
get_variants "13:32315474-32400266" > "$OUTDIR/BRCA2_all.txt"

BRCA1_NOVEL=$(grep -c "\.$" "$OUTDIR/BRCA1_all.txt" 2>/dev/null || echo "0")
BRCA2_NOVEL=$(grep -c "\.$" "$OUTDIR/BRCA2_all.txt" 2>/dev/null || echo "0")

echo ""
echo "   Warianty bez rsID (potencjalnie rzadkie):"
echo "   BRCA1: $BRCA1_NOVEL | BRCA2: $BRCA2_NOVEL"
echo ""

#===============================================================================
# TP53 (Li-Fraumeni)
#===============================================================================

echo "🧬 TP53 - Li-Fraumeni Syndrome"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

TP53_VARIANTS=$(check_region "17:7668402-7687538")
echo "   TP53: $TP53_VARIANTS wariantow"
echo ""

TP53_FOUND=0

for variant in "17:7674220|R175H|Hotspot" "17:7675076|R248Q|Hotspot" "17:7675994|R273H|Hotspot" "17:7676154|R282W|Pathogenic"; do
    pos=$(echo "$variant" | cut -d'|' -f1)
    name=$(echo "$variant" | cut -d'|' -f2)
    type=$(echo "$variant" | cut -d'|' -f3)
    
    GT=$(check_snp "$pos")
    if [ -n "$GT" ] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 TP53 $name: $GT ($type)"
        TP53_FOUND=$((TP53_FOUND + 1))
    fi
done

if [ $TP53_FOUND -eq 0 ]; then
    echo "  ✅ Nie znaleziono patogennych mutacji TP53"
fi
echo ""

#===============================================================================
# Lynch Syndrome
#===============================================================================

echo "🔵 Lynch Syndrome - MMR genes"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

MLH1=$(check_region "3:36993350-37050918")
MSH2=$(check_region "2:47403067-47709917")
MSH6=$(check_region "2:47695530-47810101")
PMS2=$(check_region "7:5970925-6009106")

echo "   MLH1: $MLH1 | MSH2: $MSH2 | MSH6: $MSH6 | PMS2: $PMS2"
echo ""

LYNCH_FOUND=0

for variant in "3:37042337|MLH1_c.350C>T" "2:47641560|MSH2_c.1216C>T" "2:47703303|MSH6_c.3261dup"; do
    pos=$(echo "$variant" | cut -d'|' -f1)
    name=$(echo "$variant" | cut -d'|' -f2)
    
    GT=$(check_snp "$pos")
    if [ -n "$GT" ] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
        echo "  🔴 $name: $GT"
        LYNCH_FOUND=$((LYNCH_FOUND + 1))
    fi
done

if [ $LYNCH_FOUND -eq 0 ]; then
    echo "  ✅ Nie znaleziono patogennych mutacji Lynch"
fi
echo ""

#===============================================================================
# Inne geny rakowe
#===============================================================================

echo "🧬 INNE GENY PREDYSPOZYCJI DO RAKA"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "   GEN       | WARIANTY | NOWOTWORY"
echo "   ─────────────────────────────────────────────────────────"

for gene_info in "PTEN|10:87863625-87971930|Breast_thyroid" "APC|5:112707498-112846239|Colorectal" \
                  "CDH1|16:68771117-68869444|Gastric" "CHEK2|22:28687743-28742422|Breast" \
                  "PALB2|16:23603160-23641310|Breast_pancreatic" "ATM|11:108222484-108369102|Breast"; do
    gene=$(echo "$gene_info" | cut -d'|' -f1)
    region=$(echo "$gene_info" | cut -d'|' -f2)
    cancers=$(echo "$gene_info" | cut -d'|' -f3)
    
    count=$(check_region "$region")
    printf "   %-10s | %8d | %s\n" "$gene" "$count" "${cancers//_/ }"
done

echo ""

# CHEK2 1100delC - specjalny wariant
GT=$(check_snp "22:28695868")
if [ -n "$GT" ] && [[ "$GT" != "0/0" ]] && [[ "$GT" != "./." ]]; then
    echo "  ⚠️  CHEK2 1100delC: $GT (2x ryzyko raka piersi)"
fi

# MUTYH - recesywny
GT1=$(check_snp "1:45332879")
GT2=$(check_snp "1:45340175")
MUTYH_COUNT=0
[[ "$GT1" == *"1"* ]] && MUTYH_COUNT=$((MUTYH_COUNT + 1))
[[ "$GT2" == *"1"* ]] && MUTYH_COUNT=$((MUTYH_COUNT + 1))

if [ $MUTYH_COUNT -gt 0 ]; then
    echo "  ℹ️  MUTYH: $MUTYH_COUNT wariant(y) (recesywny - potrzebne 2 dla choroby)"
fi

echo ""

#===============================================================================
# PODSUMOWANIE
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                              PODSUMOWANIE                                 ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

TOTAL=$((BRCA1_FOUND + BRCA2_FOUND + TP53_FOUND + LYNCH_FOUND))

echo "   BRCA1 patogenne:  $BRCA1_FOUND"
echo "   BRCA2 patogenne:  $BRCA2_FOUND"
echo "   TP53 patogenne:   $TP53_FOUND"
echo "   Lynch patogenne:  $LYNCH_FOUND"
echo "   ────────────────────"
echo "   RAZEM:            $TOTAL"
echo ""

if [ $TOTAL -eq 0 ]; then
    echo "✅ NIE ZNALEZIONO znanych patogennych mutacji w genach rakowych!"
    echo ""
    echo "   Pamietaj: ~90% rakow jest sporadycznych, nie dziedzicznych."
    echo "   Regularne badania przesiewowe wazne niezaleznie od genetyki."
else
    echo "⚠️  ZNALEZIONO $TOTAL patogennych wariantow!"
    echo "   → Zalecana konsultacja z genetykiem klinicznym"
fi

echo ""
echo "📁 Pliki: $OUTDIR/BRCA1_all.txt, BRCA2_all.txt"
echo ""
