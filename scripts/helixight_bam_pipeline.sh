#!/bin/bash
#===============================================================================
# HELIXIGHT BAM PROCESSING PIPELINE
# Generuje VCF, FASTA i dane dla EVO2
#===============================================================================

set -e

# Konfiguracja
BAM_FILE="${1:-saryd.bam}"
OUTPUT_DIR="${2:-./helixight_output}"
THREADS="${3:-8}"
REFERENCE="${4:-/path/to/GRCh37.fa}"  # Musisz pobrać!

# Kolory
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}"
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║           HELIXIGHT BAM PROCESSING PIPELINE                       ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Tworzenie katalogów
mkdir -p "$OUTPUT_DIR"/{vcf,fasta,evo2,qc,logs}

#===============================================================================
# KROK 0: Sprawdzenie narzędzi
#===============================================================================
echo -e "${YELLOW}[0/7] Sprawdzanie wymaganych narzędzi...${NC}"

check_tool() {
    if command -v $1 &> /dev/null; then
        echo -e "  ${GREEN}✓${NC} $1 znaleziony"
        return 0
    else
        echo -e "  ${RED}✗${NC} $1 NIE ZNALEZIONY"
        return 1
    fi
}

MISSING_TOOLS=0
check_tool samtools || MISSING_TOOLS=1
check_tool bcftools || MISSING_TOOLS=1
check_tool bgzip || MISSING_TOOLS=1
check_tool tabix || MISSING_TOOLS=1

# Opcjonalne ale zalecane
check_tool gatk || echo -e "  ${YELLOW}!${NC} GATK opcjonalny (lepsze wyniki)"
check_tool freebayes || echo -e "  ${YELLOW}!${NC} FreeBayes opcjonalny"

if [ $MISSING_TOOLS -eq 1 ]; then
    echo -e "${RED}Brakujące narzędzia! Zainstaluj przez:${NC}"
    echo "  brew install samtools bcftools htslib  # macOS"
    echo "  apt install samtools bcftools tabix    # Linux"
    exit 1
fi

#===============================================================================
# KROK 1: QC i statystyki BAM
#===============================================================================
echo -e "\n${YELLOW}[1/7] Analiza jakości BAM...${NC}"

echo "Generowanie statystyk flagstat..."
samtools flagstat "$BAM_FILE" > "$OUTPUT_DIR/qc/flagstat.txt" &

echo "Generowanie statystyk idxstats..."
samtools idxstats "$BAM_FILE" > "$OUTPUT_DIR/qc/idxstats.txt" &

echo "Obliczanie pokrycia (może potrwać)..."
samtools depth -a "$BAM_FILE" | \
    awk '{sum+=$3; cnt++} END {print "Średnie pokrycie: " sum/cnt "x"}' \
    > "$OUTPUT_DIR/qc/coverage_summary.txt" &

wait
echo -e "${GREEN}✓ QC zakończone${NC}"

#===============================================================================
# KROK 2: Variant Calling - METODA SZYBKA (bcftools)
#===============================================================================
echo -e "\n${YELLOW}[2/7] Variant calling (bcftools mpileup)...${NC}"
echo "To zajmie 2-4 godziny dla WGS..."

# Sprawdź czy mamy referencję
if [ ! -f "$REFERENCE" ]; then
    echo -e "${YELLOW}UWAGA: Brak pliku referencyjnego!${NC}"
    echo "Pobierz GRCh37/hg19:"
    echo "  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
    echo "  gunzip human_g1k_v37.fasta.gz"
    echo "  samtools faidx human_g1k_v37.fasta"
    echo ""
    echo "Lub GRCh38/hg38:"
    echo "  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    
    # Spróbuj użyć bez referencji (ograniczone możliwości)
    echo -e "\n${YELLOW}Próbuję kontynuować bez referencji (ograniczone)...${NC}"
    SKIP_REFERENCE=1
else
    SKIP_REFERENCE=0
fi

if [ $SKIP_REFERENCE -eq 0 ]; then
    # Pełny variant calling z referencją
    bcftools mpileup \
        -Ou \
        -f "$REFERENCE" \
        --threads "$THREADS" \
        --max-depth 250 \
        --min-MQ 20 \
        --min-BQ 20 \
        "$BAM_FILE" | \
    bcftools call \
        -mv \
        --threads "$THREADS" \
        -Oz \
        -o "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz"
    
    # Indeksowanie
    tabix -p vcf "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz"
    
    echo -e "${GREEN}✓ VCF wygenerowany: $OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz${NC}"
fi

#===============================================================================
# KROK 3: Filtrowanie i podział VCF
#===============================================================================
echo -e "\n${YELLOW}[3/7] Filtrowanie wariantów...${NC}"

if [ -f "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" ]; then
    # Tylko SNPy wysokiej jakości
    bcftools view \
        -v snps \
        -i 'QUAL>20 && DP>10' \
        "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" \
        -Oz -o "$OUTPUT_DIR/vcf/saryd_snps_filtered.vcf.gz"
    
    # Tylko indele
    bcftools view \
        -v indels \
        -i 'QUAL>20 && DP>5' \
        "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" \
        -Oz -o "$OUTPUT_DIR/vcf/saryd_indels_filtered.vcf.gz"
    
    # Indeksowanie
    tabix -p vcf "$OUTPUT_DIR/vcf/saryd_snps_filtered.vcf.gz"
    tabix -p vcf "$OUTPUT_DIR/vcf/saryd_indels_filtered.vcf.gz"
    
    # Statystyki
    echo "Statystyki wariantów:" > "$OUTPUT_DIR/vcf/variant_stats.txt"
    bcftools stats "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" >> "$OUTPUT_DIR/vcf/variant_stats.txt"
    
    echo -e "${GREEN}✓ Filtrowanie zakończone${NC}"
fi

#===============================================================================
# KROK 4: Ekstrakcja sekwencji genów dla EVO2
#===============================================================================
echo -e "\n${YELLOW}[4/7] Przygotowanie danych dla EVO2...${NC}"

# Lista kluczowych genów do analizy EVO2
GENES_OF_INTEREST=(
    # Farmakogenomika
    "CYP2D6:22:42522500-42526883"
    "CYP2C19:10:96522463-96612671"
    "CYP2C9:10:96698415-96749148"
    "CYP3A4:7:99354604-99381888"
    "DPYD:1:97543299-98386615"
    "TPMT:6:18128542-18155374"
    "SLCO1B1:12:21284127-21392730"
    
    # Nowotwory
    "BRCA1:17:41196312-41277500"
    "BRCA2:13:32889617-32973809"
    "TP53:17:7571720-7590868"
    "MLH1:3:37034841-37092337"
    "MSH2:2:47630108-47710367"
    
    # Kardiologia
    "APOE:19:45409039-45412650"
    "LDLR:19:11200038-11244506"
    "PCSK9:1:55505149-55530526"
    
    # Metabolizm
    "MTHFR:1:11845787-11866160"
    "TCF7L2:10:114710009-114927437"
    "FTO:16:53737875-54155853"
    
    # Neurologia
    "APOE:19:45409039-45412650"
    "COMT:22:19929263-19957498"
    "BDNF:11:27676439-27743605"
    "DRD2:11:113280318-113346413"
    
    # Longevity
    "FOXO3:6:108881025-109005972"
    "TERT:5:1253262-1295184"
    
    # Sport
    "ACTN3:11:66546177-66563148"
    "ACE:17:61554422-61575741"
    "PPARGC1A:4:23755722-24869642"
)

mkdir -p "$OUTPUT_DIR/evo2/gene_sequences"
mkdir -p "$OUTPUT_DIR/evo2/gene_variants"

echo "Ekstrakcja sekwencji genów..."

for gene_info in "${GENES_OF_INTEREST[@]}"; do
    GENE_NAME=$(echo $gene_info | cut -d: -f1)
    CHROM=$(echo $gene_info | cut -d: -f2)
    REGION=$(echo $gene_info | cut -d: -f3)
    
    echo "  Przetwarzanie $GENE_NAME ($CHROM:$REGION)..."
    
    # Ekstrakcja consensus sequence z BAM
    if [ $SKIP_REFERENCE -eq 0 ]; then
        samtools mpileup -r "$CHROM:$REGION" -f "$REFERENCE" "$BAM_FILE" 2>/dev/null | \
        perl -lane '
            $cov = $F[3];
            $bases = uc($F[4]);
            $bases =~ s/[\^\$]//g;
            $bases =~ s/[+-]\d+[ACGTNacgtn]+//g;
            
            %counts = (A=>0, C=>0, G=>0, T=>0);
            for $b (split //, $bases) {
                $b = uc($b);
                $counts{$b}++ if exists $counts{$b};
                $counts{$F[2]}++ if $b eq "." || $b eq ",";
            }
            
            @sorted = sort {$counts{$b} <=> $counts{$a}} keys %counts;
            $consensus = $sorted[0];
            $consensus = $F[2] if $counts{$sorted[0]} == 0;
            print $consensus;
        ' | tr -d '\n' > "$OUTPUT_DIR/evo2/gene_sequences/${GENE_NAME}_consensus.txt"
        
        # Dodaj nagłówek FASTA
        echo ">${GENE_NAME}_${CHROM}:${REGION}" > "$OUTPUT_DIR/evo2/gene_sequences/${GENE_NAME}.fasta"
        fold -w 80 "$OUTPUT_DIR/evo2/gene_sequences/${GENE_NAME}_consensus.txt" >> "$OUTPUT_DIR/evo2/gene_sequences/${GENE_NAME}.fasta"
        rm "$OUTPUT_DIR/evo2/gene_sequences/${GENE_NAME}_consensus.txt"
    fi
    
    # Ekstrakcja wariantów z VCF dla tego genu
    if [ -f "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" ]; then
        bcftools view -r "$CHROM:$REGION" "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" \
            > "$OUTPUT_DIR/evo2/gene_variants/${GENE_NAME}_variants.vcf" 2>/dev/null || true
    fi
done

echo -e "${GREEN}✓ Dane EVO2 przygotowane${NC}"

#===============================================================================
# KROK 5: Generowanie osobistego genomu FASTA
#===============================================================================
echo -e "\n${YELLOW}[5/7] Generowanie spersonalizowanego genomu...${NC}"

if [ $SKIP_REFERENCE -eq 0 ] && [ -f "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" ]; then
    echo "Tworzenie consensus FASTA z Twoimi wariantami..."
    
    # Dla każdego chromosomu
    for CHR in {1..22} X Y MT; do
        echo "  Chromosom $CHR..."
        bcftools consensus \
            -f "$REFERENCE" \
            -s - \
            "$OUTPUT_DIR/vcf/saryd_all_variants.vcf.gz" \
            -r "$CHR" \
            >> "$OUTPUT_DIR/fasta/saryd_personal_genome.fasta" 2>/dev/null || true
    done
    
    echo -e "${GREEN}✓ Osobisty genom: $OUTPUT_DIR/fasta/saryd_personal_genome.fasta${NC}"
else
    echo -e "${YELLOW}Pominięto (wymaga referencji i VCF)${NC}"
fi

#===============================================================================
# KROK 6: Przygotowanie formatu dla EVO2 RunPod
#===============================================================================
echo -e "\n${YELLOW}[6/7] Formatowanie dla EVO2 (RunPod)...${NC}"

# Tworzenie JSON z metadanymi
cat > "$OUTPUT_DIR/evo2/sample_metadata.json" << EOF
{
    "sample_id": "saryd",
    "reference_genome": "GRCh37",
    "sequencing_type": "WGS",
    "analysis_date": "$(date -Iseconds)",
    "files": {
        "vcf": "saryd_all_variants.vcf.gz",
        "gene_sequences": "gene_sequences/",
        "gene_variants": "gene_variants/"
    },
    "genes_extracted": [
        $(printf '"%s",' "${GENES_OF_INTEREST[@]}" | sed 's/,$//')
    ]
}
EOF

# Skrypt do uruchomienia na RunPod
cat > "$OUTPUT_DIR/evo2/runpod_evo2_analysis.py" << 'PYTHON_SCRIPT'
#!/usr/bin/env python3
"""
EVO2 Analysis Script for RunPod
Analizuje sekwencje genów i przewiduje efekty wariantów
"""

import os
import json
import torch
from pathlib import Path

# Konfiguracja EVO2 (dostosuj do swojej instalacji)
EVO2_MODEL_PATH = "/workspace/evo2-40b"  # Ścieżka na RunPod

def load_evo2_model():
    """Ładuje model EVO2"""
    try:
        from evo2 import Evo2Model
        model = Evo2Model.from_pretrained(EVO2_MODEL_PATH)
        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
        return model
    except ImportError:
        print("EVO2 nie zainstalowany. Instaluj przez:")
        print("  pip install evo2")
        return None

def load_gene_sequence(gene_fasta_path):
    """Wczytuje sekwencję z pliku FASTA"""
    sequence = ""
    with open(gene_fasta_path, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def predict_variant_effect(model, ref_seq, alt_seq, gene_name):
    """
    Przewiduje efekt wariantu używając EVO2
    Porównuje prawdopodobieństwo sekwencji referencyjnej vs alternatywnej
    """
    if model is None:
        return {"error": "Model not loaded"}
    
    with torch.no_grad():
        # Tokenizacja
        ref_tokens = model.tokenize(ref_seq)
        alt_tokens = model.tokenize(alt_seq)
        
        if torch.cuda.is_available():
            ref_tokens = ref_tokens.cuda()
            alt_tokens = alt_tokens.cuda()
        
        # Likelihood sekwencji
        ref_logprob = model.score(ref_tokens)
        alt_logprob = model.score(alt_tokens)
        
        # Delta log-likelihood
        delta_ll = alt_logprob - ref_logprob
        
        return {
            "gene": gene_name,
            "ref_logprob": float(ref_logprob),
            "alt_logprob": float(alt_logprob),
            "delta_logprob": float(delta_ll),
            "interpretation": interpret_delta(delta_ll)
        }

def interpret_delta(delta):
    """Interpretuje różnicę log-likelihood"""
    if delta < -5:
        return "HIGHLY_DELETERIOUS"
    elif delta < -2:
        return "LIKELY_DELETERIOUS"
    elif delta < -0.5:
        return "POSSIBLY_DELETERIOUS"
    elif delta < 0.5:
        return "NEUTRAL"
    elif delta < 2:
        return "POSSIBLY_BENEFICIAL"
    else:
        return "LIKELY_BENEFICIAL"

def analyze_gene_variants(model, gene_dir, variants_dir):
    """Analizuje wszystkie warianty dla genów"""
    results = []
    
    gene_files = list(Path(gene_dir).glob("*.fasta"))
    
    for gene_fasta in gene_files:
        gene_name = gene_fasta.stem
        print(f"Analizuję {gene_name}...")
        
        # Wczytaj sekwencję
        sequence = load_gene_sequence(gene_fasta)
        
        # Wczytaj warianty VCF
        vcf_path = Path(variants_dir) / f"{gene_name}_variants.vcf"
        if vcf_path.exists():
            variants = parse_vcf(vcf_path)
            
            for var in variants:
                # Zastosuj wariant do sekwencji
                alt_sequence = apply_variant(sequence, var)
                
                # Przewiduj efekt
                effect = predict_variant_effect(
                    model, 
                    sequence, 
                    alt_sequence, 
                    gene_name
                )
                effect['variant'] = var
                results.append(effect)
    
    return results

def parse_vcf(vcf_path):
    """Parsuje prosty plik VCF"""
    variants = []
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                variants.append({
                    'chrom': parts[0],
                    'pos': int(parts[1]),
                    'ref': parts[3],
                    'alt': parts[4].split(',')[0]
                })
    return variants

def apply_variant(sequence, variant):
    """Aplikuje wariant do sekwencji (uproszczone)"""
    # W rzeczywistości potrzebujesz mapowania pozycji
    # To jest placeholder
    return sequence

def main():
    print("=" * 60)
    print("EVO2 VARIANT EFFECT PREDICTION")
    print("=" * 60)
    
    # Ładowanie modelu
    print("\nŁadowanie modelu EVO2...")
    model = load_evo2_model()
    
    if model is None:
        print("Nie można załadować modelu. Sprawdź instalację.")
        return
    
    # Ścieżki
    gene_dir = "./gene_sequences"
    variants_dir = "./gene_variants"
    
    # Analiza
    print("\nRozpoczynam analizę wariantów...")
    results = analyze_gene_variants(model, gene_dir, variants_dir)
    
    # Zapis wyników
    output_path = "evo2_variant_effects.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nWyniki zapisane do: {output_path}")
    print(f"Przeanalizowano {len(results)} wariantów")

if __name__ == "__main__":
    main()
PYTHON_SCRIPT

chmod +x "$OUTPUT_DIR/evo2/runpod_evo2_analysis.py"

echo -e "${GREEN}✓ Skrypty EVO2 przygotowane${NC}"

#===============================================================================
# KROK 7: Podsumowanie
#===============================================================================
echo -e "\n${YELLOW}[7/7] Generowanie podsumowania...${NC}"

cat > "$OUTPUT_DIR/README.md" << EOF
# Helixight BAM Processing Output

## Wygenerowane pliki

### VCF (Variant Call Format)
- \`vcf/saryd_all_variants.vcf.gz\` - Wszystkie warianty
- \`vcf/saryd_snps_filtered.vcf.gz\` - Tylko SNPy (QUAL>20, DP>10)
- \`vcf/saryd_indels_filtered.vcf.gz\` - Tylko indele

### FASTA
- \`fasta/saryd_personal_genome.fasta\` - Twój spersonalizowany genom

### EVO2 Data
- \`evo2/gene_sequences/\` - Sekwencje FASTA dla kluczowych genów
- \`evo2/gene_variants/\` - Warianty VCF per gen
- \`evo2/runpod_evo2_analysis.py\` - Skrypt do uruchomienia na RunPod
- \`evo2/sample_metadata.json\` - Metadane próbki

### QC
- \`qc/flagstat.txt\` - Statystyki alignmentu
- \`qc/idxstats.txt\` - Statystyki per chromosom
- \`qc/coverage_summary.txt\` - Średnie pokrycie

## Użycie na RunPod

1. Skopiuj folder \`evo2/\` na RunPod
2. Zainstaluj zależności: \`pip install evo2 torch\`
3. Uruchom: \`python runpod_evo2_analysis.py\`

## Statystyki

EOF

# Dodaj statystyki jeśli dostępne
if [ -f "$OUTPUT_DIR/qc/flagstat.txt" ]; then
    echo "### Alignment Statistics" >> "$OUTPUT_DIR/README.md"
    echo '```' >> "$OUTPUT_DIR/README.md"
    cat "$OUTPUT_DIR/qc/flagstat.txt" >> "$OUTPUT_DIR/README.md"
    echo '```' >> "$OUTPUT_DIR/README.md"
fi

echo -e "${GREEN}"
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║                    PRZETWARZANIE ZAKOŃCZONE!                      ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo "Wygenerowane pliki znajdują się w: $OUTPUT_DIR"
echo ""
echo "Następne kroki:"
echo "  1. Sprawdź statystyki: cat $OUTPUT_DIR/qc/flagstat.txt"
echo "  2. Przejrzyj warianty: zcat $OUTPUT_DIR/vcf/saryd_snps_filtered.vcf.gz | head -100"
echo "  3. Skopiuj dane na RunPod: scp -r $OUTPUT_DIR/evo2 user@runpod:~/"
echo ""
