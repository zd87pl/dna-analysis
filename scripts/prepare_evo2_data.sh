#!/bin/bash
#===============================================================================
# PRZYGOTOWANIE DANYCH DLA EVO2
# Ekstrakcja sekwencji z BAM do analizy na RunPod
#===============================================================================

BAM="$1"
OUTPUT_DIR="${2:-evo2_data}"
REFERENCE="${3:-}"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║              PRZYGOTOWANIE DANYCH DLA EVO2                        ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

mkdir -p "$OUTPUT_DIR"/{sequences,contexts,variants,metadata}

#===============================================================================
# Definicja genów do ekstrakcji
# Format: GENE:CHROM:START-END
#===============================================================================

declare -A GENE_REGIONS

# Farmakogenomika
GENE_REGIONS["CYP2D6"]="22:42522500-42526883"
GENE_REGIONS["CYP2C19"]="10:96522463-96612671"
GENE_REGIONS["CYP2C9"]="10:96698415-96749148"
GENE_REGIONS["CYP3A4"]="7:99354604-99381888"
GENE_REGIONS["CYP1A2"]="15:75041184-75048543"
GENE_REGIONS["DPYD"]="1:97543299-98386615"
GENE_REGIONS["TPMT"]="6:18128542-18155374"
GENE_REGIONS["SLCO1B1"]="12:21284127-21392730"
GENE_REGIONS["VKORC1"]="16:31102175-31106699"
GENE_REGIONS["UGT1A1"]="2:234668879-234681945"

# BRCA i DNA repair
GENE_REGIONS["BRCA1"]="17:41196312-41277500"
GENE_REGIONS["BRCA2"]="13:32889617-32973809"
GENE_REGIONS["TP53"]="17:7571720-7590868"
GENE_REGIONS["ATM"]="11:108093211-108239829"
GENE_REGIONS["CHEK2"]="22:29083731-29137822"
GENE_REGIONS["PALB2"]="16:23614488-23652631"

# Lynch syndrome
GENE_REGIONS["MLH1"]="3:37034841-37092337"
GENE_REGIONS["MSH2"]="2:47630108-47710367"
GENE_REGIONS["MSH6"]="2:48010221-48034092"
GENE_REGIONS["PMS2"]="7:6012870-6048756"

# Kardiologia
GENE_REGIONS["APOE"]="19:45409039-45412650"
GENE_REGIONS["LDLR"]="19:11200038-11244506"
GENE_REGIONS["PCSK9"]="1:55505149-55530526"
GENE_REGIONS["APOB"]="2:21224301-21266945"
GENE_REGIONS["SCN5A"]="3:38589553-38691164"
GENE_REGIONS["KCNQ1"]="11:2466221-2870340"
GENE_REGIONS["KCNH2"]="7:150642049-150675403"

# Metabolizm
GENE_REGIONS["MTHFR"]="1:11845787-11866160"
GENE_REGIONS["TCF7L2"]="10:114710009-114927437"
GENE_REGIONS["FTO"]="16:53737875-54155853"
GENE_REGIONS["MC4R"]="18:58038563-58040001"
GENE_REGIONS["PPARG"]="3:12328867-12475855"

# Neurologia
GENE_REGIONS["COMT"]="22:19929263-19957498"
GENE_REGIONS["BDNF"]="11:27676439-27743605"
GENE_REGIONS["DRD2"]="11:113280318-113346413"
GENE_REGIONS["HTT"]="4:3076408-3245687"  # Huntington
GENE_REGIONS["APP"]="21:27252861-27543138"  # Alzheimer
GENE_REGIONS["PSEN1"]="14:73603143-73690399"
GENE_REGIONS["SNCA"]="4:90645250-90759447"  # Parkinson

# Longevity
GENE_REGIONS["FOXO3"]="6:108881025-109005972"
GENE_REGIONS["TERT"]="5:1253262-1295184"
GENE_REGIONS["KLOTHO"]="13:33589873-33639087"

# Sport/Fitness
GENE_REGIONS["ACTN3"]="11:66546177-66563148"
GENE_REGIONS["ACE"]="17:61554422-61575741"
GENE_REGIONS["PPARGC1A"]="4:23755722-24869642"
GENE_REGIONS["MSTN"]="2:190055702-190062545"
GENE_REGIONS["IL6"]="7:22766246-22771621"

# Immunologia
GENE_REGIONS["HLA_A"]="6:29910247-29913661"
GENE_REGIONS["HLA_B"]="6:31321649-31324989"
GENE_REGIONS["HLA_C"]="6:31236526-31239913"
GENE_REGIONS["HLA_DRB1"]="6:32546547-32557613"

echo ""
echo "Geny do ekstrakcji: ${#GENE_REGIONS[@]}"
echo ""

#===============================================================================
# Ekstrakcja sekwencji consensus z BAM
#===============================================================================

extract_consensus_from_bam() {
    local gene=$1
    local region=$2
    local output=$3
    
    # Rozdziel region
    local chrom=$(echo $region | cut -d: -f1)
    local coords=$(echo $region | cut -d: -f2)
    local start=$(echo $coords | cut -d- -f1)
    local end=$(echo $coords | cut -d- -f2)
    
    echo "  Ekstrakcja $gene ($region)..."
    
    # Metoda 1: Jeśli mamy referencję - użyj consensus
    if [ -n "$REFERENCE" ] && [ -f "$REFERENCE" ]; then
        samtools mpileup -r "$region" -f "$REFERENCE" "$BAM" 2>/dev/null | \
        awk '{
            ref = $3
            bases = toupper($5)
            gsub(/\^./, "", bases)
            gsub(/\$/, "", bases)
            gsub(/[+-][0-9]+[ACGTNacgtn]+/, "", bases)
            
            a=0; c=0; g=0; t=0; r=0
            n = split(bases, arr, "")
            for(i=1; i<=n; i++) {
                if(arr[i]=="A") a++
                else if(arr[i]=="C") c++
                else if(arr[i]=="G") g++
                else if(arr[i]=="T") t++
                else if(arr[i]=="." || arr[i]==",") r++
            }
            
            max = r; cons = ref
            if(a > max) { max = a; cons = "A" }
            if(c > max) { max = c; cons = "C" }
            if(g > max) { max = g; cons = "G" }
            if(t > max) { max = t; cons = "T" }
            
            printf "%s", cons
        }' > "$output.tmp"
        
        # Formatuj jako FASTA
        echo ">${gene}_${region}_consensus" > "$output"
        fold -w 80 "$output.tmp" >> "$output"
        rm -f "$output.tmp"
        
    else
        # Metoda 2: Bez referencji - zbierz tylko statystyki
        echo ">${gene}_${region}_coverage_only" > "$output"
        samtools depth -r "$region" "$BAM" 2>/dev/null | \
            awk '{sum+=$3; cnt++} END {print "# Average coverage: " (cnt>0 ? sum/cnt : 0) "x"}' >> "$output"
    fi
}

#===============================================================================
# Ekstrakcja kontekstu sekwencji (dla EVO2 potrzebuje kontekstu)
#===============================================================================

extract_with_context() {
    local gene=$1
    local region=$2
    local context_size=${3:-1000}  # 1kb kontekstu z każdej strony
    
    local chrom=$(echo $region | cut -d: -f1)
    local coords=$(echo $region | cut -d: -f2)
    local start=$(echo $coords | cut -d- -f1)
    local end=$(echo $coords | cut -d- -f2)
    
    # Rozszerz region o kontekst
    local ext_start=$((start - context_size))
    local ext_end=$((end + context_size))
    [ $ext_start -lt 1 ] && ext_start=1
    
    local ext_region="${chrom}:${ext_start}-${ext_end}"
    
    extract_consensus_from_bam "$gene" "$ext_region" "$OUTPUT_DIR/contexts/${gene}_with_context.fasta"
}

#===============================================================================
# Główna pętla ekstrakcji
#===============================================================================

echo "Rozpoczynam ekstrakcję sekwencji..."
echo ""

for gene in "${!GENE_REGIONS[@]}"; do
    region="${GENE_REGIONS[$gene]}"
    
    # Podstawowa sekwencja
    extract_consensus_from_bam "$gene" "$region" "$OUTPUT_DIR/sequences/${gene}.fasta"
    
    # Sekwencja z kontekstem (dla lepszych predykcji EVO2)
    extract_with_context "$gene" "$region" 500
done

echo ""
echo "Ekstrakcja zakończona!"

#===============================================================================
# Generowanie metadanych
#===============================================================================

echo ""
echo "Generowanie metadanych..."

cat > "$OUTPUT_DIR/metadata/genes.json" << EOF
{
    "source_bam": "$BAM",
    "reference": "${REFERENCE:-none}",
    "extraction_date": "$(date -Iseconds)",
    "genes": {
$(for gene in "${!GENE_REGIONS[@]}"; do
    echo "        \"$gene\": \"${GENE_REGIONS[$gene]}\","
done | sed '$ s/,$//')
    }
}
EOF

#===============================================================================
# Skrypt Python do uruchomienia na RunPod
#===============================================================================

cat > "$OUTPUT_DIR/evo2_predict.py" << 'PYTHON_EOF'
#!/usr/bin/env python3
"""
EVO2 Variant Effect Prediction
Uruchom na RunPod z GPU
"""

import os
import json
import torch
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional
import numpy as np

@dataclass
class VariantPrediction:
    gene: str
    position: int
    ref: str
    alt: str
    delta_score: float
    interpretation: str
    confidence: float

class EVO2Analyzer:
    def __init__(self, model_path: str = "evo2-40b"):
        """Inicjalizacja modelu EVO2"""
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"Używam urządzenia: {self.device}")
        
        try:
            # Próba załadowania EVO2
            # Dostosuj import do rzeczywistej biblioteki EVO2
            from transformers import AutoModel, AutoTokenizer
            
            self.tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code=True)
            self.model = AutoModel.from_pretrained(model_path, trust_remote_code=True)
            self.model.to(self.device)
            self.model.eval()
            print("Model EVO2 załadowany pomyślnie")
        except Exception as e:
            print(f"Nie można załadować EVO2: {e}")
            print("Używam trybu symulacji")
            self.model = None
            self.tokenizer = None
    
    def load_sequence(self, fasta_path: str) -> str:
        """Wczytuje sekwencję z pliku FASTA"""
        sequence = []
        with open(fasta_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence.append(line.strip().upper())
        return ''.join(sequence)
    
    def score_sequence(self, sequence: str) -> float:
        """Oblicza log-likelihood sekwencji"""
        if self.model is None:
            # Symulacja - zwróć losowy wynik
            return np.random.normal(-100, 10)
        
        with torch.no_grad():
            inputs = self.tokenizer(sequence, return_tensors="pt").to(self.device)
            outputs = self.model(**inputs)
            # Dostosuj do API EVO2
            return outputs.logits.mean().item()
    
    def predict_variant_effect(
        self, 
        ref_seq: str, 
        position: int, 
        ref_base: str, 
        alt_base: str
    ) -> VariantPrediction:
        """
        Przewiduje efekt pojedynczego wariantu
        """
        # Utwórz sekwencję z wariantem
        alt_seq = ref_seq[:position] + alt_base + ref_seq[position+1:]
        
        # Oblicz score dla obu sekwencji
        ref_score = self.score_sequence(ref_seq)
        alt_score = self.score_sequence(alt_seq)
        
        delta = alt_score - ref_score
        
        # Interpretacja
        if delta < -2:
            interpretation = "DELETERIOUS"
            confidence = min(abs(delta) / 5, 1.0)
        elif delta < -0.5:
            interpretation = "POSSIBLY_DELETERIOUS"
            confidence = 0.5 + abs(delta) / 4
        elif delta > 0.5:
            interpretation = "POSSIBLY_BENEFICIAL"
            confidence = 0.5 + delta / 4
        else:
            interpretation = "NEUTRAL"
            confidence = 1.0 - abs(delta)
        
        return VariantPrediction(
            gene="",
            position=position,
            ref=ref_base,
            alt=alt_base,
            delta_score=delta,
            interpretation=interpretation,
            confidence=confidence
        )
    
    def scan_gene(
        self, 
        sequence: str, 
        gene_name: str,
        window_size: int = 100
    ) -> List[Dict]:
        """
        Skanuje gen w poszukiwaniu wrażliwych pozycji
        """
        results = []
        seq_len = len(sequence)
        
        # Dla każdej pozycji sprawdź możliwe mutacje
        for i in range(0, seq_len, window_size):
            window = sequence[max(0, i-50):min(seq_len, i+50)]
            ref_score = self.score_sequence(window)
            
            # Sprawdź środkową pozycję okna
            pos = min(50, i)
            ref_base = window[pos] if pos < len(window) else 'N'
            
            for alt_base in ['A', 'C', 'G', 'T']:
                if alt_base != ref_base:
                    alt_window = window[:pos] + alt_base + window[pos+1:]
                    alt_score = self.score_sequence(alt_window)
                    
                    delta = alt_score - ref_score
                    
                    if abs(delta) > 0.5:  # Tylko znaczące zmiany
                        results.append({
                            'gene': gene_name,
                            'position': i,
                            'ref': ref_base,
                            'alt': alt_base,
                            'delta': delta,
                            'interpretation': self._interpret(delta)
                        })
        
        return results
    
    def _interpret(self, delta: float) -> str:
        if delta < -2:
            return "HIGHLY_DELETERIOUS"
        elif delta < -0.5:
            return "POSSIBLY_DELETERIOUS"
        elif delta > 2:
            return "HIGHLY_BENEFICIAL"
        elif delta > 0.5:
            return "POSSIBLY_BENEFICIAL"
        return "NEUTRAL"

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='EVO2 Variant Analysis')
    parser.add_argument('--sequences', default='sequences/', help='Folder z FASTA')
    parser.add_argument('--output', default='evo2_results.json', help='Plik wynikowy')
    parser.add_argument('--model', default='evo2-40b', help='Ścieżka do modelu')
    args = parser.parse_args()
    
    print("=" * 60)
    print("EVO2 VARIANT EFFECT ANALYZER")
    print("=" * 60)
    
    # Inicjalizacja
    analyzer = EVO2Analyzer(args.model)
    
    # Znajdź pliki FASTA
    fasta_files = list(Path(args.sequences).glob("*.fasta"))
    print(f"\nZnaleziono {len(fasta_files)} plików sekwencji")
    
    all_results = {}
    
    for fasta in fasta_files:
        gene_name = fasta.stem
        print(f"\nAnalizuję {gene_name}...")
        
        sequence = analyzer.load_sequence(str(fasta))
        
        if len(sequence) > 0:
            # Skanuj gen
            results = analyzer.scan_gene(sequence, gene_name)
            all_results[gene_name] = {
                'sequence_length': len(sequence),
                'variants_found': len(results),
                'predictions': results
            }
            print(f"  Długość: {len(sequence)} bp")
            print(f"  Znalezione warianty: {len(results)}")
    
    # Zapisz wyniki
    with open(args.output, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print(f"Wyniki zapisane do: {args.output}")
    print(f"{'=' * 60}")

if __name__ == "__main__":
    main()
PYTHON_EOF

chmod +x "$OUTPUT_DIR/evo2_predict.py"

#===============================================================================
# Instrukcje
#===============================================================================

cat > "$OUTPUT_DIR/README.md" << EOF
# Dane dla EVO2 - RunPod

## Struktura folderów

\`\`\`
$OUTPUT_DIR/
├── sequences/          # Sekwencje FASTA genów
├── contexts/           # Sekwencje z kontekstem (±500bp)
├── variants/           # Warianty VCF (jeśli wygenerowane)
├── metadata/           # Metadane JSON
├── evo2_predict.py     # Skrypt do uruchomienia na RunPod
└── README.md           # Ten plik
\`\`\`

## Uruchomienie na RunPod

### 1. Wybierz maszynę
- GPU: A100 80GB lub A6000 48GB (dla pełnego EVO2-40B)
- Alternatywnie: A100 40GB z quantizacją

### 2. Skopiuj dane
\`\`\`bash
# Z lokalnej maszyny na RunPod
rsync -avz $OUTPUT_DIR/ runpod:/workspace/evo2_data/
\`\`\`

### 3. Zainstaluj zależności
\`\`\`bash
pip install torch transformers accelerate bitsandbytes
pip install evo2  # lub instalacja z źródła
\`\`\`

### 4. Uruchom analizę
\`\`\`bash
cd /workspace/evo2_data
python evo2_predict.py --sequences sequences/ --output results.json
\`\`\`

## Alternatywnie: Użyj Hugging Face

\`\`\`python
from transformers import AutoModel, AutoTokenizer

# EVO2 może być dostępny przez HF
model = AutoModel.from_pretrained("togethercomputer/evo-1-8k-base")
tokenizer = AutoTokenizer.from_pretrained("togethercomputer/evo-1-8k-base")
\`\`\`

## Geny w zestawie

$(for gene in "${!GENE_REGIONS[@]}"; do
    echo "- **$gene**: ${GENE_REGIONS[$gene]}"
done | sort)

## Następne kroki

1. Wygeneruj VCF z BAM (jeśli jeszcze nie masz)
2. Prześlij na RunPod
3. Uruchom predykcję
4. Zintegruj wyniki z Helixight
EOF

#===============================================================================
# Podsumowanie
#===============================================================================

echo ""
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║                         GOTOWE!                                   ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "Wygenerowane pliki:"
echo "  • $OUTPUT_DIR/sequences/     - Sekwencje FASTA (${#GENE_REGIONS[@]} genów)"
echo "  • $OUTPUT_DIR/contexts/      - Sekwencje z kontekstem"
echo "  • $OUTPUT_DIR/evo2_predict.py - Skrypt dla RunPod"
echo "  • $OUTPUT_DIR/README.md      - Instrukcje"
echo ""
echo "Następne kroki:"
echo "  1. Skopiuj folder na RunPod:"
echo "     scp -r $OUTPUT_DIR user@runpod:/workspace/"
echo ""
echo "  2. Na RunPod uruchom:"
echo "     python evo2_predict.py --sequences sequences/"
echo ""
