#!/bin/bash
#===============================================================================
# DeepVariant na M4 Mac - Setup i uruchomienie
# DeepVariant używa sieci neuronowej = może użyć Neural Engine / GPU
#===============================================================================

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║           DeepVariant Setup dla Apple Silicon                     ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

#===============================================================================
# OPCJA 1: Docker (najprostsza, ale wolniejsza na Mac)
#===============================================================================

setup_docker() {
    echo ""
    echo "=== OPCJA 1: Docker ==="
    echo ""
    
    # Sprawdź Docker
    if ! command -v docker &> /dev/null; then
        echo "Docker nie znaleziony. Zainstaluj Docker Desktop dla Mac:"
        echo "  https://www.docker.com/products/docker-desktop/"
        return 1
    fi
    
    echo "Pobieranie DeepVariant..."
    docker pull google/deepvariant:1.6.0
    
    echo ""
    echo "Uruchomienie DeepVariant:"
    echo ""
    cat << 'EOF'
# Ustaw zmienne
BAM=/path/to/saryd.bam
REF=/path/to/human_g1k_v37.fasta
OUTPUT_DIR=/path/to/output

# Uruchom DeepVariant
docker run \
  -v "${BAM%/*}":"${BAM%/*}" \
  -v "${REF%/*}":"${REF%/*}" \
  -v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
  google/deepvariant:1.6.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref="${REF}" \
  --reads="${BAM}" \
  --output_vcf="${OUTPUT_DIR}/deepvariant.vcf.gz" \
  --output_gvcf="${OUTPUT_DIR}/deepvariant.g.vcf.gz" \
  --num_shards=10
EOF
}

#===============================================================================
# OPCJA 2: Singularity (lepsza wydajność)
#===============================================================================

setup_singularity() {
    echo ""
    echo "=== OPCJA 2: Singularity ==="
    echo ""
    
    # Na macOS Singularity wymaga VM
    echo "Singularity na macOS wymaga Lima lub podobnego VM."
    echo "Dla najlepszej wydajności użyj Linux VM lub cloud."
}

#===============================================================================
# OPCJA 3: Native Python (eksperymentalne na M4)
#===============================================================================

setup_native() {
    echo ""
    echo "=== OPCJA 3: Native Python (TensorFlow Metal) ==="
    echo ""
    
    echo "Ta opcja może wykorzystać GPU M4 przez TensorFlow Metal."
    echo ""
    
    # Sprawdź Python
    if ! command -v python3 &> /dev/null; then
        echo "Python3 nie znaleziony. Zainstaluj przez:"
        echo "  brew install python@3.11"
        return 1
    fi
    
    echo "Tworzenie środowiska..."
    python3 -m venv deepvariant_env
    source deepvariant_env/bin/activate
    
    echo "Instalacja TensorFlow Metal (GPU M4)..."
    pip install tensorflow-macos tensorflow-metal
    
    echo ""
    echo "⚠️  UWAGA: Pełny DeepVariant wymaga dodatkowych kroków."
    echo "Oficjalnie wspierany jest tylko Linux."
    echo ""
    echo "Alternatywa: Użyj Google Cloud z GPU:"
    echo "  - n1-standard-16 + 1x NVIDIA T4"
    echo "  - Koszt: ~$1-2 za pełny WGS"
    echo "  - Czas: ~1-2h"
}

#===============================================================================
# OPCJA 4: Cloud (najszybsza dla WGS)
#===============================================================================

setup_cloud() {
    echo ""
    echo "=== OPCJA 4: Google Cloud / AWS (ZALECANE dla szybkości) ==="
    echo ""
    
    cat << 'EOF'
GOOGLE CLOUD z GPU:
-------------------
# 1. Utwórz VM z GPU
gcloud compute instances create deepvariant-vm \
    --zone=us-central1-a \
    --machine-type=n1-standard-16 \
    --accelerator=type=nvidia-tesla-t4,count=1 \
    --image-family=deepvariant-1-6-0-gpu \
    --image-project=deepvariant-docker \
    --boot-disk-size=200GB

# 2. Skopiuj dane
gcloud compute scp saryd.bam deepvariant-vm:~/
gcloud compute scp human_g1k_v37.fasta* deepvariant-vm:~/

# 3. SSH i uruchom
gcloud compute ssh deepvariant-vm

# Na VM:
/opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=human_g1k_v37.fasta \
    --reads=saryd.bam \
    --output_vcf=output.vcf.gz \
    --num_shards=16

# 4. Pobierz wyniki i usuń VM
gcloud compute scp deepvariant-vm:~/output.vcf.gz ./
gcloud compute instances delete deepvariant-vm

Szacowany koszt: $2-5
Szacowany czas: 1-2h

AWS z GPU:
----------
# Użyj AMI: Deep Learning AMI (Ubuntu)
# Instance: g4dn.4xlarge (1x T4 GPU)
# Koszt: ~$1.20/h

EOF
}

#===============================================================================
# Menu
#===============================================================================

echo ""
echo "Wybierz opcję instalacji:"
echo "  1) Docker (najprostsza)"
echo "  2) Singularity"
echo "  3) Native Python (eksperymentalne)"
echo "  4) Cloud - Google/AWS (najszybsza)"
echo "  5) Pokaż wszystkie opcje"
echo ""

read -p "Wybór [1-5]: " choice

case $choice in
    1) setup_docker ;;
    2) setup_singularity ;;
    3) setup_native ;;
    4) setup_cloud ;;
    5) 
        setup_docker
        setup_singularity
        setup_native
        setup_cloud
        ;;
    *) echo "Nieprawidłowy wybór" ;;
esac

echo ""
echo "═══════════════════════════════════════════════════════════════════"
echo ""
echo "PORÓWNANIE: bcftools vs DeepVariant"
echo ""
echo "┌─────────────────┬──────────────┬──────────────┬─────────────────┐"
echo "│ Aspekt          │ bcftools     │ DeepVariant  │ GATK            │"
echo "├─────────────────┼──────────────┼──────────────┼─────────────────┤"
echo "│ Dokładność SNP  │ ~99.0%       │ ~99.7%       │ ~99.5%          │"
echo "│ Dokładność Indel│ ~95%         │ ~99%         │ ~97%            │"
echo "│ Czas (WGS)      │ 1-2h         │ 2-4h (GPU)   │ 6-12h           │"
echo "│ RAM             │ 8-16GB       │ 16-32GB      │ 32-64GB         │"
echo "│ GPU             │ Nie          │ TAK (opcja)  │ Nie             │"
echo "│ Instalacja      │ Łatwa        │ Średnia      │ Trudna          │"
echo "└─────────────────┴──────────────┴──────────────┴─────────────────┘"
echo ""
echo "REKOMENDACJA dla M4 MacBook Air:"
echo "  • Szybkość: bcftools parallel (turbo_vcf.sh)"
echo "  • Dokładność: DeepVariant przez Docker lub Cloud"
echo ""
