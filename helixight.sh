#!/bin/bash
#===============================================================================
#
#    ██╗  ██╗███████╗██╗     ██╗██╗  ██╗██╗ ██████╗ ██╗  ██╗████████╗
#    ██║  ██║██╔════╝██║     ██║╚██╗██╔╝██║██╔════╝ ██║  ██║╚══██╔══╝
#    ███████║█████╗  ██║     ██║ ╚███╔╝ ██║██║  ███╗███████║   ██║   
#    ██╔══██║██╔══╝  ██║     ██║ ██╔██╗ ██║██║   ██║██╔══██║   ██║   
#    ██║  ██║███████╗███████╗██║██╔╝ ██╗██║╚██████╔╝██║  ██║   ██║   
#    ╚═╝  ╚═╝╚══════╝╚══════╝╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝ ╚═╝  ╚═╝   ╚═╝   
#
#    OPEN SOURCE GENETIC ANALYSIS TOOLKIT
#    From BAM to Actionable Insights
#
#    Version: 1.0.0
#    License: MIT
#    Repository: https://github.com/helixight/helixight-oss
#
#===============================================================================

VERSION="1.0.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color

# Global variables
VCF_FILE=""
BAM_FILE=""
REFERENCE=""
OUTPUT_DIR="helixight_results"
LANGUAGE="en"

#===============================================================================
# FUNCTIONS
#===============================================================================

print_banner() {
    clear
    echo -e "${CYAN}"
    cat << "EOF"
    ╔═══════════════════════════════════════════════════════════════════════════╗
    ║                                                                           ║
    ║     ██╗  ██╗███████╗██╗     ██╗██╗  ██╗██╗ ██████╗ ██╗  ██╗████████╗      ║
    ║     ██║  ██║██╔════╝██║     ██║╚██╗██╔╝██║██╔════╝ ██║  ██║╚══██╔══╝      ║
    ║     ███████║█████╗  ██║     ██║ ╚███╔╝ ██║██║  ███╗███████║   ██║         ║
    ║     ██╔══██║██╔══╝  ██║     ██║ ██╔██╗ ██║██║   ██║██╔══██║   ██║         ║
    ║     ██║  ██║███████╗███████╗██║██╔╝ ██╗██║╚██████╔╝██║  ██║   ██║         ║
    ║     ╚═╝  ╚═╝╚══════╝╚══════╝╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝ ╚═╝  ╚═╝   ╚═╝         ║
    ║                                                                           ║
    ║              🧬 OPEN SOURCE GENETIC ANALYSIS TOOLKIT 🧬                   ║
    ║                      From BAM to Actionable Insights                      ║
    ║                                                                           ║
    ╚═══════════════════════════════════════════════════════════════════════════╝
EOF
    echo -e "${NC}"
    echo -e "    ${WHITE}Version: $VERSION${NC}    ${BLUE}Language: $LANGUAGE${NC}"
    echo ""
}

msg_en() { [ "$LANGUAGE" == "en" ] && echo "$1" || echo "$2"; }

print_menu() {
    echo -e "${WHITE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo ""
    echo -e "  ${YELLOW}$(msg_en "SETUP & INSTALLATION" "KONFIGURACJA I INSTALACJA")${NC}"
    echo -e "    ${GREEN}1)${NC} $(msg_en "Install dependencies (bcftools, samtools, etc.)" "Zainstaluj zależności (bcftools, samtools, itp.)")"
    echo -e "    ${GREEN}2)${NC} $(msg_en "Check system requirements" "Sprawdź wymagania systemowe")"
    echo -e "    ${GREEN}3)${NC} $(msg_en "Download reference genome (GRCh38)" "Pobierz genom referencyjny (GRCh38)")"
    echo ""
    echo -e "  ${YELLOW}$(msg_en "DATA PREPARATION" "PRZYGOTOWANIE DANYCH")${NC}"
    echo -e "    ${GREEN}4)${NC} $(msg_en "Create VCF from BAM file (variant calling)" "Stwórz VCF z pliku BAM (variant calling)")"
    echo -e "    ${GREEN}5)${NC} $(msg_en "Select existing VCF file" "Wybierz istniejący plik VCF")"
    echo ""
    echo -e "  ${YELLOW}$(msg_en "GENETIC ANALYSIS" "ANALIZY GENETYCZNE")${NC}"
    echo -e "    ${GREEN}10)${NC} $(msg_en "🏃 Athletic & Fitness Genetics" "🏃 Genetyka sportowa i fitness")"
    echo -e "    ${GREEN}11)${NC} $(msg_en "🏊 Triathlon Predisposition Analysis" "🏊 Analiza predyspozycji do triathlonu")"
    echo -e "    ${GREEN}12)${NC} $(msg_en "💊 Personalized Protocol (training, supplements)" "💊 Spersonalizowany protokół (trening, suplementy)")"
    echo -e "    ${GREEN}13)${NC} $(msg_en "🧬 Mega Analysis (500+ variants)" "🧬 Mega Analiza (500+ wariantów)")"
    echo -e "    ${GREEN}14)${NC} $(msg_en "🎲 Fun Genetics (traits, ancestry hints)" "🎲 Zabawna Genetyka (cechy, pochodzenie)")"
    echo ""
    echo -e "  ${YELLOW}$(msg_en "HEALTH & CLINICAL" "ZDROWIE I KLINIKA")${NC}"
    echo -e "    ${GREEN}20)${NC} $(msg_en "🔬 ClinVar Pathogenic Variant Scan" "🔬 Skan patogennych wariantów ClinVar")"
    echo -e "    ${GREEN}21)${NC} $(msg_en "🎀 Cancer Genes Analysis (BRCA1/2, TP53, Lynch)" "🎀 Analiza genów rakowych (BRCA1/2, TP53, Lynch)")"
    echo -e "    ${GREEN}22)${NC} $(msg_en "📊 Polygenic Risk Scores" "📊 Poligeniczne skale ryzyka")"
    echo -e "    ${GREEN}23)${NC} $(msg_en "💎 Rare Variants Analysis" "💎 Analiza rzadkich wariantów")"
    echo ""
    echo -e "  ${YELLOW}$(msg_en "ADVANCED" "ZAAWANSOWANE")${NC}"
    echo -e "    ${GREEN}30)${NC} $(msg_en "🧬 STR/Repeat Analysis (requires BAM)" "🧬 Analiza STR/powtórzeń (wymaga BAM)")"
    echo -e "    ${GREEN}31)${NC} $(msg_en "🌍 mtDNA Haplogroup Analysis" "🌍 Analiza haplogrupy mtDNA")"
    echo -e "    ${GREEN}32)${NC} $(msg_en "🤖 Prepare data for EVO2 AI analysis" "🤖 Przygotuj dane do analizy EVO2 AI")"
    echo ""
    echo -e "  ${YELLOW}$(msg_en "REPORTS" "RAPORTY")${NC}"
    echo -e "    ${GREEN}40)${NC} $(msg_en "📄 Generate PDF Summary Report" "📄 Wygeneruj raport PDF")"
    echo -e "    ${GREEN}41)${NC} $(msg_en "🔄 Run ALL analyses and generate full report" "🔄 Uruchom WSZYSTKIE analizy i wygeneruj pełny raport")"
    echo ""
    echo -e "  ${YELLOW}$(msg_en "OTHER" "INNE")${NC}"
    echo -e "    ${GREEN}L)${NC}  $(msg_en "Change language (EN/PL)" "Zmień język (EN/PL)")"
    echo -e "    ${GREEN}H)${NC}  $(msg_en "Help & Documentation" "Pomoc i dokumentacja")"
    echo -e "    ${GREEN}Q)${NC}  $(msg_en "Quit" "Wyjście")"
    echo ""
    echo -e "${WHITE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    
    # Show current file status
    echo ""
    if [ -n "$VCF_FILE" ]; then
        echo -e "  ${GREEN}✓${NC} VCF: ${CYAN}$VCF_FILE${NC}"
    else
        echo -e "  ${RED}✗${NC} $(msg_en "No VCF file selected" "Nie wybrano pliku VCF")"
    fi
    if [ -n "$BAM_FILE" ]; then
        echo -e "  ${GREEN}✓${NC} BAM: ${CYAN}$BAM_FILE${NC}"
    fi
    echo ""
}

check_dependencies() {
    echo -e "\n${YELLOW}$(msg_en "Checking dependencies..." "Sprawdzam zależności...")${NC}\n"
    
    local deps=("bcftools" "samtools" "bgzip" "tabix" "wget" "gzip")
    local missing=()
    
    for dep in "${deps[@]}"; do
        if command -v "$dep" &> /dev/null; then
            echo -e "  ${GREEN}✓${NC} $dep $(command -v $dep)"
        else
            echo -e "  ${RED}✗${NC} $dep - $(msg_en "NOT FOUND" "NIE ZNALEZIONO")"
            missing+=("$dep")
        fi
    done
    
    # Check optional
    echo ""
    echo -e "${YELLOW}$(msg_en "Optional tools:" "Opcjonalne narzędzia:")${NC}"
    
    for dep in "deepvariant" "docker" "ExpansionHunter"; do
        if command -v "$dep" &> /dev/null; then
            echo -e "  ${GREEN}✓${NC} $dep"
        else
            echo -e "  ${YELLOW}○${NC} $dep - $(msg_en "not installed (optional)" "nie zainstalowany (opcjonalny)")"
        fi
    done
    
    echo ""
    if [ ${#missing[@]} -gt 0 ]; then
        echo -e "${RED}$(msg_en "Missing required dependencies!" "Brakuje wymaganych zależności!")${NC}"
        echo -e "$(msg_en "Run option 1 to install them." "Uruchom opcję 1 aby je zainstalować.")"
    else
        echo -e "${GREEN}$(msg_en "All required dependencies are installed!" "Wszystkie wymagane zależności są zainstalowane!")${NC}"
    fi
    
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

install_dependencies() {
    echo -e "\n${YELLOW}$(msg_en "Installing dependencies..." "Instaluję zależności...")${NC}\n"
    
    # Detect OS
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        OS=$ID
    elif [ -f /etc/debian_version ]; then
        OS="debian"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        OS="macos"
    else
        OS="unknown"
    fi
    
    echo -e "$(msg_en "Detected OS:" "Wykryty system:") ${CYAN}$OS${NC}"
    echo ""
    
    case $OS in
        ubuntu|debian|pop)
            echo -e "${YELLOW}$(msg_en "Installing via apt..." "Instaluję przez apt...")${NC}"
            sudo apt-get update
            sudo apt-get install -y bcftools samtools tabix wget gzip curl
            ;;
        fedora|centos|rhel)
            echo -e "${YELLOW}$(msg_en "Installing via dnf/yum..." "Instaluję przez dnf/yum...")${NC}"
            sudo dnf install -y bcftools samtools htslib wget gzip curl || \
            sudo yum install -y bcftools samtools htslib wget gzip curl
            ;;
        arch|manjaro)
            echo -e "${YELLOW}$(msg_en "Installing via pacman..." "Instaluję przez pacman...")${NC}"
            sudo pacman -S --noconfirm bcftools samtools htslib wget gzip curl
            ;;
        macos)
            echo -e "${YELLOW}$(msg_en "Installing via Homebrew..." "Instaluję przez Homebrew...")${NC}"
            if ! command -v brew &> /dev/null; then
                echo -e "${RED}$(msg_en "Homebrew not found. Install from https://brew.sh" "Homebrew nie znaleziony. Zainstaluj z https://brew.sh")${NC}"
            else
                brew install bcftools samtools htslib wget
            fi
            ;;
        *)
            echo -e "${RED}$(msg_en "Unknown OS. Please install manually:" "Nieznany system. Zainstaluj ręcznie:")${NC}"
            echo "  bcftools, samtools, tabix, wget, gzip"
            ;;
    esac
    
    echo ""
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

download_reference() {
    echo -e "\n${YELLOW}$(msg_en "Reference Genome Download" "Pobieranie genomu referencyjnego")${NC}\n"
    
    echo -e "$(msg_en "This will download GRCh38 reference genome (~3GB compressed, ~15GB uncompressed)" "To pobierze genom referencyjny GRCh38 (~3GB skompresowany, ~15GB rozpakowany)")"
    echo ""
    
    read -p "$(msg_en "Download directory [./reference]: " "Katalog pobierania [./reference]: ")" REF_DIR
    REF_DIR=${REF_DIR:-./reference}
    
    mkdir -p "$REF_DIR"
    
    echo ""
    echo -e "${YELLOW}$(msg_en "Select reference source:" "Wybierz źródło referencji:")${NC}"
    echo "  1) NCBI (recommended)"
    echo "  2) Ensembl"
    echo "  3) UCSC"
    
    read -p "$(msg_en "Choice [1]: " "Wybór [1]: ")" ref_choice
    ref_choice=${ref_choice:-1}
    
    case $ref_choice in
        1)
            URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
            ;;
        2)
            URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
            ;;
        3)
            URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
            ;;
    esac
    
    echo ""
    echo -e "${CYAN}$(msg_en "Downloading..." "Pobieram...")${NC}"
    wget -c "$URL" -O "$REF_DIR/GRCh38.fa.gz"
    
    echo -e "${CYAN}$(msg_en "Decompressing..." "Rozpakowuję...")${NC}"
    gunzip -k "$REF_DIR/GRCh38.fa.gz"
    
    echo -e "${CYAN}$(msg_en "Indexing..." "Indeksuję...")${NC}"
    samtools faidx "$REF_DIR/GRCh38.fa"
    
    REFERENCE="$REF_DIR/GRCh38.fa"
    echo ""
    echo -e "${GREEN}$(msg_en "Reference genome ready:" "Genom referencyjny gotowy:") $REFERENCE${NC}"
    
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

select_vcf() {
    echo -e "\n${YELLOW}$(msg_en "Select VCF File" "Wybierz plik VCF")${NC}\n"
    
    # List VCF files in current directory
    echo -e "$(msg_en "VCF files in current directory:" "Pliki VCF w bieżącym katalogu:")"
    echo ""
    
    local vcf_files=($(ls -1 *.vcf *.vcf.gz 2>/dev/null))
    
    if [ ${#vcf_files[@]} -eq 0 ]; then
        echo -e "${YELLOW}$(msg_en "No VCF files found in current directory." "Nie znaleziono plików VCF w bieżącym katalogu.")${NC}"
    else
        for i in "${!vcf_files[@]}"; do
            echo "  $((i+1))) ${vcf_files[$i]}"
        done
    fi
    
    echo ""
    read -p "$(msg_en "Enter VCF file path (or number from list): " "Podaj ścieżkę do pliku VCF (lub numer z listy): ")" vcf_input
    
    # Check if it's a number
    if [[ "$vcf_input" =~ ^[0-9]+$ ]] && [ "$vcf_input" -le "${#vcf_files[@]}" ] && [ "$vcf_input" -gt 0 ]; then
        VCF_FILE="${vcf_files[$((vcf_input-1))]}"
    else
        VCF_FILE="$vcf_input"
    fi
    
    # Validate
    if [ -f "$VCF_FILE" ]; then
        echo -e "${GREEN}$(msg_en "VCF file selected:" "Wybrano plik VCF:") $VCF_FILE${NC}"
        
        # Quick stats
        echo ""
        echo -e "$(msg_en "Quick statistics:" "Szybkie statystyki:")"
        local total=$(bcftools view -H "$VCF_FILE" 2>/dev/null | wc -l)
        echo "  $(msg_en "Total variants:" "Całkowita liczba wariantów:") $total"
    else
        echo -e "${RED}$(msg_en "File not found:" "Nie znaleziono pliku:") $VCF_FILE${NC}"
        VCF_FILE=""
    fi
    
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

create_vcf_from_bam() {
    echo -e "\n${YELLOW}$(msg_en "Create VCF from BAM File" "Stwórz VCF z pliku BAM")${NC}\n"
    
    read -p "$(msg_en "Enter BAM file path: " "Podaj ścieżkę do pliku BAM: ")" BAM_FILE
    
    if [ ! -f "$BAM_FILE" ]; then
        echo -e "${RED}$(msg_en "BAM file not found!" "Nie znaleziono pliku BAM!")${NC}"
        read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi
    
    # Check for reference
    if [ -z "$REFERENCE" ] || [ ! -f "$REFERENCE" ]; then
        read -p "$(msg_en "Enter reference genome path: " "Podaj ścieżkę do genomu referencyjnego: ")" REFERENCE
        if [ ! -f "$REFERENCE" ]; then
            echo -e "${RED}$(msg_en "Reference not found! Run option 3 first." "Nie znaleziono referencji! Najpierw uruchom opcję 3.")${NC}"
            read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
            return
        fi
    fi
    
    echo ""
    echo -e "${YELLOW}$(msg_en "Select variant calling method:" "Wybierz metodę variant calling:")${NC}"
    echo "  1) bcftools mpileup (fast, good quality)"
    echo "  2) DeepVariant via Docker (best quality, slower)"
    echo "  3) GATK HaplotypeCaller (classic, requires Java)"
    
    read -p "$(msg_en "Choice [1]: " "Wybór [1]: ")" method
    method=${method:-1}
    
    local output_name=$(basename "$BAM_FILE" .bam)
    output_name="${output_name}_variants.vcf.gz"
    
    echo ""
    echo -e "${CYAN}$(msg_en "Starting variant calling..." "Rozpoczynam variant calling...")${NC}"
    echo -e "$(msg_en "This may take several hours for WGS data." "To może potrwać kilka godzin dla danych WGS.")"
    echo ""
    
    case $method in
        1)
            # bcftools mpileup
            THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
            
            bcftools mpileup \
                -Ou \
                -f "$REFERENCE" \
                --threads $THREADS \
                -q 20 \
                -Q 20 \
                "$BAM_FILE" | \
            bcftools call \
                -mv \
                --threads $THREADS \
                -Oz \
                -o "$output_name"
            
            bcftools index -t "$output_name"
            ;;
        2)
            # DeepVariant
            if ! command -v docker &> /dev/null; then
                echo -e "${RED}$(msg_en "Docker not found! Install Docker first." "Docker nie znaleziony! Najpierw zainstaluj Docker.")${NC}"
                read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
                return
            fi
            
            docker run \
                -v "$(pwd):/input" \
                -v "$(pwd):/output" \
                -v "$(dirname $REFERENCE):/reference" \
                google/deepvariant:latest \
                /opt/deepvariant/bin/run_deepvariant \
                --model_type=WGS \
                --ref="/reference/$(basename $REFERENCE)" \
                --reads="/input/$(basename $BAM_FILE)" \
                --output_vcf="/output/$output_name" \
                --num_shards=$(nproc)
            ;;
        3)
            echo -e "${YELLOW}$(msg_en "GATK requires manual setup. See documentation." "GATK wymaga ręcznej konfiguracji. Zobacz dokumentację.")${NC}"
            read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
            return
            ;;
    esac
    
    if [ -f "$output_name" ]; then
        VCF_FILE="$output_name"
        echo ""
        echo -e "${GREEN}$(msg_en "VCF created successfully:" "VCF utworzony pomyślnie:") $VCF_FILE${NC}"
    else
        echo -e "${RED}$(msg_en "VCF creation failed!" "Tworzenie VCF nie powiodło się!")${NC}"
    fi
    
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

run_analysis() {
    local script="$1"
    local name="$2"
    
    if [ -z "$VCF_FILE" ]; then
        echo -e "${RED}$(msg_en "No VCF file selected! Use option 5 first." "Nie wybrano pliku VCF! Najpierw użyj opcji 5.")${NC}"
        read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi
    
    if [ ! -f "$SCRIPTS_DIR/$script" ]; then
        echo -e "${RED}$(msg_en "Script not found:" "Nie znaleziono skryptu:") $script${NC}"
        read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi
    
    echo -e "\n${CYAN}$(msg_en "Running:" "Uruchamiam:") $name${NC}\n"
    
    chmod +x "$SCRIPTS_DIR/$script"
    "$SCRIPTS_DIR/$script" "$VCF_FILE"
    
    echo ""
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

generate_pdf_report() {
    echo -e "\n${YELLOW}$(msg_en "Generating PDF Report" "Generowanie raportu PDF")${NC}\n"
    
    # Check if we have any results
    local result_dirs=("athletic_analysis" "triathlon_analysis" "personalized_protocol" "mega_analysis_results" "clinvar_analysis" "cancer_genes_analysis")
    local has_results=false
    
    for dir in "${result_dirs[@]}"; do
        if [ -d "$dir" ]; then
            has_results=true
            break
        fi
    done
    
    if [ "$has_results" = false ]; then
        echo -e "${YELLOW}$(msg_en "No analysis results found. Run some analyses first!" "Nie znaleziono wyników analiz. Najpierw uruchom analizy!")${NC}"
        read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi
    
    # Create combined report
    local report_file="helixight_report_$(date +%Y%m%d_%H%M%S).txt"
    
    {
        echo "================================================================================"
        echo "                    HELIXIGHT GENETIC ANALYSIS REPORT"
        echo "                    Generated: $(date)"
        echo "                    VCF File: $VCF_FILE"
        echo "================================================================================"
        echo ""
        
        for dir in "${result_dirs[@]}"; do
            if [ -d "$dir" ]; then
                echo ""
                echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                echo "SECTION: $dir"
                echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
                echo ""
                cat "$dir"/*.txt 2>/dev/null || echo "(No text reports found)"
                echo ""
            fi
        done
        
        echo ""
        echo "================================================================================"
        echo "                              END OF REPORT"
        echo "================================================================================"
    } > "$report_file"
    
    echo -e "${GREEN}$(msg_en "Report generated:" "Raport wygenerowany:") $report_file${NC}"
    
    # Try to convert to PDF if pandoc is available
    if command -v pandoc &> /dev/null; then
        local pdf_file="${report_file%.txt}.pdf"
        pandoc "$report_file" -o "$pdf_file" 2>/dev/null && \
            echo -e "${GREEN}$(msg_en "PDF generated:" "PDF wygenerowany:") $pdf_file${NC}"
    else
        echo -e "${YELLOW}$(msg_en "Install pandoc for PDF generation: sudo apt install pandoc" "Zainstaluj pandoc dla generowania PDF: sudo apt install pandoc")${NC}"
    fi
    
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

run_all_analyses() {
    if [ -z "$VCF_FILE" ]; then
        echo -e "${RED}$(msg_en "No VCF file selected! Use option 5 first." "Nie wybrano pliku VCF! Najpierw użyj opcji 5.")${NC}"
        read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi
    
    echo -e "\n${CYAN}$(msg_en "Running ALL analyses..." "Uruchamiam WSZYSTKIE analizy...")${NC}\n"
    echo -e "$(msg_en "This will take some time..." "To zajmie trochę czasu...")"
    echo ""
    
    local scripts=(
        "athletic_genetics.sh"
        "triathlon_genetics.sh"
        "personalized_protocol.sh"
        "mega_analysis.sh"
        "clinvar_scan.sh"
        "cancer_genes_scan.sh"
        "polygenic_risk_scores.sh"
        "rare_variants.sh"
        "fun_genetics.sh"
    )
    
    for script in "${scripts[@]}"; do
        if [ -f "$SCRIPTS_DIR/$script" ]; then
            echo -e "\n${YELLOW}━━━ Running: $script ━━━${NC}\n"
            chmod +x "$SCRIPTS_DIR/$script"
            "$SCRIPTS_DIR/$script" "$VCF_FILE"
        fi
    done
    
    echo ""
    echo -e "${GREEN}$(msg_en "All analyses complete!" "Wszystkie analizy zakończone!")${NC}"
    
    # Generate report
    generate_pdf_report
}

show_help() {
    echo -e "\n${YELLOW}$(msg_en "HELIXIGHT - Help & Documentation" "HELIXIGHT - Pomoc i dokumentacja")${NC}\n"
    
    cat << EOF
$(msg_en "QUICK START:" "SZYBKI START:")

  1. Install dependencies (option 1)
  2. Download reference genome if you have BAM files (option 3)
  3. Either:
     - Create VCF from BAM (option 4), OR
     - Select existing VCF file (option 5)
  4. Run desired analyses (options 10-32)
  5. Generate PDF report (option 40)

$(msg_en "ANALYSIS DESCRIPTIONS:" "OPISY ANALIZ:")

  🏃 Athletic Genetics
     - ACTN3, ACE, PPARGC1A - muscle fiber type, endurance
     - Injury risk, recovery speed
  
  🏊 Triathlon Analysis
     - VO2max potential, fat metabolism, lactate threshold
     - Specific recommendations for Sprint/Olympic/Ironman
  
  💊 Personalized Protocol
     - Training recommendations based on genetics
     - Supplement stack based on MTHFR, VDR, COMT, etc.
     - What to avoid (caffeine, alcohol, etc.)
  
  🧬 Mega Analysis
     - 500+ genetic variants across 16 health domains
     - Pharmacogenomics, longevity, mental health
  
  🔬 ClinVar Scan
     - Known pathogenic variants from ClinVar database
     - Cancer predisposition, carrier status
  
  📊 Polygenic Risk Scores
     - Combined risk scores for complex diseases
     - Heart disease, diabetes, Alzheimer's

$(msg_en "MORE INFO:" "WIĘCEJ INFORMACJI:")
  - GitHub: https://github.com/helixight/helixight-oss
  - Documentation: https://helixight.github.io/docs

EOF
    
    read -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

toggle_language() {
    if [ "$LANGUAGE" == "en" ]; then
        LANGUAGE="pl"
    else
        LANGUAGE="en"
    fi
}

#===============================================================================
# MAIN LOOP
#===============================================================================

main() {
    while true; do
        print_banner
        print_menu
        
        read -p "$(msg_en "Enter your choice: " "Wybierz opcję: ")" choice
        
        case $choice in
            1) install_dependencies ;;
            2) check_dependencies ;;
            3) download_reference ;;
            4) create_vcf_from_bam ;;
            5) select_vcf ;;
            
            10) run_analysis "athletic_genetics.sh" "Athletic & Fitness Genetics" ;;
            11) run_analysis "triathlon_genetics.sh" "Triathlon Predisposition Analysis" ;;
            12) run_analysis "personalized_protocol.sh" "Personalized Protocol" ;;
            13) run_analysis "mega_analysis.sh" "Mega Analysis (500+ variants)" ;;
            14) run_analysis "fun_genetics.sh" "Fun Genetics" ;;
            
            20) run_analysis "clinvar_scan.sh" "ClinVar Pathogenic Variant Scan" ;;
            21) run_analysis "cancer_genes_scan.sh" "Cancer Genes Analysis" ;;
            22) run_analysis "polygenic_risk_scores.sh" "Polygenic Risk Scores" ;;
            23) run_analysis "rare_variants.sh" "Rare Variants Analysis" ;;
            
            30) run_analysis "str_analysis.sh" "STR/Repeat Analysis" ;;
            31) run_analysis "extract_mtdna_haplogrep.sh" "mtDNA Haplogroup" ;;
            32) run_analysis "prepare_evo2_data.sh" "Prepare EVO2 Data" ;;
            
            40) generate_pdf_report ;;
            41) run_all_analyses ;;
            
            [Ll]) toggle_language ;;
            [Hh]) show_help ;;
            [Qq]) 
                echo -e "\n${CYAN}$(msg_en "Thank you for using Helixight! 🧬" "Dziękujemy za korzystanie z Helixight! 🧬")${NC}\n"
                exit 0 
                ;;
            *)
                echo -e "${RED}$(msg_en "Invalid option!" "Nieprawidłowa opcja!")${NC}"
                sleep 1
                ;;
        esac
    done
}

# Run main
main
