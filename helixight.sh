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
#    Version: 1.3.0
#    License: MIT
#    Repository: https://github.com/zd87pl/dna-analysis
#
#===============================================================================

set -o pipefail

VERSION="1.3.0"
readonly VERSION
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$SCRIPT_DIR/scripts"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
# shellcheck disable=SC2034  # kept for menu colouring parity
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color

# Global variables
VCF_FILE=""
BAM_FILE=""
REFERENCE=""
OUTPUT_DIR="helixight_results"
LANGUAGE="${HELIXIGHT_LANGUAGE:-en}"
[[ "$LANGUAGE" == "pl" ]] || LANGUAGE="en"
export HELIXIGHT_LANGUAGE="$LANGUAGE"

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

# Validate file path - prevent command injection
validate_path() {
    local path="$1"
    local type="$2"  # "file", "dir" or "any" (character check only)

    # Reject shell metacharacters, newlines and a leading '-' (would be parsed
    # as an option by bcftools/samtools). Single-quoted so the backtick is literal.
    local re='[;|&$`<>]'
    if [[ -z "$path" || "$path" =~ $re || "$path" == -* || "$path" == *$'\n'* ]]; then
        echo -e "${RED}Error: Invalid characters in path${NC}" >&2
        return 1
    fi

    # Check path exists
    if [[ "$type" == "file" ]] && [[ ! -f "$path" ]]; then
        echo -e "${RED}Error: File not found: $path${NC}" >&2
        return 1
    elif [[ "$type" == "dir" ]] && [[ ! -d "$path" ]]; then
        echo -e "${RED}Error: Directory not found: $path${NC}" >&2
        return 1
    fi

    return 0
}

# Validate VCF file format
validate_vcf() {
    local vcf="$1"

    if [[ ! -f "$vcf" ]]; then
        return 1
    fi

    # Check file extension
    if [[ ! "$vcf" =~ \.(vcf|vcf\.gz)$ ]]; then
        echo -e "${YELLOW}Warning: File does not have .vcf or .vcf.gz extension${NC}" >&2
    fi

    # Basic format check - bcftools must be able to read the header
    # (plain, gzip or BGZF). No zcat|head pipeline: that SIGPIPEs under pipefail.
    if ! bcftools view -h "$vcf" >/dev/null 2>&1; then
        echo -e "${RED}Error: File does not appear to be a valid VCF${NC}" >&2
        return 1
    fi

    return 0
}

# Validate BAM file format
validate_bam() {
    local bam="$1"

    if [[ ! -f "$bam" ]]; then
        return 1
    fi

    # Check file extension
    if [[ ! "$bam" =~ \.bam$ ]]; then
        echo -e "${YELLOW}Warning: File does not have .bam extension${NC}" >&2
    fi

    # Use samtools to validate
    if ! samtools quickcheck "$bam" 2>/dev/null; then
        echo -e "${RED}Error: File does not appear to be a valid BAM${NC}" >&2
        return 1
    fi

    return 0
}

#===============================================================================
# FUNCTIONS
#===============================================================================

print_banner() {
    [[ -t 1 ]] && clear
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
    echo -e "    ${GREEN}13)${NC} $(msg_en "🧬 Mega Analysis (~80 variants, 10 sections)" "🧬 Mega Analiza (~80 wariantów, 10 sekcji)")"
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
    echo -e "    ${GREEN}40)${NC} $(msg_en "📄 Generate summary report (text; PDF if pandoc is installed)" "📄 Wygeneruj raport zbiorczy (tekst; PDF jeśli zainstalowano pandoc)")"
    echo -e "    ${GREEN}41)${NC} $(msg_en "🔄 Run the 9 VCF analyses (10-14, 20-23) and generate the summary report" "🔄 Uruchom 9 analiz VCF (10-14, 20-23) i wygeneruj raport zbiorczy")"
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
    
    local deps=("bcftools" "samtools" "bgzip" "tabix" "wget" "gzip" "bc")
    local missing=()
    
    for dep in "${deps[@]}"; do
        if command -v "$dep" &> /dev/null; then
            echo -e "  ${GREEN}✓${NC} $dep $(command -v "$dep")"
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
    
    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

install_dependencies() {
    echo -e "\n${YELLOW}$(msg_en "Installing dependencies..." "Instaluję zależności...")${NC}\n"
    
    # Detect OS. Never source /etc/os-release: it assigns VERSION, NAME, ID ...
    # (VERSION is readonly here, so sourcing would even fail). Parse ID instead.
    local OS
    if [ -f /etc/os-release ]; then
        OS=$(sed -n 's/^ID=//p' /etc/os-release | tr -d '"' | head -1)
        OS=${OS:-unknown}
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
            sudo apt-get install -y bcftools samtools tabix wget gzip curl bc
            ;;
        fedora|centos|rhel)
            echo -e "${YELLOW}$(msg_en "Installing via dnf/yum..." "Instaluję przez dnf/yum...")${NC}"
            sudo dnf install -y bcftools samtools htslib wget gzip curl bc || \
            sudo yum install -y bcftools samtools htslib wget gzip curl bc
            ;;
        arch|manjaro)
            echo -e "${YELLOW}$(msg_en "Installing via pacman..." "Instaluję przez pacman...")${NC}"
            sudo pacman -S --noconfirm bcftools samtools htslib wget gzip curl bc
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
            echo "  bcftools, samtools, tabix, wget, gzip, bc"
            ;;
    esac
    
    echo ""
    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

download_reference() {
    echo -e "\n${YELLOW}$(msg_en "Reference Genome Download" "Pobieranie genomu referencyjnego")${NC}\n"
    
    echo -e "$(msg_en "This will download GRCh38 reference genome (~3GB compressed, ~15GB uncompressed)" "To pobierze genom referencyjny GRCh38 (~3GB skompresowany, ~15GB rozpakowany)")"
    echo ""
    
    read -r -p "$(msg_en "Download directory [./reference]: " "Katalog pobierania [./reference]: ")" REF_DIR
    REF_DIR=${REF_DIR:-./reference}

    # Character check only: the directory may not exist yet
    if ! validate_path "$REF_DIR" "any"; then
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return 1
    fi
    if ! mkdir -p "$REF_DIR"; then
        echo -e "${RED}$(msg_en "Cannot create directory:" "Nie można utworzyć katalogu:") $REF_DIR${NC}" >&2
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return 1
    fi
    
    echo ""
    echo -e "${YELLOW}$(msg_en "Select reference source:" "Wybierz źródło referencji:")${NC}"
    echo "  1) NCBI (recommended)"
    echo "  2) Ensembl"
    echo "  3) UCSC"
    
    read -r -p "$(msg_en "Choice [1]: " "Wybór [1]: ")" ref_choice
    ref_choice=${ref_choice:-1}

    local URL
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
        *)
            echo -e "${RED}$(msg_en "Invalid choice" "Nieprawidłowy wybór")${NC}" >&2
            read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
            return 1
            ;;
    esac

    local gz="$REF_DIR/GRCh38.fa.gz"
    local fa="$REF_DIR/GRCh38.fa"

    # Every step is checked: the "ready" banner must only appear when the
    # FASTA and its index really exist.
    if [[ -f "$fa" ]]; then
        echo -e "${YELLOW}$(msg_en "Reference already decompressed, skipping download:" "Referencja już rozpakowana, pomijam pobieranie:") $fa${NC}"
    else
        echo ""
        echo -e "${CYAN}$(msg_en "Downloading..." "Pobieram...")${NC}"
        if ! wget -c "$URL" -O "$gz"; then
            rm -f "$gz"
            echo -e "${RED}$(msg_en "Download failed" "Pobieranie nie powiodło się")${NC}" >&2
            read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
            return 1
        fi

        echo -e "${CYAN}$(msg_en "Decompressing..." "Rozpakowuję...")${NC}"
        if ! gunzip -k "$gz"; then
            rm -f "$fa"
            echo -e "${RED}$(msg_en "Decompression failed" "Rozpakowywanie nie powiodło się")${NC}" >&2
            read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
            return 1
        fi
    fi

    echo -e "${CYAN}$(msg_en "Indexing..." "Indeksuję...")${NC}"
    if ! samtools faidx "$fa"; then
        echo -e "${RED}$(msg_en "Indexing failed" "Indeksowanie nie powiodło się")${NC}" >&2
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return 1
    fi

    REFERENCE="$fa"
    echo ""
    echo -e "${GREEN}$(msg_en "Reference genome ready:" "Genom referencyjny gotowy:") $REFERENCE${NC}"

    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

select_vcf() {
    echo -e "\n${YELLOW}$(msg_en "Select VCF File" "Wybierz plik VCF")${NC}\n"

    # List VCF files in current directory - safe method using find
    echo -e "$(msg_en "VCF files in current directory:" "Pliki VCF w bieżącym katalogu:")"
    echo ""

    local vcf_files=()
    while IFS= read -r -d '' file; do
        vcf_files+=("$file")
    done < <(find . -maxdepth 1 -type f \( -name "*.vcf" -o -name "*.vcf.gz" \) -print0 2>/dev/null | sort -z)

    if [[ ${#vcf_files[@]} -eq 0 ]]; then
        echo -e "${YELLOW}$(msg_en "No VCF files found in current directory." "Nie znaleziono plików VCF w bieżącym katalogu.")${NC}"
    else
        for i in "${!vcf_files[@]}"; do
            # Remove leading ./ from path for display
            local display_name="${vcf_files[$i]#./}"
            echo "  $((i+1))) $display_name"
        done
    fi

    echo ""
    read -r -p "$(msg_en "Enter VCF file path (or number from list): " "Podaj ścieżkę do pliku VCF (lub numer z listy): ")" vcf_input

    # Sanitize input - remove leading/trailing whitespace
    vcf_input="${vcf_input#"${vcf_input%%[![:space:]]*}"}"
    vcf_input="${vcf_input%"${vcf_input##*[![:space:]]}"}"

    # Check if it's a number
    if [[ "$vcf_input" =~ ^[0-9]+$ ]] && [[ "$vcf_input" -le "${#vcf_files[@]}" ]] && [[ "$vcf_input" -gt 0 ]]; then
        VCF_FILE="${vcf_files[$((vcf_input-1))]}"
    else
        VCF_FILE="$vcf_input"
    fi

    # Validate path for dangerous characters
    if ! validate_path "$VCF_FILE" "file"; then
        VCF_FILE=""
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    # Validate VCF format
    if ! validate_vcf "$VCF_FILE"; then
        VCF_FILE=""
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    echo -e "${GREEN}$(msg_en "VCF file selected:" "Wybrano plik VCF:") $VCF_FILE${NC}"

    # Quick stats
    echo ""
    echo -e "$(msg_en "Quick statistics:" "Szybkie statystyki:")"
    local total
    total=$(bcftools view -H "$VCF_FILE" 2>/dev/null | wc -l)
    echo "  $(msg_en "Total variants:" "Całkowita liczba wariantów:") $total"

    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

create_vcf_from_bam() {
    echo -e "\n${YELLOW}$(msg_en "Create VCF from BAM File" "Stwórz VCF z pliku BAM")${NC}\n"

    read -r -p "$(msg_en "Enter BAM file path: " "Podaj ścieżkę do pliku BAM: ")" BAM_FILE

    # Validate path
    if ! validate_path "$BAM_FILE" "file"; then
        BAM_FILE=""
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    # Validate BAM format
    if ! validate_bam "$BAM_FILE"; then
        BAM_FILE=""
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    # Check for reference
    if [[ -z "$REFERENCE" ]] || [[ ! -f "$REFERENCE" ]]; then
        read -r -p "$(msg_en "Enter reference genome path: " "Podaj ścieżkę do genomu referencyjnego: ")" REFERENCE
        if ! validate_path "$REFERENCE" "file"; then
            echo -e "${RED}$(msg_en "Reference not found! Run option 3 first." "Nie znaleziono referencji! Najpierw uruchom opcję 3.")${NC}"
            read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
            return
        fi
    fi

    echo ""
    echo -e "${YELLOW}$(msg_en "Select variant calling method:" "Wybierz metodę variant calling:")${NC}"
    echo "  1) bcftools mpileup (fast, good quality)"
    echo "  2) DeepVariant via Docker (best quality, slower)"
    echo "  3) GATK HaplotypeCaller (classic, requires Java)"

    read -r -p "$(msg_en "Choice [1]: " "Wybór [1]: ")" method
    method=${method:-1}

    # Validate method choice
    if [[ ! "$method" =~ ^[1-3]$ ]]; then
        method=1
    fi

    local output_name
    output_name="$(basename "$BAM_FILE" .bam)_variants.vcf.gz"

    echo ""
    echo -e "${CYAN}$(msg_en "Starting variant calling..." "Rozpoczynam variant calling...")${NC}"
    echo -e "$(msg_en "This may take several hours for WGS data." "To może potrwać kilka godzin dla danych WGS.")"
    echo ""

    # The BAM's @SQ contig names must all exist in the reference: mpileup
    # exits 0 and writes an empty VCF for a contig the reference lacks.
    if [[ ! -f "$REFERENCE.fai" ]] && ! samtools faidx "$REFERENCE"; then
        echo -e "${RED}$(msg_en "Cannot index reference:" "Nie można zindeksować referencji:") $REFERENCE${NC}" >&2
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi
    local missing_contigs
    missing_contigs=$(samtools view -H "$BAM_FILE" | grep '^@SQ' | cut -f2 | sed 's/^SN://' \
        | grep -vxF -f <(cut -f1 "$REFERENCE.fai") | head -5 | tr '\n' ' ')
    if [[ -n "$missing_contigs" ]]; then
        echo -e "${RED}$(msg_en "BAM contigs not present in the reference (check chr prefix / genome build):" "Kontigi BAM nieobecne w referencji (sprawdź prefiks chr / wersję genomu):") $missing_contigs${NC}" >&2
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    case $method in
        1)
            # bcftools mpileup
            local THREADS
            THREADS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

            if ! bcftools mpileup \
                -Ou \
                -f "$REFERENCE" \
                --threads "$THREADS" \
                -q 20 \
                -Q 20 \
                "$BAM_FILE" | \
            bcftools call \
                -mv \
                --threads "$THREADS" \
                -Oz \
                -o "$output_name"; then
                echo -e "${RED}$(msg_en "Variant calling failed!" "Variant calling nie powiodło się!")${NC}"
                read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
                return
            fi

            bcftools index -t "$output_name" || true
            ;;
        2)
            # DeepVariant
            if ! command -v docker &> /dev/null; then
                echo -e "${RED}$(msg_en "Docker not found! Install Docker first." "Docker nie znaleziony! Najpierw zainstaluj Docker.")${NC}"
                read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
                return
            fi

            # Get absolute paths to avoid issues with spaces
            local abs_bam_dir abs_ref_dir abs_output_dir
            abs_bam_dir="$(cd "$(dirname "$BAM_FILE")" && pwd)"
            abs_ref_dir="$(cd "$(dirname "$REFERENCE")" && pwd)"
            abs_output_dir="$(pwd)"

            if ! docker run \
                -v "${abs_bam_dir}:/input:ro" \
                -v "${abs_output_dir}:/output" \
                -v "${abs_ref_dir}:/reference:ro" \
                google/deepvariant:latest \
                /opt/deepvariant/bin/run_deepvariant \
                --model_type=WGS \
                --ref="/reference/$(basename "$REFERENCE")" \
                --reads="/input/$(basename "$BAM_FILE")" \
                --output_vcf="/output/$output_name" \
                --num_shards="$(nproc 2>/dev/null || echo 4)"; then
                echo -e "${RED}$(msg_en "DeepVariant failed!" "DeepVariant nie powiodło się!")${NC}"
                read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
                return
            fi
            ;;
        3)
            echo -e "${YELLOW}$(msg_en "GATK requires manual setup. See documentation." "GATK wymaga ręcznej konfiguracji. Zobacz dokumentację.")${NC}"
            read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
            return
            ;;
    esac

    # A VCF with a header but zero records is a failure, not a result
    if [[ -f "$output_name" ]] && [[ -n "$(bcftools view -H "$output_name" 2>/dev/null | head -c1)" ]]; then
        VCF_FILE="$output_name"
        echo ""
        echo -e "${GREEN}$(msg_en "VCF created successfully:" "VCF utworzony pomyślnie:") $VCF_FILE${NC}"
    else
        echo -e "${RED}$(msg_en "VCF creation failed (no variants written)!" "Tworzenie VCF nie powiodło się (brak wariantów)!")${NC}"
    fi

    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

# Absolute path of an existing file (no realpath: absent on older macOS)
abs_path() {
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
}

# Per-input results directory: $OUTPUT_DIR/<input name without .vcf*/.bam>
run_dir_for() {
    local base
    base=$(basename "$1")
    base=${base%%.vcf*}
    base=${base%.bam}
    base=${base%.cram}
    echo "$OUTPUT_DIR/$base"
}

# Run scripts/$script inside the run directory, tee an ANSI-stripped copy of
# its output to <script>.txt there and return the script's exit status.
# Usage: run_script_in_dir INPUT SCRIPT [extra script args...]
run_script_in_dir() {
    local input="$1" script="$2"
    shift 2
    local input_abs run_dir rc esc
    input_abs=$(abs_path "$input")
    run_dir=$(run_dir_for "$input")
    mkdir -p "$run_dir" || return 1
    [ -x "$SCRIPTS_DIR/$script" ] || chmod +x "$SCRIPTS_DIR/$script"
    esc=$(printf '\033')
    (
        cd "$run_dir" || exit 1
        "$SCRIPTS_DIR/$script" "$input_abs" "$@" 2>&1 \
            | sed "s/${esc}\[[0-9;]*m//g" \
            | tee "${script%.sh}.txt"
        exit "${PIPESTATUS[0]}"
    )
    rc=$?
    if [ "$rc" -ne 0 ]; then
        echo -e "${RED}$(msg_en "Analysis failed (exit $rc):" "Analiza nie powiodła się (kod $rc):") $script${NC}" >&2
    else
        echo -e "${GREEN}$(msg_en "Results saved in:" "Wyniki zapisane w:") $run_dir${NC}"
    fi
    return "$rc"
}

# VCF analyses. Usage: run_analysis SCRIPT NAME [extra script args...]
run_analysis() {
    local script="$1"
    local name="$2"
    shift 2

    if [ -z "$VCF_FILE" ]; then
        echo -e "${RED}$(msg_en "No VCF file selected! Use option 5 first." "Nie wybrano pliku VCF! Najpierw użyj opcji 5.")${NC}"
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    if [ ! -f "$SCRIPTS_DIR/$script" ]; then
        echo -e "${RED}$(msg_en "Script not found:" "Nie znaleziono skryptu:") $script${NC}"
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    echo -e "\n${CYAN}$(msg_en "Running:" "Uruchamiam:") $name${NC}\n"

    run_script_in_dir "$VCF_FILE" "$script" "$@"

    echo ""
    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

# BAM-based analyses (STR, EVO2 data): these scripts take a BAM, not a VCF.
# Usage: run_bam_analysis SCRIPT NAME
run_bam_analysis() {
    local script="$1"
    local name="$2"

    if [ -z "$BAM_FILE" ]; then
        echo -e "${RED}$(msg_en "This analysis requires a BAM file. Set a BAM via option 4 first." "Ta analiza wymaga pliku BAM. Najpierw ustaw BAM w opcji 4.")${NC}"
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    if ! validate_bam "$BAM_FILE"; then
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    if [ ! -f "$SCRIPTS_DIR/$script" ]; then
        echo -e "${RED}$(msg_en "Script not found:" "Nie znaleziono skryptu:") $script${NC}"
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    echo -e "\n${CYAN}$(msg_en "Running:" "Uruchamiam:") $name${NC}\n"

    run_script_in_dir "$BAM_FILE" "$script"

    echo ""
    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

generate_pdf_report() {
    echo -e "\n${YELLOW}$(msg_en "Generating Summary Report" "Generowanie raportu zbiorczego")${NC}\n"

    if [ -z "$VCF_FILE" ]; then
        echo -e "${RED}$(msg_en "No VCF file selected! Use option 5 first." "Nie wybrano pliku VCF! Najpierw użyj opcji 5.")${NC}"
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    # Collect every text report written by run_analysis for this VCF
    local run_dir
    run_dir=$(run_dir_for "$VCF_FILE")
    local result_files=()
    if [ -d "$run_dir" ]; then
        while IFS= read -r -d '' file; do
            result_files+=("$file")
        done < <(find "$run_dir" -type f -name '*.txt' ! -name 'helixight_report_*' -print0 2>/dev/null | sort -z)
    fi

    if [ ${#result_files[@]} -eq 0 ]; then
        echo -e "${YELLOW}$(msg_en "No analysis results found. Run some analyses first!" "Nie znaleziono wyników analiz. Najpierw uruchom analizy!")${NC}"
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
        return
    fi

    # Combined report lives OUTSIDE the run directory so it is never
    # concatenated into the next report.
    local report_file
    report_file="$OUTPUT_DIR/helixight_report_$(basename "$run_dir")_$(date +%Y%m%d_%H%M%S).txt"

    {
        echo "================================================================================"
        echo "                    HELIXIGHT GENETIC ANALYSIS REPORT"
        echo "                    Generated: $(date)"
        echo "                    VCF File: $VCF_FILE"
        echo "================================================================================"
        echo ""
        
        for file in "${result_files[@]}"; do
            echo ""
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "SECTION: ${file#"$run_dir"/}"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            cat "$file" 2>/dev/null || echo "(No text reports found)"
            echo ""
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
        if pandoc "$report_file" -o "$pdf_file"; then
            echo -e "${GREEN}$(msg_en "PDF generated:" "PDF wygenerowany:") $pdf_file${NC}"
        else
            echo -e "${YELLOW}$(msg_en "pandoc could not produce a PDF (a PDF engine such as LaTeX is required); the text report above is complete" "pandoc nie zdołał utworzyć PDF (wymagany silnik PDF, np. LaTeX); raport tekstowy powyżej jest kompletny")${NC}"
        fi
    else
        echo -e "${YELLOW}$(msg_en "Install pandoc for PDF generation: sudo apt install pandoc" "Zainstaluj pandoc dla generowania PDF: sudo apt install pandoc")${NC}"
    fi
    
    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

run_all_analyses() {
    if [ -z "$VCF_FILE" ]; then
        echo -e "${RED}$(msg_en "No VCF file selected! Use option 5 first." "Nie wybrano pliku VCF! Najpierw użyj opcji 5.")${NC}"
        read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
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
    
    local failed=()
    for script in "${scripts[@]}"; do
        if [ -f "$SCRIPTS_DIR/$script" ]; then
            echo -e "\n${YELLOW}━━━ Running: $script ━━━${NC}\n"
            run_script_in_dir "$VCF_FILE" "$script" || failed+=("$script")
        fi
    done

    echo ""
    if [ ${#failed[@]} -gt 0 ]; then
        echo -e "${RED}$(msg_en "Analyses that failed:" "Analizy zakończone błędem:") ${failed[*]}${NC}"
    else
        echo -e "${GREEN}$(msg_en "All analyses complete!" "Wszystkie analizy zakończone!")${NC}"
    fi
    
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
  5. Generate summary report (option 40; text, PDF only if pandoc is installed)

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
     - ~80 genetic variants in 10 sections (Mega Analysis)
     - Pharmacogenomics, longevity, mental health
  
  🔬 ClinVar Scan
     - Known pathogenic variants from ClinVar database
     - Cancer predisposition, carrier status
  
  📊 Polygenic Risk Scores
     - Combined risk scores for complex diseases
     - Heart disease, diabetes, Alzheimer's

$(msg_en "MORE INFO:" "WIĘCEJ INFORMACJI:")
  - GitHub: https://github.com/zd87pl/dna-analysis
  - Documentation: https://helixight.github.io/docs

EOF
    
    read -r -p "$(msg_en "Press Enter to continue..." "Naciśnij Enter aby kontynuować...")"
}

toggle_language() {
    if [ "$LANGUAGE" == "en" ]; then
        LANGUAGE="pl"
    else
        LANGUAGE="en"
    fi
    export HELIXIGHT_LANGUAGE="$LANGUAGE"
}

#===============================================================================
# MAIN LOOP
#===============================================================================

main() {
    while true; do
        print_banner
        print_menu
        
        # EOF (closed/non-interactive stdin) exits instead of looping forever;
        # [[ -n $choice ]] keeps a last choice without trailing newline working.
        read -r -p "$(msg_en "Enter your choice: " "Wybierz opcję: ")" choice || [[ -n $choice ]] || { echo; exit 0; }
        
        case $choice in
            1) install_dependencies ;;
            2) check_dependencies ;;
            3) download_reference ;;
            4) create_vcf_from_bam ;;
            5) select_vcf ;;
            
            10) run_analysis "athletic_genetics.sh" "Athletic & Fitness Genetics" ;;
            11) run_analysis "triathlon_genetics.sh" "Triathlon Predisposition Analysis" ;;
            12) run_analysis "personalized_protocol.sh" "Personalized Protocol" ;;
            13) run_analysis "mega_analysis.sh" "Mega Analysis (~80 variants, 10 sections)" ;;
            14) run_analysis "fun_genetics.sh" "Fun Genetics" ;;
            
            20) run_analysis "clinvar_scan.sh" "ClinVar Pathogenic Variant Scan" ;;
            21) run_analysis "cancer_genes_scan.sh" "Cancer Genes Analysis" ;;
            22) run_analysis "polygenic_risk_scores.sh" "Polygenic Risk Scores" ;;
            23) run_analysis "rare_variants.sh" "Rare Variants Analysis" ;;
            
            30) run_bam_analysis "str_analysis.sh" "STR/Repeat Analysis" ;;
            31) run_analysis "extract_mtdna_haplogrep.sh" "mtDNA Haplogroup" "mtdna_analysis" ;;
            32) run_bam_analysis "prepare_evo2_data.sh" "Prepare EVO2 Data" ;;
            
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
