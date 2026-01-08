#!/bin/bash
#===============================================================================
# HELIXIGHT QUICK INSTALLER
# Run this script to set up Helixight on your system
#===============================================================================

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

error_exit() {
    echo -e "${RED}Error: $1${NC}" >&2
    exit 1
}

warn() {
    echo -e "${YELLOW}Warning: $1${NC}" >&2
}

success() {
    echo -e "${GREEN}$1${NC}"
}

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    HELIXIGHT QUICK INSTALLER                              ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

# Check we're in the correct directory
if [[ ! -f "helixight.sh" ]]; then
    error_exit "helixight.sh not found. Please run this script from the Helixight directory."
fi

if [[ ! -d "scripts" ]]; then
    error_exit "scripts directory not found. Please run this script from the Helixight directory."
fi

# Detect OS
OS="unknown"
if [[ -f /etc/os-release ]]; then
    # shellcheck source=/dev/null
    . /etc/os-release
    OS="$ID"
elif [[ -f /etc/debian_version ]]; then
    OS="debian"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
fi

echo "Detected OS: $OS"
echo ""

# Install dependencies
echo "Installing dependencies..."

install_success=true

case $OS in
    ubuntu|debian|pop|linuxmint)
        if ! sudo apt-get update; then
            warn "apt-get update failed, continuing anyway..."
        fi
        if ! sudo apt-get install -y bcftools samtools tabix wget gzip curl; then
            install_success=false
        fi
        ;;
    fedora)
        if ! sudo dnf install -y bcftools samtools htslib wget gzip curl; then
            install_success=false
        fi
        ;;
    centos|rhel|rocky|almalinux)
        sudo yum install -y epel-release 2>/dev/null || true
        if ! sudo yum install -y bcftools samtools htslib wget gzip curl; then
            install_success=false
        fi
        ;;
    arch|manjaro)
        if ! sudo pacman -S --noconfirm bcftools samtools htslib wget gzip curl; then
            install_success=false
        fi
        ;;
    macos)
        if ! command -v brew &> /dev/null; then
            echo "Installing Homebrew..."
            if ! /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"; then
                error_exit "Failed to install Homebrew"
            fi
        fi
        if ! brew install bcftools samtools htslib wget; then
            install_success=false
        fi
        ;;
    *)
        echo "Unknown OS. Please install manually:"
        echo "  - bcftools"
        echo "  - samtools"
        echo "  - tabix (htslib)"
        echo "  - wget"
        install_success=false
        ;;
esac

if [[ "$install_success" == "false" ]]; then
    warn "Some dependencies may not have been installed correctly."
fi

# Make all scripts executable
echo ""
echo "Setting up Helixight..."

if [[ -f "helixight.sh" ]]; then
    chmod +x helixight.sh
    success "Made helixight.sh executable"
else
    warn "helixight.sh not found"
fi

# Make scripts executable safely
script_count=0
for script in scripts/*.sh; do
    if [[ -f "$script" ]]; then
        chmod +x "$script"
        ((script_count++)) || true
    fi
done

if [[ $script_count -gt 0 ]]; then
    success "Made $script_count analysis scripts executable"
else
    warn "No scripts found in scripts/ directory"
fi

# Add to PATH (optional)
echo ""
read -r -p "Add Helixight to PATH? (y/n) [n]: " add_path
if [[ "${add_path,,}" =~ ^y ]]; then
    INSTALL_DIR="$(pwd)"

    # Detect shell config file
    SHELL_RC=""
    if [[ -n "${ZSH_VERSION:-}" ]] || [[ "$SHELL" == *"zsh"* ]]; then
        SHELL_RC="$HOME/.zshrc"
    elif [[ -f "$HOME/.bashrc" ]]; then
        SHELL_RC="$HOME/.bashrc"
    elif [[ -f "$HOME/.bash_profile" ]]; then
        SHELL_RC="$HOME/.bash_profile"
    fi

    if [[ -n "$SHELL_RC" ]] && [[ -f "$SHELL_RC" ]]; then
        # Check if already added
        if grep -q "Helixight Genetic Analysis Toolkit" "$SHELL_RC" 2>/dev/null; then
            echo "Helixight already in PATH configuration"
        else
            {
                echo ""
                echo "# Helixight Genetic Analysis Toolkit"
                echo "export PATH=\"\$PATH:$INSTALL_DIR\""
            } >> "$SHELL_RC"
            success "Added to $SHELL_RC"
            echo "Run 'source $SHELL_RC' or restart terminal"
        fi
    else
        warn "Could not detect shell config file. Add manually:"
        echo "  export PATH=\"\$PATH:$INSTALL_DIR\""
    fi
fi

# Verify installation
echo ""
echo "Verifying installation..."
missing_deps=()
for dep in bcftools samtools tabix wget; do
    if command -v "$dep" &> /dev/null; then
        success "  ✓ $dep found"
    else
        echo -e "  ${RED}✗ $dep not found${NC}"
        missing_deps+=("$dep")
    fi
done

echo ""
if [[ ${#missing_deps[@]} -gt 0 ]]; then
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║             INSTALLATION INCOMPLETE - Missing dependencies               ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Please install: ${missing_deps[*]}"
else
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║                    INSTALLATION COMPLETE! ✅                              ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"
fi

echo ""
echo "To start Helixight, run:"
echo "  ./helixight.sh"
echo ""
echo "Or run individual analyses:"
echo "  ./scripts/athletic_genetics.sh your_variants.vcf.gz"
echo "  ./scripts/triathlon_genetics.sh your_variants.vcf.gz"
echo ""
