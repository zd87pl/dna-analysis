#!/bin/bash
#===============================================================================
# HELIXIGHT QUICK INSTALLER
# Run this script to set up Helixight on your system
#===============================================================================

set -e

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    HELIXIGHT QUICK INSTALLER                              ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""

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

echo "Detected OS: $OS"
echo ""

# Install dependencies
echo "Installing dependencies..."

case $OS in
    ubuntu|debian|pop|linuxmint)
        sudo apt-get update
        sudo apt-get install -y bcftools samtools tabix wget gzip curl
        ;;
    fedora)
        sudo dnf install -y bcftools samtools htslib wget gzip curl
        ;;
    centos|rhel|rocky|almalinux)
        sudo yum install -y epel-release
        sudo yum install -y bcftools samtools htslib wget gzip curl
        ;;
    arch|manjaro)
        sudo pacman -S --noconfirm bcftools samtools htslib wget gzip curl
        ;;
    macos)
        if ! command -v brew &> /dev/null; then
            echo "Installing Homebrew..."
            /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        fi
        brew install bcftools samtools htslib wget
        ;;
    *)
        echo "Unknown OS. Please install manually:"
        echo "  - bcftools"
        echo "  - samtools"
        echo "  - tabix (htslib)"
        echo "  - wget"
        exit 1
        ;;
esac

# Make all scripts executable
echo ""
echo "Setting up Helixight..."
chmod +x helixight.sh
chmod +x scripts/*.sh

# Add to PATH (optional)
read -p "Add Helixight to PATH? (y/n) [n]: " add_path
if [[ "$add_path" =~ ^[Yy]$ ]]; then
    INSTALL_DIR="$(pwd)"
    
    # Detect shell
    if [ -n "$ZSH_VERSION" ]; then
        SHELL_RC="$HOME/.zshrc"
    else
        SHELL_RC="$HOME/.bashrc"
    fi
    
    echo "" >> "$SHELL_RC"
    echo "# Helixight Genetic Analysis Toolkit" >> "$SHELL_RC"
    echo "export PATH=\"\$PATH:$INSTALL_DIR\"" >> "$SHELL_RC"
    echo "Added to $SHELL_RC"
    echo "Run 'source $SHELL_RC' or restart terminal"
fi

echo ""
echo "╔═══════════════════════════════════════════════════════════════════════════╗"
echo "║                    INSTALLATION COMPLETE! ✅                              ║"
echo "╚═══════════════════════════════════════════════════════════════════════════╝"
echo ""
echo "To start Helixight, run:"
echo "  ./helixight.sh"
echo ""
echo "Or run individual analyses:"
echo "  ./scripts/athletic_genetics.sh your_variants.vcf.gz"
echo "  ./scripts/triathlon_genetics.sh your_variants.vcf.gz"
echo ""
