# 🧬 Helixight - Open Source Genetic Analysis Toolkit

**Version 1.0.1**

<p align="center">
  <strong>From BAM to Actionable Insights</strong><br>
  Analyze your whole genome sequencing data locally, privately, and for free.
</p>

<p align="center">
  <a href="#-important-disclaimer">⚠️ Disclaimer</a> •
  <a href="#-features">Features</a> •
  <a href="#-quick-start">Quick Start</a> •
  <a href="#-analyses">Analyses</a> •
  <a href="#-contributing">Contributing</a>
</p>

---

## ⚠️ IMPORTANT DISCLAIMER

> **🚨 THIS SOFTWARE IS NOT FOR CLINICAL OR MEDICAL USE 🚨**
>
> Helixight is provided **for educational, research, and entertainment purposes only**.
>
> - ❌ **NOT** a medical device
> - ❌ **NOT** clinically validated
> - ❌ **NOT** a substitute for professional genetic counseling
> - ❌ **NOT** intended for diagnosing, treating, or preventing any disease
> - ❌ **NOT** approved by FDA, EMA, or any regulatory body
>
> **The results produced by this software:**
> - May contain errors, inaccuracies, or false positives/negatives
> - Should NOT be used for any medical or healthcare decisions
> - Do NOT account for all genetic factors, environmental influences, or gene-gene interactions
> - Are based on simplified interpretations of complex genetic science
>
> **If you discover concerning results:**
> - Do NOT panic - these are not clinical-grade results
> - Consult a certified genetic counselor or medical geneticist
> - Get proper clinical testing through accredited laboratories
>
> **By using this software, you acknowledge that:**
> - You understand this is for fun and learning only
> - You will not make any health decisions based on these results
> - The authors bear no responsibility for any actions taken based on the output
>
> **Genetic predisposition ≠ destiny.** Most traits and diseases are influenced by lifestyle, environment, and complex interactions between hundreds of genes.

---

## 🌟 Features

- **🔒 100% Private** - All analysis runs locally on your machine. Your genetic data never leaves your computer.
- **💰 Free & Open Source** - MIT License, no hidden costs, no subscriptions
- **🖥️ Interactive CLI** - Easy-to-use menu-driven interface
- **🌍 Bilingual** - English and Polish language support
- **📊 Comprehensive** - 500+ genetic variants analyzed across 16 categories
- **🎓 Educational** - Learn about your genetics in a fun, accessible way

---

## 📋 What Can You Explore?

| Category | Description | Examples |
|----------|-------------|----------|
| 🏃 **Athletic Genetics** | Muscle fiber type, endurance potential | ACTN3, ACE, PPARGC1A |
| 🏊 **Triathlon Predisposition** | VO2max, fat metabolism, injury risk | Sprint vs Ironman tendencies |
| 💊 **Personalized Insights** | Training style, supplement ideas, nutrition | Based on MTHFR, COMT, CYP1A2 |
| 🧬 **Comprehensive Analysis** | 500+ variants across 16 domains | Pharmacogenomics, traits |
| 🔬 **Variant Exploration** | Known significant variants | Research-grade exploration |
| 📊 **Polygenic Scores** | Combined variant effects | Educational risk exploration |
| 🎲 **Fun Genetics** | Traits and ancestry hints | Eye color, caffeine, cilantro |

---

## 🚀 Quick Start

### Option 1: Docker (Recommended) 🐳

The easiest way to run Helixight is with Docker - no installation required!

```bash
# Clone the repository
git clone https://github.com/helixight/helixight-oss.git
cd helixight-oss

# Start the web interface
docker-compose up -d

# Open http://localhost:8501 in your browser
```

That's it! The web interface will be available at **http://localhost:8501**

<p align="center">
  <img src="docs/screenshot.png" alt="Helixight Web Interface" width="600">
</p>

#### Docker Usage Tips

```bash
# View logs
docker-compose logs -f

# Stop the container
docker-compose down

# Rebuild after updates
docker-compose build --no-cache
docker-compose up -d
```

#### Adding Your Data

Place your VCF files in the `data/` directory, or upload them through the web interface:

```bash
cp your_variants.vcf.gz ./data/
```

**Supported file formats:**
- `.vcf` - Variant Call Format
- `.vcf.gz` - Compressed VCF (gzip)

Analysis results are saved to the `results/` directory.

---

### Option 2: Command Line Installation

```bash
# Clone the repository
git clone https://github.com/helixight/helixight-oss.git
cd helixight-oss

# Install dependencies (Ubuntu/Debian)
sudo apt install bcftools samtools tabix wget

# Or use the installer
./install.sh

# Run the interactive menu
./helixight.sh
```

---

## 📦 Requirements

### For Docker (Recommended)
- Docker Engine 20.10+
- Docker Compose v2.0+
- 4GB RAM minimum (8GB recommended for large VCF files)

### For CLI Installation

| Tool | Purpose | Install |
|------|---------|---------|
| bcftools | VCF manipulation | `apt install bcftools` |
| samtools | BAM/SAM processing | `apt install samtools` |
| tabix | VCF indexing | `apt install tabix` |
| wget | File downloads | `apt install wget` |

### One-Line Install

**Ubuntu/Debian:**
```bash
sudo apt update && sudo apt install -y bcftools samtools tabix wget gzip curl
```

**macOS (Homebrew):**
```bash
brew install bcftools samtools htslib wget
```

---

## 📖 Usage

### Interactive Mode (Recommended)

```bash
./helixight.sh
```

This launches an interactive menu where you can:
1. ✅ Install dependencies
2. ✅ Download reference genome (if you have BAM files)
3. ✅ Create VCF from BAM files
4. ✅ Run various genetic explorations
5. ✅ Generate summary reports

### Direct Script Execution

```bash
# Athletic genetics exploration
./scripts/athletic_genetics.sh your_variants.vcf.gz

# Triathlon predisposition analysis
./scripts/triathlon_genetics.sh your_variants.vcf.gz

# Personalized insights
./scripts/personalized_protocol.sh your_variants.vcf.gz

# Comprehensive exploration (500+ variants)
./scripts/mega_analysis.sh your_variants.vcf.gz

# Fun traits
./scripts/fun_genetics.sh your_variants.vcf.gz
```

---

## 🧬 Creating VCF from BAM

If you have a BAM file from whole genome sequencing:

### Option 1: Interactive Menu
```bash
./helixight.sh
# Select option 4: Create VCF from BAM file
```

### Option 2: Manual bcftools
```bash
bcftools mpileup -Ou -f GRCh38.fa your_sample.bam | \
bcftools call -mv -Oz -o your_variants.vcf.gz

bcftools index -t your_variants.vcf.gz
```

---

## 📊 Available Analyses

### 🏃 Athletic Genetics
Explore genes associated with athletic performance:
- **ACTN3** - Muscle fiber type tendencies
- **ACE** - Cardiovascular efficiency markers
- **PPARGC1A** - Mitochondrial biogenesis
- **IL6** - Recovery and inflammation markers
- **COMT** - Stress response (Warrior vs Worrier)

### 🏊 Triathlon Genetics
Comprehensive triathlon-oriented exploration:
- VO2max potential markers
- Fat oxidation capacity
- Lactate threshold genetics
- Injury susceptibility markers
- Mental resilience indicators
- Per-discipline scoring (Swim/Bike/Run)

### 💊 Personalized Insights
Educational recommendations based on genetics:
- Training style suggestions
- Supplement considerations (MTHFR, VDR, etc.)
- Nutritional tendencies
- Lifestyle insights

### 🧬 Mega Analysis
Comprehensive exploration of 500+ variants:
- Ancestral markers (Neanderthal DNA)
- Carrier status exploration
- Chronotype (owl vs lark)
- Sensory genetics (taste, smell)
- Trait predictions

---

## 📁 Project Structure

```
helixight-oss/
├── helixight.sh          # Main interactive CLI launcher
├── install.sh            # Quick installer for CLI
├── Dockerfile            # Docker container definition
├── docker-compose.yml    # Docker Compose configuration
├── README.md             # English documentation
├── README_PL.md          # Polish documentation
├── LICENSE               # MIT License
├── scripts/              # Analysis scripts (22 scripts)
│   ├── athletic_genetics.sh
│   ├── triathlon_genetics.sh
│   ├── personalized_protocol.sh
│   ├── mega_analysis.sh
│   ├── clinvar_scan.sh
│   ├── cancer_genes_scan.sh
│   └── ... (more scripts)
├── frontend/             # Web interface (Streamlit)
│   ├── app.py            # Main Streamlit application
│   ├── analysis.py       # Python wrapper for scripts
│   └── requirements.txt  # Python dependencies
├── data/                 # Place your VCF files here
├── results/              # Analysis results output
└── docs/
    └── variant_database.md
```

---

## 🔬 Scientific Note

The variants analyzed by this toolkit are based on published research and publicly available databases including:
- [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [PharmGKB](https://www.pharmgkb.org/)
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
- [SNPedia](https://www.snpedia.com/)

However, genetic science is complex and constantly evolving. Many associations have small effect sizes, limited replication, or population-specific effects. This tool provides a simplified educational view and should not be confused with clinical interpretation.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Ideas for Contributions
- [ ] Additional genetic variants with citations
- [ ] More analysis scripts
- [ ] Improved visualizations
- [x] Web interface ✅
- [x] Docker container ✅
- [ ] Additional language support
- [ ] PDF report generation
- [ ] Interactive result charts

---

## 📝 Changelog

### v1.0.1
- **Security**: Fixed command injection vulnerabilities in shell scripts
- **Security**: Added input validation for all file paths and user inputs
- **Security**: Implemented lock files to prevent race conditions in ClinVar downloads
- **Bug Fix**: Fixed file upload validation for `.vcf.gz` files in web interface
- **Bug Fix**: Fixed download button nesting issue in Streamlit
- **Bug Fix**: Fixed potential KeyError in results display
- **Bug Fix**: Fixed unreachable code in analysis cancellation
- **Bug Fix**: Added missing `RESULTS_DIR` environment variable in Docker
- **Improvement**: Added proper error handling with `set -euo pipefail` in all scripts
- **Improvement**: Replaced hardcoded values with configurable constants
- **Improvement**: Added Python 3.8+ compatibility for type hints

### v1.0.0
- Initial release with CLI and Docker web interface
- 22 analysis scripts covering athletic, health, and fun genetics
- Streamlit-based web frontend
- Docker containerization with bioinformatics tools

---

## 📜 License

MIT License - see [LICENSE](LICENSE) file.

**Additional Terms:** By using this software, you agree to the disclaimer above and acknowledge this is not for clinical or medical use.

---

## ⚠️ Final Reminder

> **This is a hobby project for genetic exploration and education.**
>
> 🎓 Use it to learn about genetics
> 🎮 Use it for fun and curiosity
> 🔬 Use it for personal research
>
> 🚫 Do NOT use it for medical decisions
> 🚫 Do NOT use it as a diagnostic tool
> 🚫 Do NOT skip professional genetic counseling if you have concerns

---

<p align="center">
  Made with 🧬 for the curious
</p>

<p align="center">
  <sub>Remember: Your genes are not your destiny. Lifestyle, environment, and choices matter too!</sub>
</p>
