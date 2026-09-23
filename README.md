# 🧬 Helixight - Open Source Genetic Analysis Toolkit

**Version 1.3.0**

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
- **🌍 Bilingual CLI menu** - English/Polish menu (`L` in the menu or `HELIXIGHT_LANGUAGE=pl`); analysis reports are currently a mix of Polish and English, and the web UI is English-only
- **📊 Comprehensive** - 13 analysis modules in the web UI; the CLI menu offers 12 (it adds the BAM-based STR and EVO2 tools but lacks pharmacogenomics, carrier screening and blood type, which you can run directly from `scripts/`); Mega Analysis alone covers ~80 variants in 10 sections
- **🎓 Educational** - Learn about your genetics in a fun, accessible way

---

## 📋 What Can You Explore?

| Category | Description | Examples |
|----------|-------------|----------|
| 🏃 **Athletic Genetics** | Muscle fiber type, endurance potential | ACTN3, ACE, PPARGC1A |
| 🏊 **Triathlon Predisposition** | VO2max, fat metabolism, injury risk | Sprint vs Ironman tendencies |
| 💊 **Personalized Insights** | Training style, supplement ideas, nutrition | Based on MTHFR, COMT, CYP1A2 |
| 💉 **Pharmacogenomics** | Drug metabolism, medication response | CYP2D6, CYP2C19, VKORC1 |
| 🧪 **Carrier Screening** | Recessive disease carrier status | CF, sickle cell, Tay-Sachs |
| 🩸 **Blood Type Prediction** | ABO and Rh factor from genetics | Blood type, transfusion compatibility |
| 🧬 **Comprehensive Analysis** | Mega Analysis: ~80 variants in 10 sections | Ancestry, carrier status, chronotype, senses, longevity |
| 🔬 **Variant Exploration** | Known significant variants | Research-grade exploration |
| 📊 **Polygenic Scores** | Combined variant effects | Educational risk exploration |
| 🎲 **Fun Genetics** | Traits and ancestry hints | Eye color, caffeine, cilantro |

---

## 🚀 Quick Start

### Option 1: Docker (Recommended) 🐳

The easiest way to run Helixight is with Docker - no installation required!

```bash
# Clone the repository
git clone https://github.com/zd87pl/dna-analysis.git
cd dna-analysis

# The container runs as an unprivileged user; run it as YOUR user so that it
# can write into the bind-mounted ./data and ./results (bash does not export
# UID, and a new terminal or `sudo docker compose` loses it). Storing the ids
# in .env (git-ignored) makes this permanent:
echo "UID=$(id -u)" >> .env; echo "GID=$(id -g)" >> .env

# Start the web interface
docker compose up -d

# Open http://localhost:8501 in your browser
```

That's it! The web interface will be available at **http://localhost:8501**

Try it right away with the bundled sample: `data/sample_genome.vcf.gz` appears in the
"Available VCF files" picker through the `./data` mount (the image also ships a copy in
`/data` for plain `docker run`).

#### Docker Usage Tips

```bash
# View logs
docker compose logs -f

# Stop the container
docker compose down

# Rebuild after updates
docker compose build --no-cache
docker compose up -d
```

- **Network exposure**: the app has no authentication, so `docker-compose.yml` binds the
  port to `127.0.0.1:8501` only. To reach it from another machine set
  `HELIXIGHT_BIND=0.0.0.0` (or a specific interface address) **and** put a reverse proxy
  with authentication and TLS in front of it; do not expose it to a network you do not
  trust.
- **File ownership**: `docker-compose.yml` runs the container as `${UID}:${GID}` and falls
  back to uid 10001 when those are not set. On Linux the bind-mounted `./data` and
  `./results` belong to you, so the container can only write to them when it runs as your
  uid: put `UID=...`/`GID=...` in `.env` next to `docker-compose.yml` (see the quick start)
  or run `export UID; export GID=$(id -g)` in every shell before `docker compose up`
  (a new terminal or `sudo docker compose up` loses the export). The symptom of the
  fallback is `Upload failed: [Errno 13] Permission denied` on upload and every analysis,
  including the bundled sample, ending with `Analysis failed: [Errno 13] Permission denied:
  /app/results/...`. A separate cause with the same symptom: if `./data` or `./results`
  were created by the previous root-running image, fix them once with
  `sudo chown -R "$(id -u):$(id -g)" data results`.
- **Result cache**: every analysis result (full report text and findings) is cached as
  plain JSON under `results/.cache/` so re-running on the same file is instant. Use the
  **Clear Cache** button in the sidebar (or delete `results/.cache/`) to remove it; the
  cache toggle in the sidebar disables it for new runs.
- **Upload size**: web uploads are limited to 1000 MB
  (`STREAMLIT_SERVER_MAX_UPLOAD_SIZE` in `docker-compose.yml`) because Streamlit holds
  the whole upload in memory. For larger files copy them into `./data` and pick them
  from "Available VCF files" instead.
- **Reference cache**: `clinvar_scan.sh` downloads ClinVar (~70 MB) once into
  `results/.reference/`; mount `./reference` there (see the commented volume in
  `docker-compose.yml`) if you want to keep it outside `results/`. The download is
  verified against NCBI's `.md5`; set `HELIXIGHT_SKIP_CHECKSUM=1` to proceed without it,
  and `HELIXIGHT_CLINVAR_MAX_AGE_DAYS` (default 90) to control when the cached copy is
  refreshed. The script needs `flock` (util-linux) and `md5sum` (coreutils) besides
  `bcftools` and `wget`; run from the CLI it caches under `$HELIXIGHT_REFERENCE_DIR`
  (default `./clinvar_analysis/`).
- **Result layout**: the web UI writes each run to
  `results/<vcf stem>_<hash12>/<analysis_id>/`, the CLI menu to
  `helixight_results/<input-name>/` (combined report in `helixight_results/`).
- **Environment variables** set by the web UI for the scripts: `HELIXIGHT_WORK_DIR`
  (`results/.work`, bgzip/index working copies), `HELIXIGHT_REFERENCE_DIR`
  (`results/.reference`, ClinVar cache) and `HELIXIGHT_ANALYSIS_TIMEOUT` (seconds per
  analysis, default 3600, read from the container environment).
- **Memory**: the compose file limits the container to 8 GB (`deploy.resources` for
  `docker compose`, `mem_limit` for the legacy `docker-compose` v1 binary).

#### Adding Your Data

Place your VCF files in the `data/` directory, or upload them through the web interface:

```bash
cp your_variants.vcf.gz ./data/
```

**Supported file formats:**
- `.vcf` - plain-text Variant Call Format
- `.vcf.gz` - gzip- or bgzip-compressed VCF (an index is not required)
- `chr`-prefixed contig names (`chr1`, `chrM`) are accepted

Every script converts its input to a sorted, bgzip-compressed, tabix-indexed working
copy with bare contig names before querying it (kept under `results/.work/` in the web UI,
and under `$HELIXIGHT_WORK_DIR` or `$TMPDIR/helixight-<uid>` (default `/tmp/helixight-<uid>`,
created mode 700) when a script is run directly). The copy, its index and the sort log are a
full copy of your genome and persist until you delete that directory. You
do not have to prepare anything. If you want
to do it yourself once, the equivalent command is:

```bash
bcftools sort -Oz -o your_variants.vcf.gz input.vcf && bcftools index -t your_variants.vcf.gz
```

**The VCF must be called against GRCh38.** GRCh37/hg19 files are refused with a clear
message (exit code 3) instead of being silently misread, because all coordinates in the
scripts are GRCh38.

Analysis results are saved to the `results/` directory.

---

### Option 2: Command Line Installation

```bash
# Clone the repository
git clone https://github.com/zd87pl/dna-analysis.git
cd dna-analysis

# Install dependencies (Ubuntu/Debian)
sudo apt install bcftools samtools tabix wget bc

# Or use the installer
./install.sh

# Run the interactive menu
./helixight.sh
```

---

## 📦 Requirements

### For Docker (Recommended)
- Docker Engine 20.10+
- Docker Compose v2.0+ (`docker compose`; the legacy `docker-compose` v1 binary works but ignores `deploy.resources`, which is why `mem_limit` is set as well)
- 4GB RAM minimum (8GB recommended for large VCF files)

### For CLI Installation

| Tool | Purpose | Install |
|------|---------|---------|
| bcftools | VCF manipulation | `apt install bcftools` |
| samtools | BAM/SAM processing | `apt install samtools` |
| tabix | VCF indexing (bgzip/tabix) | `apt install tabix` |
| wget | File downloads | `apt install wget` |
| bc | Score arithmetic in the scripts | `apt install bc` |
| flock, md5sum | ClinVar cache locking and checksum (`clinvar_scan.sh`) | util-linux / coreutils (preinstalled on Linux) |

### One-Line Install

**Ubuntu/Debian:**
```bash
sudo apt update && sudo apt install -y bcftools samtools tabix wget gzip curl bc
```

**macOS (Homebrew):**
```bash
brew install bcftools samtools htslib wget
```

`install.sh` uses these package managers; on macOS it requires [Homebrew](https://brew.sh) to be installed first (it stops with a message otherwise). macOS ships `bc`.

macOS note: `clinvar_scan.sh` (CLI options `20` and `41`) verifies the ClinVar download
with `md5sum` and serialises concurrent downloads with `flock`, neither of which macOS
ships. `brew install coreutils` provides `md5sum` (as `gmd5sum`, or `md5sum` with gnubin on
PATH); without it the script falls back to the built-in `md5`. Without `flock` the scan
still runs but two simultaneous ClinVar downloads are not serialised (a warning is
printed). The Docker image has both.

---

## 📖 Usage

### Interactive Mode (Recommended)

```bash
./helixight.sh
```

This launches an interactive menu with these options:
- `1` Install dependencies, `2` Check dependencies, `3` Download reference genome (if you have BAM files)
- `4` Create VCF from BAM files, `5` Select an existing VCF file
- `10`-`14` Genetic analyses (athletic, triathlon, personalized protocol, mega analysis, fun genetics)
- `20`-`23` Health & clinical (ClinVar scan, cancer genes, polygenic risk scores, rare variants)
- `30`-`32` Advanced (STR/repeat analysis and EVO2 data preparation need a BAM; mtDNA haplogroup)
- `40` Generate summary report (plain text; also converted to PDF if pandoc and a PDF engine such as LaTeX are installed - the guaranteed PDF export is in the web UI), `41` Run the 9 VCF analyses of options 10-14 and 20-23, then generate the summary report
- `L` Change language (EN/PL), `H` Help, `Q` Quit

### Direct Script Execution

```bash
# Athletic genetics exploration
./scripts/athletic_genetics.sh your_variants.vcf.gz

# Triathlon predisposition analysis
./scripts/triathlon_genetics.sh your_variants.vcf.gz

# Personalized insights
./scripts/personalized_protocol.sh your_variants.vcf.gz

# Comprehensive exploration (~80 variants in 10 sections)
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
Comprehensive exploration of ~80 variants in 10 sections:
- Ancestral markers (Neanderthal DNA)
- Carrier status exploration
- Chronotype (owl vs lark)
- Sensory genetics (taste, smell)
- Trait predictions

---

## 📁 Project Structure

```
dna-analysis/
├── helixight.sh          # Main interactive CLI launcher
├── install.sh            # Quick installer for CLI
├── Dockerfile            # Docker container definition
├── docker-compose.yml    # Docker Compose configuration
├── pytest.ini            # Pytest configuration
├── README.md             # English documentation
├── README_PL.md          # Polish documentation
├── LICENSE               # MIT License
├── .github/
│   ├── workflows/ci.yml  # CI: pytest, shellcheck, script smoke run, docker build
│   └── dependabot.yml    # Dependency update checks (pip, docker, actions)
├── scripts/              # Analysis scripts (bash + bcftools)
│   ├── lib/
│   │   ├── vcf_common.sh # Shared helpers (input conversion, genotype lookups)
│   │   └── GRCh38_COORDINATES.md  # Canonical GRCh38 coordinate table
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
│   ├── pdf_report.py     # PDF report generator
│   └── requirements.txt  # Python dependencies
├── tests/                # Test suite
│   ├── conftest.py       # Test environment (RESULTS_DIR, sys.path)
│   ├── test_analysis.py  # Analysis module tests
│   ├── test_pdf_report.py# PDF generation tests
│   ├── shell/test_vcf_common.sh  # Shared shell library tests
│   ├── shell/test_scripts.sh     # Per-script regression tests (synthetic VCFs)
│   └── requirements.txt  # Test dependencies
├── data/                 # Place your VCF files here
│   ├── sample_genome.vcf        # Sample test data (GRCh38, synthetic genotypes)
│   ├── sample_genome.vcf.gz     # Same, bgzip-compressed
│   └── sample_genome.vcf.gz.tbi # ... and its tabix index
├── results/              # Analysis results output (.cache/, .work/, .reference/ live here)
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

## 🛠️ Development / Running Tests

Python 3.11 is used in CI and in the Docker image.

```bash
# Python tests (frontend wrappers and PDF report)
python -m venv .venv && . .venv/bin/activate
pip install -r tests/requirements.txt
python -m pytest

# Shared shell library tests (need bcftools, bgzip, tabix)
bash tests/shell/test_vcf_common.sh

# Per-script regression tests (synthetic VCFs; need bcftools)
bash tests/shell/test_scripts.sh

# Lint the shell scripts
shellcheck -S error scripts/*.sh scripts/lib/*.sh helixight.sh install.sh
```

The scripts can be exercised without the web UI; run them from a scratch directory
(they write their output directories into the current working directory):

```bash
mkdir -p /tmp/hx && cd /tmp/hx
HELIXIGHT_WORK_DIR=/tmp/hx/work bash /path/to/dna-analysis/scripts/athletic_genetics.sh /path/to/dna-analysis/data/sample_genome.vcf.gz
```

Every script ends its report with `HX_DATA_COVERAGE: <found>/<queried> markers`; the web
UI uses it to tell "no usable data in this file" apart from a genuine result.

CI (`.github/workflows/ci.yml`) runs on every push to `main` and every pull request: the
pytest suite with coverage, `shellcheck`, a smoke run of every web-UI analysis script
against the sample VCF (bare and `chr`-prefixed contigs), and a Docker build.

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. Please update
`README.md` and `README_PL.md` together.

### Ideas for Contributions
- [ ] Additional genetic variants with citations
- [ ] More analysis scripts
- [ ] Improved visualizations
- [x] Web interface ✅
- [x] Docker container ✅
- [ ] Additional language support
- [x] PDF report generation ✅
- [x] VCF quality metrics ✅
- [x] Result caching ✅
- [x] Sample test data ✅
- [x] Interactive result charts ✅

---

## 📝 Changelog

### v1.3.0
- **Bug Fix**: Input normalisation - every script converts plain, gzip, bgzip and `chr`-prefixed VCFs to a sorted, indexed, bare-contig working copy (`scripts/lib/vcf_common.sh`) instead of failing or silently returning "no data"; GRCh37 files are refused with a clear message
- **Bug Fix**: Corrected GRCh38 coordinates and interpretations (COMT, CYP1A2, APOE, LCT, ALDH2, FTO, TCF7L2, MTHFR, VDR, HFE, CYP2D6, SLCO1B1, TPMT, GBA and others; canonical table in `scripts/lib/GRCh38_COORDINATES.md`); genotype lookups now verify REF/ALT before scoring
- **Bug Fix**: Scripts report `HX_DATA_COVERAGE` so the web UI distinguishes "no usable data" from a clean result; count bugs that aborted `pharmacogenomics.sh`, `carrier_screening.sh` and `clinvar_scan.sh` under `set -e` fixed
- **Bug Fix**: PDF export no longer crashes on emoji/Unicode in findings or on long lines
- **Bug Fix**: Upload validation in the web UI (file name sanitisation, VCF header check, no double in-memory copy); per-file working directories so concurrent runs do not overwrite each other
- **Bug Fix**: Sample VCF regenerated - parseable by bcftools, unique IDs/positions, GRCh38 coordinates, ACE I/D and ABO/CFTR deletions modelled, a benign BRCA2 record and six mtDNA (rCRS) records so the cancer and mtDNA modules complete on the demo; shipped bgzip-compressed with a tabix index
- **Security**: Container runs as an unprivileged user; web port bound to `127.0.0.1` by default; upload limit lowered to 1000 MB; ClinVar download verified against its `.md5`
- **Improvement**: CI workflow (pytest + coverage, shellcheck, script smoke run on the sample VCF, Docker build) and Dependabot
- **Improvement**: Honest feature claims - variant counts, menu numbering and language support in the documentation match the code
- **Bug Fix**: Blood type - rs8176746 (B allele) read in plus-strand orientation (G>T) so B/B and A/B genomes are no longer reported as "A or O"; A/O and A/A genomes in variants-only VCFs (no rs7853989 record) are called instead of "Unknown"
- **Bug Fix**: Unverified GRCh37-range probes no longer produce calls: BRCA1/2 named-mutation positions are informational notes (not counted as pathogenic), Tay-Sachs in `mega_analysis.sh` uses the same TATC-insertion probe as `carrier_screening.sh`, COL5A1 rs12722 and NFE2L2 are skipped, PPARA rs4253778 is scored in the same (endurance = G) orientation by both fitness scripts, HBB HbC/IVS-I-1 rows moved next to rs334; CYP2D6*2 (normal function) is informational, not a red finding
- **Bug Fix**: Shared library - VCFs whose records use undeclared INFO/FORMAT tags are accepted, SNP+indel mixed records keep their SNP call, working copies live in a per-user mode-700 directory, no bash-4-only expansions (macOS bash 3.2), and a `##contig` line without `length=` triggers the "cannot verify build" warning
- **Bug Fix**: Web UI - an empty (0-record) VCF cannot be analysed, cache entries written by 1.2.0 are ignored (format version), a timed-out or cancelled script that leaves a child holding the pipe is killed after 5 s, `output_files` lists only the current run's files

### v1.2.0
- **Feature**: Pharmacogenomics (PGx) analysis - drug metabolism genetics (CYP2D6, CYP2C19, CYP2C9, VKORC1, etc.)
- **Feature**: Carrier Screening - recessive disease carrier status (CF, sickle cell, Tay-Sachs, etc.)
- **Feature**: Blood Type Prediction - ABO and Rh factor from genetic variants
- **Improvement**: Expanded sample VCF with 70+ variants for testing all analyses
- **Improvement**: Now 14 analysis modules available

### v1.1.0
- **Feature**: PDF report generation with professional formatting and charts
- **Feature**: VCF quality metrics display (variant counts, Ti/Tv ratio, genotype distribution)
- **Feature**: Result caching to avoid re-running analyses on the same file
- **Feature**: Sample VCF test data (`data/sample_genome.vcf`) with 40+ variants
- **Feature**: Unit test suite with pytest for analysis and PDF modules
- **Improvement**: Cache toggle in settings UI
- **Improvement**: Clear cache button for manual cache management

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
