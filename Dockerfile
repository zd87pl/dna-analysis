#===============================================================================
# HELIXIGHT - DNA Analysis Toolkit
# Docker container with web frontend
#===============================================================================

FROM python:3.11-slim-bookworm

LABEL maintainer="Helixight <hello@helixight.com>"
LABEL description="Open-source genetic analysis toolkit with web interface"
LABEL version="1.3.0"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Install system dependencies and bioinformatics tools.
# Runtime packages only: nothing is compiled in this image (every Python
# dependency installs from a wheel), so no compilers or -dev headers.
RUN apt-get update && apt-get install -y --no-install-recommends \
    # curl: HEALTHCHECK; wget: ClinVar download in clinvar_scan.sh
    curl \
    wget \
    ca-certificates \
    # Compression tools
    gzip \
    bzip2 \
    xz-utils \
    # Bioinformatics tools
    bcftools \
    samtools \
    tabix \
    # Utilities (bc: score arithmetic in the scripts; procps: ps for the UI)
    bc \
    procps \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Create app directory
WORKDIR /app

# Copy requirements first for better caching
COPY frontend/requirements.txt /app/frontend/requirements.txt

# Install Python dependencies
RUN pip install --no-cache-dir -r /app/frontend/requirements.txt

# Copy application code
COPY scripts/ /app/scripts/
COPY frontend/ /app/frontend/
COPY helixight.sh /app/

# Make scripts executable
RUN chmod +x /app/helixight.sh \
    && chmod +x /app/scripts/*.sh

# Create data and results directories for volume mounts
RUN mkdir -p /data /app/results

# Bundled sample so a plain `docker run` has something to analyse
# (docker compose bind-mounts ./data over it, which also contains the sample)
COPY data/sample_genome.vcf.gz data/sample_genome.vcf.gz.tbi /data/

# Run as an unprivileged user: the app accepts uploads and shells out to bash.
# /data and /app/results are usually bind mounts; docker-compose.yml overrides
# the uid/gid with the invoking user's so host files stay owned by them.
RUN useradd -r -u 10001 -m helixight \
    && chown -R helixight:helixight /app /data

# Set environment variables
ENV SCRIPTS_DIR=/app/scripts
ENV DATA_DIR=/data
ENV RESULTS_DIR=/app/results

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

USER helixight

# Run Streamlit
ENTRYPOINT ["streamlit", "run", "/app/frontend/app.py", \
    "--server.port=8501", \
    "--server.address=0.0.0.0", \
    "--server.headless=true", \
    "--browser.gatherUsageStats=false"]
