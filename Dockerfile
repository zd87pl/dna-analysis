#===============================================================================
# HELIXIGHT - DNA Analysis Toolkit
# Docker container with web frontend
#===============================================================================

FROM python:3.11-slim-bookworm

LABEL maintainer="Helixight <hello@helixight.com>"
LABEL description="Open-source genetic analysis toolkit with web interface"
LABEL version="1.0.1"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Install system dependencies and bioinformatics tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    # Build essentials
    build-essential \
    curl \
    wget \
    ca-certificates \
    # Compression tools
    gzip \
    bzip2 \
    xz-utils \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    # Bioinformatics tools
    bcftools \
    samtools \
    tabix \
    # Utilities
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
COPY docs/ /app/docs/

# Make scripts executable
RUN chmod +x /app/helixight.sh \
    && chmod +x /app/scripts/*.sh

# Create data and results directories for volume mounts
RUN mkdir -p /data /app/results

# Set environment variables
ENV HELIXIGHT_HOME=/app
ENV SCRIPTS_DIR=/app/scripts
ENV DATA_DIR=/data
ENV RESULTS_DIR=/app/results

# Expose Streamlit port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Run Streamlit
ENTRYPOINT ["streamlit", "run", "/app/frontend/app.py", \
    "--server.port=8501", \
    "--server.address=0.0.0.0", \
    "--server.headless=true", \
    "--browser.gatherUsageStats=false"]
