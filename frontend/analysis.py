"""
Helixight Analysis Module
Python wrapper for bash analysis scripts with progress tracking
"""

import subprocess
import os
import re
import json
import threading
import hashlib
import gzip
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Callable, Dict, List, Any, Tuple
from enum import Enum
from datetime import datetime


class AnalysisStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


@dataclass
class VCFMetrics:
    """Stores VCF quality metrics"""
    total_variants: int = 0
    snps: int = 0
    indels: int = 0
    heterozygous: int = 0
    homozygous_alt: int = 0
    homozygous_ref: int = 0
    transitions: int = 0
    transversions: int = 0
    titv_ratio: float = 0.0
    quality_pass: int = 0
    file_size_mb: float = 0.0


@dataclass
class AnalysisResult:
    """Stores results from an analysis script"""
    name: str
    status: AnalysisStatus
    output: str = ""
    error: str = ""
    scores: Dict[str, Any] = field(default_factory=dict)
    findings: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    output_files: List[str] = field(default_factory=list)


class AnalysisRunner:
    """Runs Helixight analysis scripts and parses their output"""

    SCRIPTS_DIR = os.environ.get("SCRIPTS_DIR", "/app/scripts")
    DATA_DIR = os.environ.get("DATA_DIR", "/data")
    RESULTS_DIR = os.environ.get("RESULTS_DIR", "/app/results")

    # Available analyses with metadata
    ANALYSES = {
        "athletic_genetics": {
            "script": "athletic_genetics.sh",
            "name": "Athletic & Fitness Genetics",
            "description": "Analyzes power, endurance, recovery, and injury risk genetics",
            "icon": "🏃",
            "category": "fitness"
        },
        "triathlon_genetics": {
            "script": "triathlon_genetics.sh",
            "name": "Triathlon Predisposition",
            "description": "Sport-specific analysis for triathlon performance",
            "icon": "🏊",
            "category": "fitness"
        },
        "personalized_protocol": {
            "script": "personalized_protocol.sh",
            "name": "Personalized Protocol",
            "description": "Training, supplement, and lifestyle recommendations",
            "icon": "💊",
            "category": "fitness"
        },
        "mega_analysis": {
            "script": "mega_analysis.sh",
            "name": "Mega Analysis (500+ variants)",
            "description": "Comprehensive analysis across 16 health domains",
            "icon": "🧬",
            "category": "comprehensive"
        },
        "fun_genetics": {
            "script": "fun_genetics.sh",
            "name": "Fun Genetics",
            "description": "Traits, ancestry hints, and interesting findings",
            "icon": "🎲",
            "category": "fun"
        },
        "clinvar_scan": {
            "script": "clinvar_scan.sh",
            "name": "ClinVar Pathogenic Scan",
            "description": "Screens for known pathogenic variants",
            "icon": "🔬",
            "category": "health"
        },
        "cancer_genes_scan": {
            "script": "cancer_genes_scan.sh",
            "name": "Cancer Genes Analysis",
            "description": "BRCA1/2, TP53, Lynch syndrome analysis",
            "icon": "🎀",
            "category": "health"
        },
        "polygenic_risk_scores": {
            "script": "polygenic_risk_scores.sh",
            "name": "Polygenic Risk Scores",
            "description": "Combined risk scores for complex diseases",
            "icon": "📊",
            "category": "health"
        },
        "rare_variants": {
            "script": "rare_variants.sh",
            "name": "Rare Variants Analysis",
            "description": "Identification of rare genetic variants",
            "icon": "💎",
            "category": "health"
        },
        "str_analysis": {
            "script": "str_analysis.sh",
            "name": "STR/Repeat Analysis",
            "description": "Short tandem repeat expansion analysis",
            "icon": "🔁",
            "category": "advanced"
        },
        "mtdna_haplogroup": {
            "script": "extract_mtdna_haplogrep.sh",
            "name": "mtDNA Haplogroup",
            "description": "Mitochondrial DNA ancestry analysis",
            "icon": "🌍",
            "category": "advanced"
        }
    }

    def __init__(self, vcf_path: str, use_cache: bool = True):
        self.vcf_path = vcf_path
        self.results: Dict[str, AnalysisResult] = {}
        self._progress_callback: Optional[Callable] = None
        self._cancel_flag = threading.Event()
        self.use_cache = use_cache

    def set_progress_callback(self, callback: Callable[[str, float, str], None]):
        """Set callback for progress updates: (analysis_name, progress, message)"""
        self._progress_callback = callback

    def cancel(self):
        """Cancel running analyses"""
        self._cancel_flag.set()

    def run_analysis(self, analysis_id: str, force_rerun: bool = False) -> AnalysisResult:
        """Run a single analysis script"""
        if analysis_id not in self.ANALYSES:
            return AnalysisResult(
                name=analysis_id,
                status=AnalysisStatus.FAILED,
                error=f"Unknown analysis: {analysis_id}"
            )

        analysis = self.ANALYSES[analysis_id]

        # Check cache first (unless force_rerun is True)
        if self.use_cache and not force_rerun:
            cached_result = _result_cache.get(self.vcf_path, analysis_id)
            if cached_result is not None:
                if self._progress_callback:
                    self._progress_callback(analysis_id, 1.0, f"Loaded from cache: {analysis['name']}")
                self.results[analysis_id] = cached_result
                return cached_result

        script_path = os.path.join(self.SCRIPTS_DIR, analysis["script"])

        if not os.path.exists(script_path):
            return AnalysisResult(
                name=analysis["name"],
                status=AnalysisStatus.FAILED,
                error=f"Script not found: {script_path}"
            )

        result = AnalysisResult(
            name=analysis["name"],
            status=AnalysisStatus.RUNNING
        )

        if self._progress_callback:
            self._progress_callback(analysis_id, 0.0, f"Starting {analysis['name']}...")

        try:
            # Run the script
            process = subprocess.Popen(
                ["/bin/bash", script_path, self.vcf_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=self.RESULTS_DIR
            )

            # Check cancel flag during execution with timeout polling
            try:
                stdout, stderr = process.communicate(timeout=1)
            except subprocess.TimeoutExpired:
                if self._cancel_flag.is_set():
                    process.kill()
                    process.wait()
                    result.status = AnalysisStatus.FAILED
                    result.error = "Analysis cancelled"
                    return result
                # Continue waiting if not cancelled
                stdout, stderr = process.communicate()

            # Final check after completion
            if self._cancel_flag.is_set():
                result.status = AnalysisStatus.FAILED
                result.error = "Analysis cancelled"
                return result

            result.output = stdout
            result.error = stderr

            if process.returncode == 0:
                result.status = AnalysisStatus.COMPLETED
                # Parse output for scores and findings
                self._parse_output(result, stdout)
                # Cache successful results
                if self.use_cache:
                    _result_cache.set(self.vcf_path, analysis_id, result)
            else:
                result.status = AnalysisStatus.FAILED

            if self._progress_callback:
                self._progress_callback(analysis_id, 1.0, f"Completed {analysis['name']}")

        except Exception as e:
            result.status = AnalysisStatus.FAILED
            result.error = str(e)

        self.results[analysis_id] = result
        return result

    def run_multiple(self, analysis_ids: List[str]) -> Dict[str, AnalysisResult]:
        """Run multiple analyses sequentially"""
        results = {}
        total = len(analysis_ids)

        for i, analysis_id in enumerate(analysis_ids):
            if self._cancel_flag.is_set():
                break

            if self._progress_callback:
                progress = i / total
                self._progress_callback(
                    "overall",
                    progress,
                    f"Running {i+1}/{total}: {self.ANALYSES.get(analysis_id, {}).get('name', analysis_id)}"
                )

            results[analysis_id] = self.run_analysis(analysis_id)

        return results

    def _parse_output(self, result: AnalysisResult, output: str):
        """Parse script output for scores, findings, and warnings"""
        lines = output.split('\n')

        for line in lines:
            # Parse score patterns like "POWER SCORE: 4 / 6" or "Power: 80%"
            score_match = re.search(r'([\w\s]+)(?:SCORE)?:\s*(\d+)\s*/\s*(\d+)', line, re.IGNORECASE)
            if score_match:
                name = score_match.group(1).strip()
                score = int(score_match.group(2))
                max_score = int(score_match.group(3))
                result.scores[name] = {
                    "value": score,
                    "max": max_score,
                    "percentage": round(score / max_score * 100) if max_score > 0 else 0
                }
                continue

            # Parse percentage patterns
            pct_match = re.search(r'([\w\s]+):\s*(\d+)%', line)
            if pct_match:
                name = pct_match.group(1).strip()
                if name not in result.scores:
                    result.scores[name] = {"percentage": int(pct_match.group(2))}

            # Parse findings (lines with specific markers)
            if any(marker in line for marker in ['🔴', '⚠️', 'FOUND', 'Pathogenic']):
                result.findings.append(line.strip())

            # Parse warnings
            if any(marker in line for marker in ['Warning', 'UWAGA', '⚠️']):
                result.warnings.append(line.strip())

            # Parse positive findings
            if '✅' in line or 'No pathogenic' in line.lower():
                result.findings.append(line.strip())

    @classmethod
    def get_analyses_by_category(cls) -> Dict[str, List[Dict]]:
        """Group analyses by category"""
        categories = {}
        for aid, info in cls.ANALYSES.items():
            cat = info.get("category", "other")
            if cat not in categories:
                categories[cat] = []
            categories[cat].append({"id": aid, **info})
        return categories

    @classmethod
    def validate_vcf(cls, file_path: str) -> Tuple[bool, str]:
        """Validate a VCF file"""
        if not os.path.exists(file_path):
            return False, "File not found"

        try:
            # Check file header
            if file_path.endswith('.gz'):
                import gzip
                with gzip.open(file_path, 'rt') as f:
                    first_line = f.readline()
            else:
                with open(file_path, 'r') as f:
                    first_line = f.readline()

            if not first_line.startswith('##fileformat=VCF'):
                return False, "Not a valid VCF file (missing header)"

            # Get variant count
            result = subprocess.run(
                ['bcftools', 'view', '-H', file_path],
                capture_output=True,
                text=True
            )
            variant_count = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0

            return True, f"Valid VCF with {variant_count:,} variants"

        except Exception as e:
            return False, f"Validation error: {str(e)}"


def get_system_info() -> Dict[str, str]:
    """Get system information for the about page"""
    info = {}

    # bcftools version
    try:
        result = subprocess.run(['bcftools', '--version'], capture_output=True, text=True)
        info['bcftools'] = result.stdout.split('\n')[0] if result.returncode == 0 else 'Not found'
    except (OSError, subprocess.SubprocessError):
        info['bcftools'] = 'Not found'

    # samtools version
    try:
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True)
        info['samtools'] = result.stdout.split('\n')[0] if result.returncode == 0 else 'Not found'
    except (OSError, subprocess.SubprocessError):
        info['samtools'] = 'Not found'

    return info


def get_vcf_metrics(file_path: str) -> VCFMetrics:
    """
    Calculate quality metrics for a VCF file using bcftools stats

    Args:
        file_path: Path to VCF file

    Returns:
        VCFMetrics dataclass with quality statistics
    """
    metrics = VCFMetrics()

    if not os.path.exists(file_path):
        return metrics

    # Get file size
    metrics.file_size_mb = os.path.getsize(file_path) / (1024 * 1024)

    try:
        # Run bcftools stats
        result = subprocess.run(
            ['bcftools', 'stats', file_path],
            capture_output=True,
            text=True,
            timeout=60
        )

        if result.returncode != 0:
            return metrics

        # Parse stats output
        for line in result.stdout.split('\n'):
            if line.startswith('SN\t'):
                parts = line.split('\t')
                if len(parts) >= 4:
                    key = parts[2].rstrip(':')
                    try:
                        value = int(parts[3])
                    except ValueError:
                        continue

                    if 'number of records' in key:
                        metrics.total_variants = value
                    elif 'number of SNPs' in key:
                        metrics.snps = value
                    elif 'number of indels' in key:
                        metrics.indels = value

            # Parse transitions/transversions
            elif line.startswith('TSTV\t'):
                parts = line.split('\t')
                if len(parts) >= 5:
                    try:
                        metrics.transitions = int(parts[2])
                        metrics.transversions = int(parts[3])
                        metrics.titv_ratio = float(parts[4])
                    except (ValueError, IndexError):
                        pass

        # Count genotypes using bcftools query
        gt_result = subprocess.run(
            ['bcftools', 'query', '-f', '%GT\n', file_path],
            capture_output=True,
            text=True,
            timeout=60
        )

        if gt_result.returncode == 0:
            for gt in gt_result.stdout.strip().split('\n'):
                if gt in ('0/1', '1/0', '0|1', '1|0'):
                    metrics.heterozygous += 1
                elif gt in ('1/1', '1|1'):
                    metrics.homozygous_alt += 1
                elif gt in ('0/0', '0|0'):
                    metrics.homozygous_ref += 1

        # Count PASS variants
        pass_result = subprocess.run(
            ['bcftools', 'view', '-f', 'PASS', '-H', file_path],
            capture_output=True,
            text=True,
            timeout=60
        )

        if pass_result.returncode == 0:
            metrics.quality_pass = len(pass_result.stdout.strip().split('\n')) if pass_result.stdout.strip() else 0

    except (subprocess.TimeoutExpired, subprocess.SubprocessError, OSError):
        pass

    return metrics


def get_file_hash(file_path: str) -> str:
    """
    Calculate MD5 hash of a file for caching purposes

    Args:
        file_path: Path to file

    Returns:
        MD5 hash string
    """
    hash_md5 = hashlib.md5()

    try:
        # For large files, only hash first and last 1MB + file size
        file_size = os.path.getsize(file_path)
        chunk_size = 1024 * 1024  # 1MB

        with open(file_path, 'rb') as f:
            # Hash first chunk
            hash_md5.update(f.read(chunk_size))

            # Hash last chunk if file is larger
            if file_size > chunk_size * 2:
                f.seek(-chunk_size, 2)
                hash_md5.update(f.read(chunk_size))

        # Include file size in hash
        hash_md5.update(str(file_size).encode())

        return hash_md5.hexdigest()

    except (OSError, IOError):
        return ""


class ResultCache:
    """
    Simple file-based cache for analysis results
    """

    def __init__(self, cache_dir: str = None):
        self.cache_dir = cache_dir or os.path.join(
            os.environ.get("RESULTS_DIR", "/app/results"),
            ".cache"
        )
        os.makedirs(self.cache_dir, exist_ok=True)

    def _get_cache_path(self, vcf_hash: str, analysis_id: str) -> str:
        """Get path for cached result"""
        return os.path.join(self.cache_dir, f"{vcf_hash}_{analysis_id}.json")

    def get(self, vcf_path: str, analysis_id: str) -> Optional[AnalysisResult]:
        """
        Retrieve cached result if available

        Args:
            vcf_path: Path to VCF file
            analysis_id: Analysis identifier

        Returns:
            Cached AnalysisResult or None
        """
        vcf_hash = get_file_hash(vcf_path)
        if not vcf_hash:
            return None

        cache_path = self._get_cache_path(vcf_hash, analysis_id)

        if not os.path.exists(cache_path):
            return None

        try:
            with open(cache_path, 'r') as f:
                data = json.load(f)

            return AnalysisResult(
                name=data['name'],
                status=AnalysisStatus(data['status']),
                output=data.get('output', ''),
                error=data.get('error', ''),
                scores=data.get('scores', {}),
                findings=data.get('findings', []),
                warnings=data.get('warnings', []),
                output_files=data.get('output_files', [])
            )

        except (json.JSONDecodeError, KeyError, OSError):
            return None

    def set(self, vcf_path: str, analysis_id: str, result: AnalysisResult) -> bool:
        """
        Cache an analysis result

        Args:
            vcf_path: Path to VCF file
            analysis_id: Analysis identifier
            result: AnalysisResult to cache

        Returns:
            True if cached successfully
        """
        vcf_hash = get_file_hash(vcf_path)
        if not vcf_hash:
            return False

        cache_path = self._get_cache_path(vcf_hash, analysis_id)

        try:
            data = {
                'name': result.name,
                'status': result.status.value,
                'output': result.output,
                'error': result.error,
                'scores': result.scores,
                'findings': result.findings,
                'warnings': result.warnings,
                'output_files': result.output_files,
                'cached_at': datetime.now().isoformat() if 'datetime' in dir() else None
            }

            with open(cache_path, 'w') as f:
                json.dump(data, f)

            return True

        except (OSError, TypeError):
            return False

    def clear(self, vcf_path: str = None):
        """
        Clear cache entries

        Args:
            vcf_path: If provided, only clear cache for this file
        """
        if vcf_path:
            vcf_hash = get_file_hash(vcf_path)
            if vcf_hash:
                for f in os.listdir(self.cache_dir):
                    if f.startswith(vcf_hash):
                        try:
                            os.remove(os.path.join(self.cache_dir, f))
                        except OSError:
                            pass
        else:
            # Clear all cache
            for f in os.listdir(self.cache_dir):
                if f.endswith('.json'):
                    try:
                        os.remove(os.path.join(self.cache_dir, f))
                    except OSError:
                        pass


# Global cache instance
_result_cache = ResultCache()
