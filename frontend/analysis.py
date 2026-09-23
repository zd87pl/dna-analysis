"""
Helixight Analysis Module
Python wrapper for bash analysis scripts with progress tracking
"""

import subprocess
import os
import re
import json
import signal
import threading
import time
import hashlib
import gzip
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Callable, Dict, List, Any, Tuple
from enum import Enum
from datetime import datetime
from uuid import uuid4

# Wall-clock limit for a single analysis script (seconds)
ANALYSIS_TIMEOUT = int(os.environ.get("HELIXIGHT_ANALYSIS_TIMEOUT", "3600"))

# Line every analysis script prints once at the end (see scripts/lib/vcf_common.sh)
_COVERAGE_RE = re.compile(r'HX_DATA_COVERAGE:\s*(\d+)\s*/\s*(\d+)\s*markers')
_ANSI_RE = re.compile(r'\x1b\[[0-9;]*m')


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
    error: str = ""


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
    returncode: Optional[int] = None
    stderr: str = ""
    # HX_DATA_COVERAGE from the script: markers with a usable genotype / markers queried
    markers_found: Optional[int] = None
    markers_queried: Optional[int] = None
    from_cache: bool = False  # runtime only, not persisted


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
            "name": "Mega Analysis (~80 variants, 10 sections)",
            "description": "Ancestry, carrier status, chronotype, senses, longevity and more",
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
        "mtdna_haplogroup": {
            "script": "extract_mtdna_haplogrep.sh",
            "name": "mtDNA Haplogroup",
            "description": "Mitochondrial DNA ancestry analysis",
            "icon": "🌍",
            "category": "advanced"
        },
        "pharmacogenomics": {
            "script": "pharmacogenomics.sh",
            "name": "Pharmacogenomics (PGx)",
            "description": "Drug metabolism and medication response genetics",
            "icon": "💉",
            "category": "health"
        },
        "carrier_screening": {
            "script": "carrier_screening.sh",
            "name": "Carrier Screening",
            "description": "Recessive disease carrier status analysis",
            "icon": "🧪",
            "category": "health",
            # Screens rare alleles: a variants-only VCF of a non-carrier has no
            # record at any of them, so 0/N found is the (negative) result, not
            # 'no usable data'.
            "absence_is_result": True
        },
        "blood_type": {
            "script": "blood_type.sh",
            "name": "Blood Type Prediction",
            "description": "ABO and Rh blood type prediction from genetics",
            "icon": "🩸",
            "absence_is_result": True,
            "category": "fun"
        }
    }

    def __init__(self, vcf_path: str, use_cache: bool = True):
        self.vcf_path = vcf_path
        self.results: Dict[str, AnalysisResult] = {}
        self._progress_callback: Optional[Callable] = None
        self._cancel_flag = threading.Event()
        self.use_cache = use_cache
        # Hash the VCF once per runner; ResultCache reuses it for every analysis
        self._vcf_hash = get_file_hash(vcf_path)
        # Lazily computed by _vcf_has_records(): None = not checked yet / unknown
        self._has_records: Optional[bool] = None

    def _vcf_has_records(self) -> Optional[bool]:
        """True/False when the VCF has at least one / no variant record; None if
        this could not be determined (bcftools missing or unreadable file).

        Reads only the first record so a multi-GB VCF is not streamed."""
        if self._has_records is None:
            self._has_records = _vcf_has_records(self.vcf_path)
        return self._has_records

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
            cached_result = get_result_cache().get(self.vcf_path, analysis_id, vcf_hash=self._vcf_hash)
            if cached_result is not None:
                cached_result.from_cache = True
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
            # Per-VCF, per-analysis working directory so concurrent runs on
            # different files never overwrite each other's output files.
            run_dir = os.path.join(
                self.RESULTS_DIR,
                f"{Path(self.vcf_path).stem}_{(self._vcf_hash or uuid4().hex)[:12]}",
                analysis_id
            )
            os.makedirs(run_dir, exist_ok=True)
            # Snapshot what is already there (earlier runs of the same VCF and
            # analysis) so output_files lists only what THIS run wrote.
            before = _snapshot_files(run_dir)

            env = {
                **os.environ,
                # Shared bgzip/index working copies (scripts/lib/vcf_common.sh)
                "HELIXIGHT_WORK_DIR": os.path.join(self.RESULTS_DIR, ".work"),
                # clinvar_scan.sh caches its ClinVar download here
                "HELIXIGHT_REFERENCE_DIR": os.path.join(self.RESULTS_DIR, ".reference"),
            }

            # Run the script (own process group so cancel/timeout kills children too)
            process = subprocess.Popen(
                ["/bin/bash", script_path, self.vcf_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=run_dir,
                env=env,
                start_new_session=True
            )

            # Wait, checking the cancel flag and the wall-clock limit every second
            deadline = time.monotonic() + ANALYSIS_TIMEOUT
            while True:
                try:
                    stdout, stderr = process.communicate(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    cancelled = self._cancel_flag.is_set()
                    if not cancelled and time.monotonic() <= deadline:
                        continue
                    _terminate_process(process)
                    try:
                        stdout, stderr = process.communicate(timeout=5)
                    except subprocess.TimeoutExpired:
                        # bash is gone but a grandchild that ignored SIGTERM still
                        # holds the stdout/stderr pipe: kill the whole group.
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except (ProcessLookupError, PermissionError, OSError):
                            process.kill()
                        stdout, stderr = process.communicate()
                    result.output = _ANSI_RE.sub('', stdout or "")
                    result.stderr = _ANSI_RE.sub('', stderr or "")
                    result.returncode = process.returncode
                    result.status = AnalysisStatus.FAILED
                    result.error = (
                        "Analysis cancelled" if cancelled
                        else f"Timed out after {ANALYSIS_TIMEOUT}s"
                    )
                    self.results[analysis_id] = result
                    return result

            stderr = _ANSI_RE.sub('', stderr or "")
            # Scripts colour some stdout lines (echo -e "${GREEN}..."); strip the
            # escapes so findings, the cache, the UI and the PDF get plain text.
            result.output = _ANSI_RE.sub('', stdout or "")
            result.stderr = stderr
            result.returncode = process.returncode

            # Collect every file the script wrote (new or modified) under its
            # working directory
            after = _snapshot_files(run_dir)
            result.output_files = sorted(
                path for path, mtime in after.items() if before.get(path) != mtime
            )

            if process.returncode == 0:
                result.status = AnalysisStatus.COMPLETED
                # Parse output for scores and findings
                self._parse_output(result, result.output)
                if result.markers_queried and result.markers_found == 0:
                    # The absence carve-out only holds for a VCF that has records
                    # elsewhere; a file with no variant records at all is no
                    # evidence of anything.
                    if analysis.get("absence_is_result") and self._vcf_has_records() is not False:
                        # Rare-allele screens: no record at any queried site is the
                        # expected negative result for a variants-only VCF. Keep it,
                        # say what it rests on, and do not cache it (a re-index
                        # with more data should re-run).
                        result.warnings.insert(0, (
                            f"⚠️ None of the {result.markers_queried} markers this analysis "
                            "looks for have a record in this VCF; the result is based on "
                            "their absence only"
                        ))
                    else:
                        # Exit 0 but the script found none of the markers it looked
                        # for: the file carries no usable data for this analysis, so
                        # do not present (or cache) the default 'no variant' text.
                        result.status = AnalysisStatus.FAILED
                        result.error = (
                            "No usable variant data found for this analysis "
                            f"(0/{result.markers_queried} markers present in the VCF"
                            + ("; the file contains no variant records)"
                               if self._has_records is False else ")")
                        )
                elif self.use_cache:
                    # Cache successful results
                    get_result_cache().set(self.vcf_path, analysis_id, result, vcf_hash=self._vcf_hash)
            else:
                result.status = AnalysisStatus.FAILED
                result.error = stderr.strip() or (
                    f"Script exited with code {process.returncode} "
                    "(no stderr output; see Raw Output)"
                )

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
            # Data coverage line printed once by scripts/lib/vcf_common.sh
            cov_match = _COVERAGE_RE.search(line)
            if cov_match:
                result.markers_found = int(cov_match.group(1))
                result.markers_queried = int(cov_match.group(2))
                continue

            # '⚠' matches the marker with or without the U+FE0F presentation selector
            is_marker_line = any(m in line for m in ('🔴', '⚠', '✅'))

            if not is_marker_line:
                # Parse score patterns like "POWER SCORE: 4 / 6" or "Wynik CAD PRS: 5 / 20".
                # A keyword (SCORE/WYNIK, or the triathlon discipline labels
                # PŁYWANIE/KOLARSTWO/BIEGANIE from triathlon_genetics.sh) is
                # mandatory so genotype lines such as "BRCA2 4075delGT: 0/1" are
                # never mistaken for scores.
                score_match = re.search(
                    r'([\w ]*?(?:SCORE|WYNIK|PŁYWANIE|KOLARSTWO|BIEGANIE)[\w ]*?):\s*(\d+)\s*/\s*(\d+)\b',
                    line, re.IGNORECASE
                )
                if score_match:
                    name = score_match.group(1).strip()
                    score = int(score_match.group(2))
                    max_score = int(score_match.group(3))
                    if name and name not in result.scores:
                        pct = round(score / max_score * 100) if max_score > 0 else 0
                        result.scores[name] = {
                            "value": score,
                            "max": max_score,
                            "percentage": max(0, min(100, pct))
                        }
                    continue

                # Parse percentage patterns like "Power: 80%"
                pct_match = re.search(r'([^:\n]+?):\s*(\d+)%', line)
                if pct_match:
                    name = pct_match.group(1).strip()
                    if name and name not in result.scores:
                        result.scores[name] = {"percentage": max(0, min(100, int(pct_match.group(2))))}

            # Classify findings / warnings (each line lands in at most one list)
            low = line.lower()
            text = line.strip()
            if '🔴' in line or re.search(r'\bFOUND\b', line):
                # Real hits only
                result.findings.append(text)
            elif '⚠' in line or 'warning' in low or 'uwaga' in low:
                # Warnings, minus boilerplate disclaimers
                if not re.search(r'disclaimer|not for clinical use|educational', low):
                    result.warnings.append(text)
            elif '✅' in line or 'no pathogenic' in low:
                # Positive findings
                result.findings.append(text)

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
                with gzip.open(file_path, 'rt') as f:
                    first_line = f.readline()
            else:
                with open(file_path, 'r') as f:
                    first_line = f.readline()

            if not first_line.startswith('##fileformat=VCF'):
                return False, "Not a valid VCF file (missing header)"

            # Get variant count (no timeout: bcftools streams the whole file)
            result = subprocess.run(
                ['bcftools', 'view', '-H', file_path],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                # The first stderr line carries the diagnostic; the last is generic
                lines = result.stderr.strip().splitlines()
                return False, f"bcftools cannot parse this file: {lines[0] if lines else 'unknown error'}"

            variant_count = len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0
            if variant_count == 0:
                return True, "Valid VCF but it contains 0 variant records"

            return True, f"Valid VCF with {variant_count:,} variants"

        except Exception as e:
            return False, f"Validation error: {str(e)}"


def get_system_info() -> Dict[str, str]:
    """Get system information for the about page"""
    info = {}

    # errors='replace': some builds (e.g. Ubuntu's samtools) embed non-UTF-8
    # bytes in their --version banner, which must not abort the page render.
    # bcftools version
    try:
        result = subprocess.run(['bcftools', '--version'], capture_output=True, text=True,
                                errors='replace', timeout=10)
        info['bcftools'] = result.stdout.split('\n')[0] if result.returncode == 0 else 'Not found'
    except (OSError, subprocess.SubprocessError):
        info['bcftools'] = 'Not found'

    # samtools version
    try:
        result = subprocess.run(['samtools', '--version'], capture_output=True, text=True,
                                errors='replace', timeout=10)
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
            metrics.error = result.stderr.strip() or f"bcftools stats exited with code {result.returncode}"
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

        # bcftools stats exits 0 even when it could parse no record at all;
        # its stderr then explains why (e.g. malformed body lines)
        if metrics.total_variants == 0 and result.stderr.strip():
            metrics.error = result.stderr.strip()

        # Count genotypes using bcftools query (FORMAT fields need [..] brackets)
        gt_result = subprocess.run(
            ['bcftools', 'query', '-f', '[%GT\n]', file_path],
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
        elif not metrics.error:
            metrics.error = gt_result.stderr.strip() or f"bcftools query exited with code {gt_result.returncode}"

        # Count PASS variants
        pass_result = subprocess.run(
            ['bcftools', 'view', '-f', 'PASS', '-H', file_path],
            capture_output=True,
            text=True,
            timeout=60
        )

        if pass_result.returncode == 0:
            metrics.quality_pass = len(pass_result.stdout.strip().split('\n')) if pass_result.stdout.strip() else 0

    except (subprocess.TimeoutExpired, subprocess.SubprocessError, OSError) as e:
        if not metrics.error:
            metrics.error = str(e)

    return metrics


def _snapshot_files(directory: str) -> Dict[str, int]:
    """{absolute path: st_mtime_ns} of every file under directory"""
    snapshot: Dict[str, int] = {}
    for root, _dirs, files in os.walk(directory):
        for fname in files:
            path = os.path.join(root, fname)
            try:
                snapshot[path] = os.stat(path).st_mtime_ns
            except OSError:
                continue
    return snapshot


def _vcf_has_records(file_path: str) -> Optional[bool]:
    """True/False when the VCF has at least one / no variant record, None when
    this cannot be determined. Only the first record is read."""
    try:
        process = subprocess.Popen(
            ['bcftools', 'view', '-H', file_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            errors='replace',
        )
    except (OSError, subprocess.SubprocessError):
        return None
    try:
        first = process.stdout.read(1)
    except (OSError, ValueError):
        first = ""
    if first:
        # There is at least one record: stop streaming the rest of the file
        try:
            process.kill()
        except OSError:
            pass
    try:
        process.stdout.close()
    except OSError:
        pass
    process.wait()
    if first:
        return True
    # No output: either the file really has no records or bcftools failed
    return False if process.returncode == 0 else None


def _terminate_process(process: subprocess.Popen) -> None:
    """Terminate a script started with start_new_session=True and its children"""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError, OSError):
        process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            process.kill()


# (path, size, mtime_ns) -> md5 hex digest; avoids re-reading a multi-GB VCF per analysis
_file_hash_memo: Dict[Tuple[str, int, int], str] = {}


def get_file_hash(file_path: str) -> str:
    """
    Calculate MD5 hash of a file for caching purposes

    The whole file is hashed (8 MB chunks) so a same-size edit anywhere in the
    file changes the hash. Digests are memoised per (path, size, mtime_ns).

    Args:
        file_path: Path to file

    Returns:
        MD5 hash string
    """
    try:
        st = os.stat(file_path)
        memo_key = (os.path.abspath(file_path), st.st_size, st.st_mtime_ns)
        cached = _file_hash_memo.get(memo_key)
        if cached is not None:
            return cached

        hash_md5 = hashlib.md5()
        chunk_size = 8 * 1024 * 1024  # 8MB
        with open(file_path, 'rb') as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                hash_md5.update(chunk)

        # Include file size in hash
        hash_md5.update(str(st.st_size).encode())

        digest = hash_md5.hexdigest()
        _file_hash_memo[memo_key] = digest
        return digest

    except (OSError, IOError):
        return ""


# Cache JSON format version. Entries without it (written by Helixight <= 1.2.0,
# whose file hash is identical for VCFs <= 1 MiB) carry results of the old
# parser and scripts and are ignored so the analysis re-runs and overwrites them.
CACHE_FORMAT = 2


class ResultCache:
    """
    Simple file-based cache for analysis results
    """

    def __init__(self, cache_dir: str = None):
        self.cache_dir = cache_dir or os.path.join(
            os.environ.get("RESULTS_DIR", "/app/results"),
            ".cache"
        )
        # The directory is created lazily by set(); importing this module or
        # constructing the cache must not touch the filesystem.

    def _get_cache_path(self, vcf_hash: str, analysis_id: str) -> str:
        """Get path for cached result"""
        return os.path.join(self.cache_dir, f"{vcf_hash}_{analysis_id}.json")

    def get(self, vcf_path: str, analysis_id: str, vcf_hash: str = None) -> Optional[AnalysisResult]:
        """
        Retrieve cached result if available

        Args:
            vcf_path: Path to VCF file
            analysis_id: Analysis identifier
            vcf_hash: Precomputed get_file_hash(vcf_path), if the caller has it

        Returns:
            Cached AnalysisResult or None
        """
        vcf_hash = vcf_hash or get_file_hash(vcf_path)
        if not vcf_hash:
            return None

        cache_path = self._get_cache_path(vcf_hash, analysis_id)

        if not os.path.exists(cache_path):
            return None

        try:
            with open(cache_path, 'r') as f:
                data = json.load(f)

            if data.get('format') != CACHE_FORMAT:
                return None

            return AnalysisResult(
                name=data['name'],
                status=AnalysisStatus(data['status']),
                output=data.get('output', ''),
                error=data.get('error', ''),
                scores=data.get('scores', {}),
                findings=data.get('findings', []),
                warnings=data.get('warnings', []),
                output_files=data.get('output_files', []),
                returncode=data.get('returncode'),
                stderr=data.get('stderr', ''),
                markers_found=data.get('markers_found'),
                markers_queried=data.get('markers_queried')
            )

        except (json.JSONDecodeError, KeyError, ValueError, TypeError, AttributeError, OSError):
            return None

    def set(self, vcf_path: str, analysis_id: str, result: AnalysisResult, vcf_hash: str = None) -> bool:
        """
        Cache an analysis result

        Args:
            vcf_path: Path to VCF file
            analysis_id: Analysis identifier
            result: AnalysisResult to cache
            vcf_hash: Precomputed get_file_hash(vcf_path), if the caller has it

        Returns:
            True if cached successfully
        """
        vcf_hash = vcf_hash or get_file_hash(vcf_path)
        if not vcf_hash:
            return False

        cache_path = self._get_cache_path(vcf_hash, analysis_id)

        try:
            data = {
                'format': CACHE_FORMAT,
                'name': result.name,
                'status': result.status.value,
                'output': result.output,
                'error': result.error,
                'scores': result.scores,
                'findings': result.findings,
                'warnings': result.warnings,
                'output_files': result.output_files,
                'returncode': result.returncode,
                'stderr': result.stderr,
                'markers_found': result.markers_found,
                'markers_queried': result.markers_queried,
                'cached_at': datetime.now().isoformat()
            }

            os.makedirs(self.cache_dir, exist_ok=True)
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
        if not os.path.isdir(self.cache_dir):
            return

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


# Global cache instance, created on first use (never at import time)
_result_cache: Optional[ResultCache] = None


def get_result_cache() -> ResultCache:
    """Return the process-wide ResultCache, creating it on first use"""
    global _result_cache
    if _result_cache is None:
        _result_cache = ResultCache()
    return _result_cache
