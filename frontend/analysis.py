"""
Helixight Analysis Module
Python wrapper for bash analysis scripts with progress tracking
"""

import subprocess
import os
import re
import json
import threading
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Callable, Dict, List, Any, Tuple
from enum import Enum


class AnalysisStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


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

    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path
        self.results: Dict[str, AnalysisResult] = {}
        self._progress_callback: Optional[Callable] = None
        self._cancel_flag = threading.Event()

    def set_progress_callback(self, callback: Callable[[str, float, str], None]):
        """Set callback for progress updates: (analysis_name, progress, message)"""
        self._progress_callback = callback

    def cancel(self):
        """Cancel running analyses"""
        self._cancel_flag.set()

    def run_analysis(self, analysis_id: str) -> AnalysisResult:
        """Run a single analysis script"""
        if analysis_id not in self.ANALYSES:
            return AnalysisResult(
                name=analysis_id,
                status=AnalysisStatus.FAILED,
                error=f"Unknown analysis: {analysis_id}"
            )

        analysis = self.ANALYSES[analysis_id]
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
