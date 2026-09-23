"""
Unit tests for the Helixight analysis module
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import json
import pytest

from analysis import (
    AnalysisStatus, AnalysisResult, VCFMetrics,
    AnalysisRunner, get_file_hash, ResultCache,
    get_vcf_metrics, get_result_cache, get_system_info, CACHE_FORMAT
)

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SAMPLE_VCF = os.path.join(REPO_ROOT, 'data', 'sample_genome.vcf')
SAMPLE_VCF_GZ = os.path.join(REPO_ROOT, 'data', 'sample_genome.vcf.gz')

needs_bcftools = pytest.mark.skipif(shutil.which('bcftools') is None, reason="bcftools not installed")

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=248956422>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE001\n"
)
VCF_TWO_RECORDS = VCF_HEADER + (
    "1\t100\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
    "1\t200\trs2\tC\tT\t50\tPASS\t.\tGT\t1/1\n"
)


def _write(path: str, content: str) -> str:
    with open(path, 'w') as f:
        f.write(content)
    return path


def _process_alive(pid: int) -> bool:
    """True when pid exists and is not a zombie"""
    try:
        with open(f"/proc/{pid}/stat") as f:
            state = f.read().rsplit(')', 1)[1].split()[0]
        return state != 'Z'
    except (OSError, IndexError):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True


class TestAnalysisStatus:
    """Tests for AnalysisStatus enum"""

    def test_status_values(self):
        assert AnalysisStatus.PENDING.value == "pending"
        assert AnalysisStatus.RUNNING.value == "running"
        assert AnalysisStatus.COMPLETED.value == "completed"
        assert AnalysisStatus.FAILED.value == "failed"


class TestAnalysisResult:
    """Tests for AnalysisResult dataclass"""

    def test_default_values(self):
        result = AnalysisResult(name="Test", status=AnalysisStatus.PENDING)
        assert result.name == "Test"
        assert result.status == AnalysisStatus.PENDING
        assert result.output == ""
        assert result.error == ""
        assert result.scores == {}
        assert result.findings == []
        assert result.warnings == []
        assert result.output_files == []
        assert result.returncode is None
        assert result.stderr == ""
        assert result.markers_found is None
        assert result.markers_queried is None
        assert result.from_cache is False

    def test_custom_values(self):
        result = AnalysisResult(
            name="Custom Test",
            status=AnalysisStatus.COMPLETED,
            output="test output",
            scores={"test": {"percentage": 75}},
            findings=["finding1"],
            warnings=["warning1"]
        )
        assert result.name == "Custom Test"
        assert result.status == AnalysisStatus.COMPLETED
        assert result.output == "test output"
        assert result.scores == {"test": {"percentage": 75}}
        assert result.findings == ["finding1"]
        assert result.warnings == ["warning1"]


class TestVCFMetrics:
    """Tests for VCFMetrics dataclass"""

    def test_default_values(self):
        metrics = VCFMetrics()
        assert metrics.total_variants == 0
        assert metrics.snps == 0
        assert metrics.indels == 0
        assert metrics.heterozygous == 0
        assert metrics.homozygous_alt == 0
        assert metrics.titv_ratio == 0.0
        assert metrics.file_size_mb == 0.0
        assert metrics.error == ""

    def test_custom_values(self):
        metrics = VCFMetrics(
            total_variants=1000,
            snps=800,
            indels=200,
            heterozygous=400,
            homozygous_alt=300,
            titv_ratio=2.1
        )
        assert metrics.total_variants == 1000
        assert metrics.snps == 800
        assert metrics.indels == 200
        assert metrics.titv_ratio == 2.1


class TestAnalysisRunner:
    """Tests for AnalysisRunner class"""

    def test_analyses_defined(self):
        """Verify all expected analyses are defined"""
        expected_analyses = [
            "athletic_genetics",
            "triathlon_genetics",
            "personalized_protocol",
            "mega_analysis",
            "fun_genetics",
            "clinvar_scan",
            "cancer_genes_scan"
        ]
        for analysis_id in expected_analyses:
            assert analysis_id in AnalysisRunner.ANALYSES
            assert "script" in AnalysisRunner.ANALYSES[analysis_id]
            assert "name" in AnalysisRunner.ANALYSES[analysis_id]
            assert "description" in AnalysisRunner.ANALYSES[analysis_id]

    def test_str_analysis_not_offered(self):
        """str_analysis.sh needs a BAM, which the web UI never has"""
        assert "str_analysis" not in AnalysisRunner.ANALYSES

    def test_every_script_exists(self):
        """Every analysis the UI offers must have its script on disk"""
        for analysis_id, info in AnalysisRunner.ANALYSES.items():
            script = os.path.join(REPO_ROOT, 'scripts', info['script'])
            assert os.path.isfile(script), f"{analysis_id}: missing {script}"

    def test_get_analyses_by_category(self):
        """Test analysis grouping by category"""
        categories = AnalysisRunner.get_analyses_by_category()
        assert isinstance(categories, dict)
        assert len(categories) > 0

        # Check that each category has analyses
        for cat, analyses in categories.items():
            assert isinstance(analyses, list)
            for analysis in analyses:
                assert "id" in analysis
                assert "name" in analysis

    def test_validate_vcf_file_not_found(self):
        """Test validation with non-existent file"""
        is_valid, message = AnalysisRunner.validate_vcf("/nonexistent/file.vcf")
        assert is_valid is False
        assert "not found" in message.lower()

    def test_validate_vcf_missing_header(self, tmp_path):
        path = _write(str(tmp_path / "nope.vcf"), "this is not a vcf\n")
        is_valid, message = AnalysisRunner.validate_vcf(path)
        assert is_valid is False
        assert "header" in message.lower()

    def test_unknown_analysis(self):
        """Test running unknown analysis"""
        runner = AnalysisRunner("/fake/path.vcf")
        result = runner.run_analysis("nonexistent_analysis")
        assert result.status == AnalysisStatus.FAILED
        assert "Unknown analysis" in result.error


@needs_bcftools
class TestValidateVCF:
    """validate_vcf / get_vcf_metrics on synthetic files (needs bcftools)"""

    def test_garbage_body_is_rejected(self, tmp_path):
        path = _write(str(tmp_path / "garbage.vcf"), VCF_HEADER + "this is not a vcf record\n")
        is_valid, message = AnalysisRunner.validate_vcf(path)
        assert is_valid is False
        assert "bcftools cannot parse" in message

    def test_header_only_is_valid_with_zero_records(self, tmp_path):
        path = _write(str(tmp_path / "empty.vcf"), VCF_HEADER)
        is_valid, message = AnalysisRunner.validate_vcf(path)
        assert is_valid is True
        assert "0 variant records" in message

    def test_two_record_file(self, tmp_path):
        path = _write(str(tmp_path / "two.vcf"), VCF_TWO_RECORDS)
        is_valid, message = AnalysisRunner.validate_vcf(path)
        assert is_valid is True
        assert message == "Valid VCF with 2 variants"

    def test_two_record_gzip_file(self, tmp_path):
        path = str(tmp_path / "two.vcf.gz")
        with gzip.open(path, 'wt') as f:
            f.write(VCF_TWO_RECORDS)
        is_valid, message = AnalysisRunner.validate_vcf(path)
        assert is_valid is True
        assert "2 variants" in message

    def test_metrics_count_genotypes(self, tmp_path):
        path = _write(str(tmp_path / "two.vcf"), VCF_TWO_RECORDS)
        metrics = get_vcf_metrics(path)
        assert metrics.total_variants == 2
        assert metrics.snps == 2
        assert metrics.heterozygous == 1
        assert metrics.homozygous_alt == 1
        assert metrics.homozygous_ref == 0
        assert metrics.quality_pass == 2
        assert metrics.error == ""

    def test_metrics_report_bcftools_error(self, tmp_path):
        path = _write(str(tmp_path / "garbage.vcf"), VCF_HEADER + "this is not a vcf record\n")
        metrics = get_vcf_metrics(path)
        assert metrics.total_variants == 0
        assert metrics.error != ""


class TestParseOutput:
    """Tests for AnalysisRunner._parse_output"""

    def _parse(self, text: str) -> AnalysisResult:
        result = AnalysisResult(name="t", status=AnalysisStatus.COMPLETED)
        AnalysisRunner("/fake/path.vcf")._parse_output(result, text)
        return result

    def test_pathogenic_genotype_lines_are_findings_not_scores(self):
        result = self._parse(
            "  🔴 BRCA2 4075delGT: 0/1 (Ashkenazi)\n"
            "  🔴 MLH1_c.350C>T: 0/1\n"
            "  POWER SCORE: 4 / 6\n"
        )
        assert result.findings == ["🔴 BRCA2 4075delGT: 0/1 (Ashkenazi)", "🔴 MLH1_c.350C>T: 0/1"]
        assert result.scores == {"POWER SCORE": {"value": 4, "max": 6, "percentage": 67}}

    def test_genotype_lines_without_keyword_are_not_scores(self):
        result = self._parse(
            "      Your genotype: 0/1\n"
            "  Date: 12 / 2024\n"
            "  SNPs / Indels: 5 / 3\n"
        )
        assert result.scores == {}

    def test_score_keyword_variants(self):
        result = self._parse(
            "  PGx SCORE: 85 / 100\n"
            "  CARRIER SCORE: 2 / 100\n"
            "Wynik CAD PRS: 5 / 20 punktow (~25 percentyl)\n"
            "WYNIK DLUGOWIECZNOSCI: 5 / 9 punktow\n"
            "Wynik dlugowiecznosci: 3 / 7\n"
        )
        assert result.scores["PGx SCORE"] == {"value": 85, "max": 100, "percentage": 85}
        assert result.scores["CARRIER SCORE"]["percentage"] == 2
        assert result.scores["Wynik CAD PRS"] == {"value": 5, "max": 20, "percentage": 25}
        assert result.scores["WYNIK DLUGOWIECZNOSCI"]["max"] == 9
        assert result.scores["Wynik dlugowiecznosci"]["max"] == 7

    def test_triathlon_discipline_scores(self):
        result = self._parse("   🏊 PŁYWANIE:  3 / 15\n   🚴 KOLARSTWO: 3 / 20\n   🏃 BIEGANIE:  4 / 20\n")
        assert result.scores["PŁYWANIE"] == {"value": 3, "max": 15, "percentage": 20}
        assert result.scores["KOLARSTWO"] == {"value": 3, "max": 20, "percentage": 15}
        assert result.scores["BIEGANIE"] == {"value": 4, "max": 20, "percentage": 20}

    def test_percentage_lines(self):
        result = self._parse("   Procent potencjalnie rzadkich/nowych: 12%\n  Bogus: 250%\n")
        assert result.scores["Procent potencjalnie rzadkich/nowych"] == {"percentage": 12}
        assert result.scores["Bogus"]["percentage"] == 100

    def test_first_score_wins(self):
        result = self._parse("POWER SCORE: 4 / 6\nPOWER SCORE: 1 / 6\n")
        assert result.scores["POWER SCORE"]["value"] == 4

    def test_warning_lands_only_in_warnings(self):
        result = self._parse("⚠️  FOUND 1 pathogenic variants!\n⚠️ Low coverage in region X\nUWAGA: cos\n")
        assert "⚠️  FOUND 1 pathogenic variants!" in result.findings   # real hit (FOUND)
        assert result.warnings == ["⚠️ Low coverage in region X", "UWAGA: cos"]
        assert "⚠️ Low coverage in region X" not in result.findings

    def test_disclaimers_are_dropped(self):
        result = self._parse(
            "⚠️  DISCLAIMER: This is a genetic PREDICTION only.\n"
            "⚠️  NOT FOR CLINICAL USE - EDUCATIONAL PURPOSES ONLY\n"
        )
        assert result.findings == []
        assert result.warnings == []

    def test_clinvar_progress_lines_are_ignored(self):
        result = self._parse(
            "🔬 Extracting Pathogenic/Likely_pathogenic variants from ClinVar...\n"
            "   ClinVar contains 5 Pathogenic/Likely_pathogenic variants\n"
        )
        assert result.findings == []
        assert result.warnings == []

    def test_clinvar_status_lines_are_not_findings(self):
        """clinvar_scan.sh housekeeping lines carry no ✅ marker any more"""
        result = self._parse("ℹ️  Using existing ClinVar database\n"
                            "   md5 checksum verified\n"
                            "ClinVar downloaded successfully!\n"
                            "✅ No pathogenic variants found in your genome!\n")
        assert result.findings == ["✅ No pathogenic variants found in your genome!"]
        assert result.warnings == []

    def test_positive_finding_once(self):
        result = self._parse("✅ No pathogenic variants found\n")
        assert result.findings == ["✅ No pathogenic variants found"]
        assert result.warnings == []

    def test_coverage_line(self):
        result = self._parse("some text\nHX_DATA_COVERAGE: 7/20 markers\n")
        assert (result.markers_found, result.markers_queried) == (7, 20)
        assert result.scores == {}

    def test_warning_marker_without_presentation_selector(self):
        result = self._parse("\u26a0 plain warning sign\n")
        assert result.warnings == ["\u26a0 plain warning sign"]
        assert result.findings == []


class TestFileHash:
    """Tests for file hashing functionality"""

    def test_hash_nonexistent_file(self):
        """Test hashing non-existent file returns empty string"""
        result = get_file_hash("/nonexistent/file")
        assert result == ""

    def test_hash_consistency(self):
        """Test that same file produces same hash"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("test content for hashing")
            temp_path = f.name

        try:
            hash1 = get_file_hash(temp_path)
            hash2 = get_file_hash(temp_path)
            assert hash1 == hash2
            assert len(hash1) == 32  # MD5 hash length
        finally:
            os.unlink(temp_path)

    def test_different_files_different_hash(self):
        """Test that different files produce different hashes"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f1:
            f1.write("content one")
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f2:
            f2.write("content two")
            path2 = f2.name

        try:
            hash1 = get_file_hash(path1)
            hash2 = get_file_hash(path2)
            assert hash1 != hash2
        finally:
            os.unlink(path1)
            os.unlink(path2)

    @staticmethod
    def _pair(tmp_path, size: int, diff_at: int):
        data = bytearray(b'A' * size)
        a = tmp_path / "a.vcf"
        b = tmp_path / "b.vcf"
        a.write_bytes(bytes(data))
        data[diff_at] = ord('B')
        b.write_bytes(bytes(data))
        return str(a), str(b)

    def test_middle_edit_changes_hash(self, tmp_path):
        """3 MB files differing only at byte 1.5 MB hash differently"""
        a, b = self._pair(tmp_path, 3 * 1024 * 1024, 1536 * 1024)
        assert get_file_hash(a) != get_file_hash(b)

    def test_last_byte_edit_changes_hash(self, tmp_path):
        """2 MB files differing only in the last byte hash differently"""
        size = 2 * 1024 * 1024
        a, b = self._pair(tmp_path, size, size - 1)
        assert get_file_hash(a) != get_file_hash(b)


class TestResultCache:
    """Tests for ResultCache class"""

    def test_cache_initialization_is_lazy(self):
        """Constructing a cache does not touch the filesystem; set() creates the dir"""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = os.path.join(tmpdir, "test_cache")
            cache = ResultCache(cache_dir)
            assert cache.cache_dir == cache_dir
            assert not os.path.exists(cache_dir)
            cache.clear()  # no dir yet: must not raise

            test_file = os.path.join(tmpdir, "test.vcf")
            with open(test_file, 'w') as f:
                f.write("##fileformat=VCFv4.2\ntest content")
            assert cache.set(test_file, "a", AnalysisResult(name="A", status=AnalysisStatus.COMPLETED))
            assert os.path.isdir(cache_dir)

    def test_import_has_no_side_effects(self):
        """`import analysis` must succeed with an unwritable RESULTS_DIR"""
        env = {**os.environ, 'RESULTS_DIR': '/proc/nope'}
        proc = subprocess.run(
            [sys.executable, '-c', 'import analysis; assert analysis._result_cache is None'],
            cwd=os.path.join(REPO_ROOT, 'frontend'), env=env,
            capture_output=True, text=True, timeout=60
        )
        assert proc.returncode == 0, proc.stderr

    def test_get_result_cache_singleton(self):
        assert get_result_cache() is get_result_cache()

    def test_cache_miss(self):
        """Test cache miss returns None"""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = ResultCache(os.path.join(tmpdir, "cache"))
            result = cache.get("/fake/path.vcf", "test_analysis")
            assert result is None

    def test_cache_set_and_get(self):
        """Test caching and retrieving results"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a test file to hash
            test_file = os.path.join(tmpdir, "test.vcf")
            with open(test_file, 'w') as f:
                f.write("##fileformat=VCFv4.2\ntest content")

            cache = ResultCache(os.path.join(tmpdir, "cache"))

            # Create a result to cache
            result = AnalysisResult(
                name="Test Analysis",
                status=AnalysisStatus.COMPLETED,
                output="test output",
                scores={"score1": {"percentage": 80}},
                findings=["finding1"],
                returncode=0,
                stderr="some diagnostics",
                markers_found=3,
                markers_queried=5,
            )

            # Cache the result
            success = cache.set(test_file, "test_analysis", result)
            assert success is True

            # Retrieve from cache
            cached = cache.get(test_file, "test_analysis")
            assert cached is not None
            assert cached.name == result.name
            assert cached.status == result.status
            assert cached.output == result.output
            assert cached.scores == result.scores
            assert cached.findings == result.findings
            assert cached.returncode == 0
            assert cached.stderr == "some diagnostics"
            assert (cached.markers_found, cached.markers_queried) == (3, 5)
            assert cached.from_cache is False  # runtime flag, set by the runner

            # cached_at is a real ISO timestamp
            cache_files = os.listdir(cache.cache_dir)
            assert len(cache_files) == 1
            with open(os.path.join(cache.cache_dir, cache_files[0])) as f:
                data = json.load(f)
            assert isinstance(data['cached_at'], str)
            assert data['cached_at'][:2] == "20"
            assert 'T' in data['cached_at']

    def test_pre_1_3_cache_entries_are_ignored(self):
        """Cache JSON without a format version (Helixight <= 1.2.0) is a miss:
        for VCFs <= 1 MiB the old hash equals the new one, so those entries
        would otherwise be served with the old parser's scores and findings"""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = os.path.join(tmpdir, "test.vcf")
            with open(test_file, 'w') as f:
                f.write("##fileformat=VCFv4.2\ntest content")
            cache = ResultCache(os.path.join(tmpdir, "cache"))
            os.makedirs(cache.cache_dir)
            with open(cache._get_cache_path(get_file_hash(test_file), "old"), 'w') as f:
                json.dump({'name': 'Old', 'status': 'completed', 'output': 'x'}, f)
            assert cache.get(test_file, "old") is None
            # A wrong version is a miss too
            with open(cache._get_cache_path(get_file_hash(test_file), "v1"), 'w') as f:
                json.dump({'format': 1, 'name': 'Old', 'status': 'completed', 'output': 'x'}, f)
            assert cache.get(test_file, "v1") is None

    def test_current_format_entries_tolerate_missing_optional_fields(self):
        """Within the current format, optional fields may be absent"""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = os.path.join(tmpdir, "test.vcf")
            with open(test_file, 'w') as f:
                f.write("##fileformat=VCFv4.2\ntest content")
            cache = ResultCache(os.path.join(tmpdir, "cache"))
            os.makedirs(cache.cache_dir)
            with open(cache._get_cache_path(get_file_hash(test_file), "cur"), 'w') as f:
                json.dump({'format': CACHE_FORMAT, 'name': 'Old', 'status': 'completed', 'output': 'x'}, f)
            cached = cache.get(test_file, "cur")
            assert cached is not None
            assert cached.returncode is None
            assert cached.markers_found is None

    def test_cache_set_writes_format_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = os.path.join(tmpdir, "test.vcf")
            with open(test_file, 'w') as f:
                f.write("##fileformat=VCFv4.2\ntest content")
            cache = ResultCache(os.path.join(tmpdir, "cache"))
            cache.set(test_file, "a", AnalysisResult(name="A", status=AnalysisStatus.COMPLETED))
            with open(cache._get_cache_path(get_file_hash(test_file), "a")) as f:
                assert json.load(f)['format'] == CACHE_FORMAT == 2

    def test_cache_miss_after_same_length_edit(self):
        """An in-place edit that keeps the file size must not hit the cache"""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = os.path.join(tmpdir, "test.vcf")
            with open(test_file, 'wb') as f:
                f.write(b"##fileformat=VCFv4.2\n" + b"A" * (1536 * 1024))
            cache = ResultCache(os.path.join(tmpdir, "cache"))
            cache.set(test_file, "a", AnalysisResult(name="A", status=AnalysisStatus.COMPLETED))
            assert cache.get(test_file, "a") is not None

            with open(test_file, 'r+b') as f:
                f.seek(1024 * 1024 + 123)
                f.write(b"B")
            # A later upload has a later mtime; make that explicit so the
            # (path, size, mtime) memo cannot mask the content change.
            st = os.stat(test_file)
            os.utime(test_file, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))

            assert cache.get(test_file, "a") is None

    def test_cache_clear(self):
        """Test clearing cache"""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = os.path.join(tmpdir, "test.vcf")
            with open(test_file, 'w') as f:
                f.write("##fileformat=VCFv4.2\ntest content")

            cache = ResultCache(os.path.join(tmpdir, "cache"))

            result = AnalysisResult(
                name="Test",
                status=AnalysisStatus.COMPLETED
            )
            cache.set(test_file, "test_analysis", result)

            # Verify it's cached
            assert cache.get(test_file, "test_analysis") is not None

            # Clear cache for specific file
            cache.clear(test_file)

            # Verify it's cleared
            assert cache.get(test_file, "test_analysis") is None


class TestRunAnalysis:
    """run_analysis against stub scripts in a temporary SCRIPTS_DIR"""

    @pytest.fixture
    def runner_env(self, tmp_path, monkeypatch):
        scripts_dir = tmp_path / "scripts"
        results_dir = tmp_path / "results"
        scripts_dir.mkdir()
        vcf = tmp_path / "genome.vcf"
        vcf.write_text(VCF_TWO_RECORDS)

        analyses = {}

        def add_script(analysis_id: str, body: str, **meta):
            path = scripts_dir / f"{analysis_id}.sh"
            path.write_text("#!/bin/bash\n" + body)
            analyses[analysis_id] = {"script": f"{analysis_id}.sh", "name": analysis_id.title(),
                                     "description": "", "icon": "", "category": "test", **meta}

        monkeypatch.setattr(AnalysisRunner, "ANALYSES", analyses)
        monkeypatch.setattr(AnalysisRunner, "SCRIPTS_DIR", str(scripts_dir))
        monkeypatch.setattr(AnalysisRunner, "RESULTS_DIR", str(results_dir))
        import analysis as analysis_module
        monkeypatch.setattr(analysis_module, "_result_cache", ResultCache(str(results_dir / ".cache")))
        return {"add_script": add_script, "vcf": str(vcf), "results_dir": str(results_dir)}

    def test_completed_then_cached(self, runner_env):
        runner_env["add_script"]("ok", 'echo "  POWER SCORE: 4 / 6"\necho "✅ fine"\n'
                                        'echo "$HELIXIGHT_WORK_DIR" > env.txt\n'
                                        'echo "HX_DATA_COVERAGE: 3/5 markers"\n')
        runner = AnalysisRunner(runner_env["vcf"])
        result = runner.run_analysis("ok")
        assert result.status == AnalysisStatus.COMPLETED
        assert result.returncode == 0
        assert result.from_cache is False
        assert result.scores == {"POWER SCORE": {"value": 4, "max": 6, "percentage": 67}}
        assert result.findings == ["✅ fine"]
        assert (result.markers_found, result.markers_queried) == (3, 5)
        # Per-VCF, per-analysis working directory; files written there are collected
        assert len(result.output_files) == 1
        env_file = result.output_files[0]
        assert env_file.startswith(os.path.join(runner_env["results_dir"], "genome_"))
        assert os.sep + "ok" + os.sep in env_file
        with open(env_file) as f:
            assert f.read().strip() == os.path.join(runner_env["results_dir"], ".work")

        second = AnalysisRunner(runner_env["vcf"]).run_analysis("ok")
        assert second.status == AnalysisStatus.COMPLETED
        assert second.from_cache is True
        assert second.scores == result.scores

    def test_output_files_lists_only_this_runs_files(self, runner_env):
        """run_dir is reused across runs of the same VCF/analysis; earlier
        runs' files must not be reported as outputs of the current run"""
        runner_env["add_script"]("writer", 'echo "run" > "report_$RANDOM$RANDOM.txt"\n')
        first = AnalysisRunner(runner_env["vcf"], use_cache=False).run_analysis("writer")
        assert len(first.output_files) == 1
        second = AnalysisRunner(runner_env["vcf"], use_cache=False).run_analysis("writer")
        assert len(second.output_files) == 1
        assert second.output_files != first.output_files
        assert os.path.exists(first.output_files[0])
        third = AnalysisRunner(runner_env["vcf"]).run_analysis("writer", force_rerun=True)
        assert len(third.output_files) == 1

    def test_reference_dir_in_environment(self, runner_env):
        runner_env["add_script"]("envcheck", 'echo "REF=$HELIXIGHT_REFERENCE_DIR"\n')
        result = AnalysisRunner(runner_env["vcf"], use_cache=False).run_analysis("envcheck")
        assert result.status == AnalysisStatus.COMPLETED
        assert f"REF={os.path.join(runner_env['results_dir'], '.reference')}" in result.output

    def test_exit_1_without_stderr(self, runner_env):
        runner_env["add_script"]("bad", 'echo "partial output"\nexit 1\n')
        result = AnalysisRunner(runner_env["vcf"]).run_analysis("bad")
        assert result.status == AnalysisStatus.FAILED
        assert result.returncode == 1
        assert "exited with code 1" in result.error
        assert "partial output" in result.output
        # Failures are never cached
        assert AnalysisRunner(runner_env["vcf"]).run_analysis("bad").from_cache is False

    def test_exit_1_with_stderr(self, runner_env):
        runner_env["add_script"]("bad2", 'echo "ERROR: boom" >&2\nexit 2\n')
        result = AnalysisRunner(runner_env["vcf"]).run_analysis("bad2")
        assert result.status == AnalysisStatus.FAILED
        assert result.error == "ERROR: boom"
        assert result.stderr.strip() == "ERROR: boom"

    def test_zero_coverage_is_failure_and_not_cached(self, runner_env):
        runner_env["add_script"]("nodata", 'echo "✅ No variants found"\necho "HX_DATA_COVERAGE: 0/12 markers"\n')
        result = AnalysisRunner(runner_env["vcf"]).run_analysis("nodata")
        assert result.status == AnalysisStatus.FAILED
        assert result.returncode == 0
        assert "No usable variant data" in result.error
        assert "0/12 markers" in result.error
        assert (result.markers_found, result.markers_queried) == (0, 12)
        again = AnalysisRunner(runner_env["vcf"]).run_analysis("nodata")
        assert again.from_cache is False
        assert again.status == AnalysisStatus.FAILED

    def test_zero_coverage_is_a_result_when_absence_is_result(self, runner_env):
        """Rare-allele screens (carrier_screening, blood_type): 0/N is the negative result"""
        runner_env["add_script"]("carrier", 'echo "✅ NO CARRIER VARIANTS DETECTED"\n'
                                            'echo "HX_DATA_COVERAGE: 0/29 markers"\n',
                                 absence_is_result=True)
        result = AnalysisRunner(runner_env["vcf"]).run_analysis("carrier")
        assert result.status == AnalysisStatus.COMPLETED
        assert result.error == ""
        assert result.findings == ["✅ NO CARRIER VARIANTS DETECTED"]
        assert (result.markers_found, result.markers_queried) == (0, 29)
        assert result.warnings[0].startswith("⚠️ None of the 29 markers")
        # Not cached: a later VCF with more data must re-run
        again = AnalysisRunner(runner_env["vcf"]).run_analysis("carrier")
        assert again.from_cache is False
        assert again.status == AnalysisStatus.COMPLETED

    @needs_bcftools
    def test_zero_coverage_on_empty_vcf_is_failure_even_when_absence_is_result(self, runner_env, tmp_path):
        """A header-only VCF (0 records) is no evidence of anything: the
        absence_is_result carve-out must not turn it into a negative result"""
        empty = tmp_path / "empty.vcf"
        empty.write_text(VCF_HEADER)
        runner_env["add_script"]("carrier", 'echo "✅ NO CARRIER VARIANTS DETECTED"\n'
                                            'echo "HX_DATA_COVERAGE: 0/29 markers"\n',
                                 absence_is_result=True)
        result = AnalysisRunner(str(empty)).run_analysis("carrier")
        assert result.status == AnalysisStatus.FAILED
        assert "No usable variant data" in result.error
        assert "no variant records" in result.error
        assert AnalysisRunner(str(empty)).run_analysis("carrier").from_cache is False

    def test_absence_is_result_flag_on_real_analyses(self):
        for analysis_id in ("carrier_screening", "blood_type"):
            assert AnalysisRunner.ANALYSES[analysis_id].get("absence_is_result") is True
        assert not AnalysisRunner.ANALYSES["athletic_genetics"].get("absence_is_result")

    def test_ansi_codes_stripped_from_stdout(self, runner_env):
        runner_env["add_script"]("colour", "echo -e '\\033[0;32m  ✅ ok\\033[0m'\n"
                                           "echo -e '\\033[0;31mwarning: x\\033[0m' >&2\n")
        result = AnalysisRunner(runner_env["vcf"], use_cache=False).run_analysis("colour")
        assert result.status == AnalysisStatus.COMPLETED
        assert result.findings == ["✅ ok"]
        assert "\x1b" not in result.output
        assert "\x1b" not in result.stderr

    def test_zero_of_zero_coverage_is_not_a_failure(self, runner_env):
        """A script that queried no markers at all is not 'no data'"""
        runner_env["add_script"]("noquery", 'echo "hello"\necho "HX_DATA_COVERAGE: 0/0 markers"\n')
        result = AnalysisRunner(runner_env["vcf"], use_cache=False).run_analysis("noquery")
        assert result.status == AnalysisStatus.COMPLETED
        assert (result.markers_found, result.markers_queried) == (0, 0)

    def test_cancel_kills_process_group(self, runner_env, tmp_path):
        pid_file = tmp_path / "sleep.pid"
        runner_env["add_script"]("slow", f'sleep 4 &\necho $! > "{pid_file}"\nwait\n')
        runner = AnalysisRunner(runner_env["vcf"], use_cache=False)
        holder = {}

        def run():
            holder["result"] = runner.run_analysis("slow")

        t = threading.Thread(target=run)
        start = time.monotonic()
        t.start()
        time.sleep(1)
        runner.cancel()
        t.join(timeout=10)
        elapsed = time.monotonic() - start

        assert not t.is_alive()
        assert elapsed < 3.5
        result = holder["result"]
        assert result.status == AnalysisStatus.FAILED
        assert result.error == "Analysis cancelled"
        sleep_pid = int(pid_file.read_text().strip())
        time.sleep(0.2)
        assert not _process_alive(sleep_pid)

    def test_timeout(self, runner_env, monkeypatch):
        import analysis as analysis_module
        monkeypatch.setattr(analysis_module, "ANALYSIS_TIMEOUT", 1)
        runner_env["add_script"]("hang", 'sleep 5\n')
        start = time.monotonic()
        result = AnalysisRunner(runner_env["vcf"], use_cache=False).run_analysis("hang")
        assert time.monotonic() - start < 4
        assert result.status == AnalysisStatus.FAILED
        assert "Timed out after 1s" in result.error


    def test_timeout_kills_grandchild_that_ignores_sigterm(self, runner_env, monkeypatch, tmp_path):
        """bash dies on SIGTERM at once, but a child that ignores it still holds
        the stdout pipe: the deadline must not wait for it to finish"""
        import analysis as analysis_module
        monkeypatch.setattr(analysis_module, "ANALYSIS_TIMEOUT", 1)
        pid_file = tmp_path / "orphan.pid"
        runner_env["add_script"]("orphan", "bash -c 'trap \"\" TERM; sleep 15' &\n"
                                           f'echo $! > "{pid_file}"\nwait\n')
        start = time.monotonic()
        result = AnalysisRunner(runner_env["vcf"], use_cache=False).run_analysis("orphan")
        assert time.monotonic() - start < 9
        assert result.status == AnalysisStatus.FAILED
        assert "Timed out after 1s" in result.error
        orphan_pid = int(pid_file.read_text().strip())
        time.sleep(0.2)
        assert not _process_alive(orphan_pid)


class TestSystemInfo:

    def test_non_utf8_version_banner_does_not_raise(self, monkeypatch):
        """samtools --version can contain non-UTF-8 bytes; the sidebar must not die"""
        import analysis as analysis_module
        real_run = subprocess.run

        def fake_run(cmd, *args, **kwargs):
            if cmd[0] == 'samtools':
                assert kwargs.get('errors') == 'replace'
                return subprocess.CompletedProcess(cmd, 0, stdout="samtools 1.19.2\nBUILD \ufffdX\ufffd\n", stderr="")
            return real_run(cmd, *args, **kwargs)

        monkeypatch.setattr(analysis_module.subprocess, "run", fake_run)
        info = get_system_info()
        assert info['samtools'] == "samtools 1.19.2"


class TestSampleVCF:
    """Tests for sample VCF file"""

    @pytest.fixture
    def sample_vcf_path(self):
        """Path to sample VCF file"""
        return SAMPLE_VCF

    def test_sample_vcf_exists(self, sample_vcf_path):
        """Test that sample VCF file exists"""
        assert os.path.exists(sample_vcf_path), "Sample VCF file should exist"
        assert os.path.exists(SAMPLE_VCF_GZ), "Compressed sample VCF file should exist"

    def test_sample_vcf_valid_header(self, sample_vcf_path):
        """Test that sample VCF has valid header"""
        with open(sample_vcf_path, 'r') as f:
            first_line = f.readline()
        assert first_line.startswith('##fileformat=VCF')

    def test_sample_vcf_has_variants(self, sample_vcf_path):
        """Test that sample VCF contains variants"""
        variant_count = 0
        with open(sample_vcf_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    variant_count += 1
        assert variant_count > 0, "Sample VCF should contain variants"

    @staticmethod
    def _records(path):
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as f:
            return [line for line in f if line.strip() and not line.startswith('#')]

    def test_sample_vcf_no_duplicate_positions(self):
        seen = set()
        for line in self._records(SAMPLE_VCF_GZ):
            chrom, pos = line.split('\t')[:2]
            assert (chrom, pos) not in seen, f"duplicate record at {chrom}:{pos}"
            seen.add((chrom, pos))

    @needs_bcftools
    def test_bcftools_parses_sample(self):
        """bcftools must read every record of the shipped sample"""
        proc = subprocess.run(['bcftools', 'view', '-H', SAMPLE_VCF_GZ],
                              capture_output=True, text=True, timeout=120)
        assert proc.returncode == 0, proc.stderr
        parsed = [line for line in proc.stdout.splitlines() if line.strip()]
        assert len(parsed) == len(self._records(SAMPLE_VCF_GZ))
        assert len(parsed) > 0

    @needs_bcftools
    def test_validate_and_metrics_on_sample(self):
        is_valid, message = AnalysisRunner.validate_vcf(SAMPLE_VCF_GZ)
        assert is_valid is True, message
        assert "0 variant records" not in message
        metrics = get_vcf_metrics(SAMPLE_VCF_GZ)
        assert metrics.total_variants == len(self._records(SAMPLE_VCF_GZ))
        assert metrics.heterozygous + metrics.homozygous_alt + metrics.homozygous_ref > 0
        assert metrics.error == ""


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
