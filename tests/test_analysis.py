"""
Unit tests for the Helixight analysis module
"""

import os
import sys
import tempfile
import json
import pytest

# Add frontend to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'frontend'))

from analysis import (
    AnalysisStatus, AnalysisResult, VCFMetrics,
    AnalysisRunner, get_file_hash, ResultCache,
    get_vcf_metrics
)


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

    def test_unknown_analysis(self):
        """Test running unknown analysis"""
        runner = AnalysisRunner("/fake/path.vcf")
        result = runner.run_analysis("nonexistent_analysis")
        assert result.status == AnalysisStatus.FAILED
        assert "Unknown analysis" in result.error


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


class TestResultCache:
    """Tests for ResultCache class"""

    def test_cache_initialization(self):
        """Test cache directory creation"""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = os.path.join(tmpdir, "test_cache")
            cache = ResultCache(cache_dir)
            assert os.path.exists(cache_dir)

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
                findings=["finding1"]
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


class TestSampleVCF:
    """Tests for sample VCF file"""

    @pytest.fixture
    def sample_vcf_path(self):
        """Path to sample VCF file"""
        return os.path.join(
            os.path.dirname(__file__), '..', 'data', 'sample_genome.vcf'
        )

    def test_sample_vcf_exists(self, sample_vcf_path):
        """Test that sample VCF file exists"""
        assert os.path.exists(sample_vcf_path), "Sample VCF file should exist"

    def test_sample_vcf_valid_header(self, sample_vcf_path):
        """Test that sample VCF has valid header"""
        if os.path.exists(sample_vcf_path):
            with open(sample_vcf_path, 'r') as f:
                first_line = f.readline()
            assert first_line.startswith('##fileformat=VCF')

    def test_sample_vcf_has_variants(self, sample_vcf_path):
        """Test that sample VCF contains variants"""
        if os.path.exists(sample_vcf_path):
            variant_count = 0
            with open(sample_vcf_path, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        variant_count += 1
            assert variant_count > 0, "Sample VCF should contain variants"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
