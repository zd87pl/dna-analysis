"""
Unit tests for the PDF report generation module
"""

import os
import sys
import pytest

# Add frontend to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'frontend'))

from analysis import AnalysisStatus, AnalysisResult
from pdf_report import generate_pdf_report, HelixightReport


class TestHelixightReport:
    """Tests for HelixightReport PDF class"""

    def test_report_initialization(self):
        """Test PDF report can be initialized"""
        report = HelixightReport()
        assert report is not None

    def test_add_page(self):
        """Test adding a page"""
        report = HelixightReport()
        report.add_page()
        assert report.page_no() == 1

    def test_chapter_title(self):
        """Test adding chapter title"""
        report = HelixightReport()
        report.add_page()
        # Should not raise
        report.chapter_title("Test Chapter", "🧬")

    def test_add_score_bar(self):
        """Test adding score bar"""
        report = HelixightReport()
        report.add_page()
        # Should not raise
        report.add_score_bar("Test Score", 75)

    def test_add_finding(self):
        """Test adding findings"""
        report = HelixightReport()
        report.add_page()
        # Should not raise
        report.add_finding("Test finding", "info")
        report.add_finding("Positive finding", "positive")
        report.add_finding("Warning finding", "warning")


class TestGeneratePDFReport:
    """Tests for generate_pdf_report function"""

    @pytest.fixture
    def sample_results(self):
        """Create sample analysis results"""
        return {
            "athletic_genetics": AnalysisResult(
                name="Athletic Genetics",
                status=AnalysisStatus.COMPLETED,
                output="Test output",
                scores={"Power": {"percentage": 75}, "Endurance": {"percentage": 60}},
                findings=["Good power potential", "Average endurance"],
                warnings=[]
            ),
            "fun_genetics": AnalysisResult(
                name="Fun Genetics",
                status=AnalysisStatus.COMPLETED,
                output="Fun output",
                scores={"Caffeine": {"percentage": 90}},
                findings=["Fast caffeine metabolism"],
                warnings=["Consider limiting intake"]
            )
        }

    @pytest.fixture
    def sample_metadata(self):
        """Create sample analysis metadata"""
        return {
            "athletic_genetics": {
                "name": "Athletic Genetics",
                "icon": "🏃",
                "description": "Athletic performance analysis"
            },
            "fun_genetics": {
                "name": "Fun Genetics",
                "icon": "🎲",
                "description": "Fun trait analysis"
            }
        }

    def test_generate_pdf_basic(self, sample_results, sample_metadata):
        """Test basic PDF generation"""
        pdf_bytes = generate_pdf_report(
            vcf_filename="test.vcf",
            results=sample_results,
            analyses_metadata=sample_metadata
        )

        assert pdf_bytes is not None
        assert len(pdf_bytes) > 0
        # PDF files start with %PDF
        assert pdf_bytes[:4] == b'%PDF'

    def test_generate_pdf_with_metrics(self, sample_results, sample_metadata):
        """Test PDF generation with VCF metrics"""
        vcf_metrics = {
            'total_variants': 1000,
            'snps': 800,
            'indels': 200,
            'heterozygous': 400,
            'homozygous_alt': 300,
            'titv_ratio': 2.1,
            'file_size_mb': 5.5
        }

        pdf_bytes = generate_pdf_report(
            vcf_filename="test.vcf",
            results=sample_results,
            analyses_metadata=sample_metadata,
            vcf_metrics=vcf_metrics
        )

        assert pdf_bytes is not None
        assert len(pdf_bytes) > 0
        assert pdf_bytes[:4] == b'%PDF'

    def test_generate_pdf_empty_results(self, sample_metadata):
        """Test PDF generation with empty results"""
        pdf_bytes = generate_pdf_report(
            vcf_filename="empty.vcf",
            results={},
            analyses_metadata=sample_metadata
        )

        assert pdf_bytes is not None
        assert pdf_bytes[:4] == b'%PDF'

    def test_generate_pdf_failed_analysis(self, sample_metadata):
        """Test PDF generation with failed analysis"""
        results = {
            "athletic_genetics": AnalysisResult(
                name="Athletic Genetics",
                status=AnalysisStatus.FAILED,
                error="Script not found"
            )
        }

        pdf_bytes = generate_pdf_report(
            vcf_filename="test.vcf",
            results=results,
            analyses_metadata=sample_metadata
        )

        assert pdf_bytes is not None
        assert pdf_bytes[:4] == b'%PDF'

    def test_generate_pdf_special_characters(self, sample_metadata):
        """Test PDF generation with special characters in findings"""
        results = {
            "fun_genetics": AnalysisResult(
                name="Fun Genetics",
                status=AnalysisStatus.COMPLETED,
                findings=[
                    "✅ No pathogenic variants found",
                    "🔴 FOUND: Risk variant",
                    "⚠️ Warning message"
                ],
                scores={}
            )
        }

        pdf_bytes = generate_pdf_report(
            vcf_filename="test.vcf",
            results=results,
            analyses_metadata=sample_metadata
        )

        assert pdf_bytes is not None
        assert pdf_bytes[:4] == b'%PDF'


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
