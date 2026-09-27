"""
Unit tests for the PDF report generation module
"""

import pytest

from analysis import AnalysisStatus, AnalysisResult
from pdf_report import generate_pdf_report, HelixightReport, _to_latin1


class TestHelixightReport:
    """Tests for HelixightReport PDF class"""

    def test_report_initialization(self):
        """Test PDF report can be initialized"""
        report = HelixightReport()
        assert report.page_no() == 0

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
        assert report.page_no() == 1

    def test_add_score_bar(self):
        """Test adding score bar"""
        report = HelixightReport()
        report.add_page()
        # Should not raise
        report.add_score_bar("Test Score", 75)
        assert report.page_no() == 1

    def test_add_score_bar_non_latin1_name(self):
        """Polish score names from triathlon_genetics.sh must not crash the core font"""
        report = HelixightReport()
        report.add_page()
        report.add_score_bar("PŁYWANIE", 60)
        assert report.page_no() == 1

    def test_add_score_bar_out_of_range(self):
        """Percentages outside 0-100 are clamped instead of drawing off the page"""
        report = HelixightReport()
        report.add_page()
        # fpdf2 happily draws a 250-wide or negative-width rect, so record the
        # widths actually passed to rect(): background, fill, background, fill.
        widths = []
        orig_rect = report.rect

        def spy_rect(x, y, w, h, *args, **kwargs):
            widths.append(w)
            return orig_rect(x, y, w, h, *args, **kwargs)

        report.rect = spy_rect
        report.add_score_bar("Too big", 250)
        report.add_score_bar("Negative", -5)
        assert report.page_no() == 1
        assert widths == [100, 100, 100, 0]

    def test_add_finding(self):
        """Test adding findings"""
        report = HelixightReport()
        report.add_page()
        # Should not raise
        report.add_finding("Test finding", "info")
        report.add_finding("Positive finding", "positive")
        report.add_finding("Warning finding", "warning")
        assert report.page_no() == 1

    def test_consecutive_findings_do_not_crash(self):
        """Two consecutive multi_cell(0, ...) calls used to leave zero width for the second"""
        report = HelixightReport()
        report.add_page()
        for i in range(5):
            report.add_finding(f"⚠️ finding number {i} with some longer text to wrap " * 3, "warning")
        assert report.page_no() >= 1


class TestToLatin1:
    """Tests for the core-font sanitiser"""

    def test_markers_become_ascii_tags(self):
        assert _to_latin1("✅ ok") == "[OK] ok"
        assert _to_latin1("🔴 hit") == "[!] hit"
        assert _to_latin1("⚠️ warn") == "[!] warn"
        assert _to_latin1("❌ no") == "[X] no"
        # Control characters (ANSI escapes that slipped through) never reach the PDF
        assert _to_latin1("\x1b[0;32m ✅ ok\x1b[0m") == "[0;32m [OK] ok[0m"
        assert "\x1b" not in _to_latin1("\x1b[0;32m ok")

    def test_polish_letters(self):
        assert _to_latin1("PŁYWANIE żółć") == "PLYWANIE zolc"

    def test_result_is_latin1_encodable(self):
        text = _to_latin1("🧬 genom_łukasza → wynik ≈ 5 — ß 🏃")
        text.encode('latin-1')  # must not raise
        assert "genom_lukasza -> wynik ~ 5 - ss" in text


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

        assert isinstance(pdf_bytes, bytes)
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

        assert isinstance(pdf_bytes, bytes)
        assert len(pdf_bytes) > 0
        assert pdf_bytes[:4] == b'%PDF'

    def test_generate_pdf_empty_results(self, sample_metadata):
        """Test PDF generation with empty results"""
        pdf_bytes = generate_pdf_report(
            vcf_filename="empty.vcf",
            results={},
            analyses_metadata=sample_metadata
        )

        assert isinstance(pdf_bytes, bytes)
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

        assert isinstance(pdf_bytes, bytes)
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

        assert isinstance(pdf_bytes, bytes)
        assert pdf_bytes[:4] == b'%PDF'

    def test_generate_pdf_findings_and_warnings(self, sample_metadata):
        """Three findings plus two warnings render as consecutive multi_cells"""
        results = {
            "fun_genetics": AnalysisResult(
                name="Fun Genetics",
                status=AnalysisStatus.COMPLETED,
                findings=["⚠️ a", "⚠️ b", "✅ No pathogenic variants found"],
                warnings=["⚠️ first warning", "⚠️ second warning"],
                scores={"PŁYWANIE": {"value": 3, "max": 5, "percentage": 60}},
            )
        }

        pdf_bytes = generate_pdf_report(
            vcf_filename="test.vcf",
            results=results,
            analyses_metadata=sample_metadata
        )

        assert isinstance(pdf_bytes, bytes)
        assert pdf_bytes[:4] == b'%PDF'

    def test_generate_pdf_non_latin1_filename(self, sample_results, sample_metadata):
        """A VCF file name outside latin-1 must not crash the title page"""
        pdf_bytes = generate_pdf_report(
            vcf_filename="genom_łukasza_🧬.vcf.gz",
            results=sample_results,
            analyses_metadata=sample_metadata
        )

        assert isinstance(pdf_bytes, bytes)
        assert pdf_bytes[:4] == b'%PDF'


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
