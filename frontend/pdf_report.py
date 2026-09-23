"""
Helixight PDF Report Generator
Generates professional PDF reports from analysis results
"""

from fpdf import FPDF
from fpdf.enums import XPos, YPos
from datetime import datetime
from typing import Dict, Any
import unicodedata


# Characters the core (latin-1) fonts cannot render, mapped to ASCII stand-ins.
_MARKER_MAP = {
    '✅': '[OK]',
    '🔴': '[!]',
    '⚠': '[!]',
    '❌': '[X]',
    '•': '-',
    '→': '->',
    '≈': '~',
    '–': '-',
    '—': '-',
}

# Letters that NFKD cannot decompose into a base letter + combining mark.
_LETTER_MAP = {
    'ł': 'l',
    'Ł': 'L',
    'ß': 'ss',
}


def _to_latin1(text: str) -> str:
    """
    Reduce arbitrary text to something the core latin-1 PDF fonts can draw.

    Markers (emoji) become ASCII tags, accented letters lose their accents,
    and anything that still cannot be encoded is dropped.
    """
    if not text:
        return text
    text = text.replace('️', '')  # variation selector after emoji
    for src, dst in _MARKER_MAP.items():
        text = text.replace(src, dst)
    for src, dst in _LETTER_MAP.items():
        text = text.replace(src, dst)
    decomposed = unicodedata.normalize('NFKD', text)
    kept = []
    for ch in decomposed:
        cat = unicodedata.category(ch)
        if cat == 'Mn':          # combining marks (accents)
            continue
        if cat in ('So', 'Sk', 'Cs', 'Co'):  # symbols, surrogates, private use
            continue
        if cat == 'Cc' and ch not in '\n\t':   # control characters (e.g. ESC)
            continue
        kept.append(ch)
    return ''.join(kept).encode('latin-1', 'ignore').decode('latin-1')


class HelixightReport(FPDF):
    """Custom PDF class for Helixight reports"""

    def __init__(self):
        super().__init__()
        self.set_auto_page_break(auto=True, margin=15)

    def normalize_text(self, text):
        """Sanitise text for the core fonts before fpdf2 encodes it"""
        if self.current_font and getattr(self.current_font, 'type', 'core') == 'core':
            text = _to_latin1(text)
        return super().normalize_text(text)

    def header(self):
        """Page header"""
        self.set_font('Helvetica', 'B', 12)
        self.set_text_color(102, 126, 234)  # Brand purple
        self.cell(0, 10, 'HELIXIGHT', 0, new_x=XPos.RIGHT, new_y=YPos.TOP, align='L')
        self.set_font('Helvetica', '', 8)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, 'Genetic Analysis Report', 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='R')
        self.line(10, 20, 200, 20)
        self.ln(5)

    def footer(self):
        """Page footer"""
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)
        self.set_text_color(128, 128, 128)
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', 0, new_x=XPos.RIGHT, new_y=YPos.TOP, align='C')

    def chapter_title(self, title: str, icon: str = ""):
        """Add a chapter/section title"""
        self.set_font('Helvetica', 'B', 14)
        self.set_text_color(51, 51, 51)
        self.set_fill_color(240, 240, 250)
        display_title = f"{icon} {title}" if icon else title
        self.cell(0, 10, display_title, 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='L', fill=True)
        self.ln(2)

    def add_score_bar(self, name: str, percentage: int, max_width: int = 100):
        """Add a visual score bar"""
        self.set_font('Helvetica', '', 10)
        self.set_text_color(51, 51, 51)

        try:
            percentage = int(percentage)
        except (TypeError, ValueError):
            percentage = 0
        percentage = max(0, min(100, percentage))
        name = str(name)[:35]

        # Name
        self.cell(60, 8, name, 0, new_x=XPos.RIGHT, new_y=YPos.TOP)

        # Bar background
        bar_x = self.get_x()
        bar_y = self.get_y() + 2
        bar_height = 4

        # Gray background
        self.set_fill_color(220, 220, 220)
        self.rect(bar_x, bar_y, max_width, bar_height, 'F')

        # Colored fill based on percentage
        if percentage >= 70:
            self.set_fill_color(40, 167, 69)  # Green
        elif percentage >= 40:
            self.set_fill_color(255, 193, 7)  # Yellow
        else:
            self.set_fill_color(220, 53, 69)  # Red

        fill_width = (percentage / 100) * max_width
        self.rect(bar_x, bar_y, fill_width, bar_height, 'F')

        # Percentage text
        self.set_x(bar_x + max_width + 5)
        self.cell(20, 8, f"{percentage}%", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    def add_finding(self, finding: str, finding_type: str = "info"):
        """Add a finding with appropriate styling"""
        self.set_font('Helvetica', '', 9)

        if finding_type == "positive" or "no pathogenic" in finding.lower():
            self.set_text_color(40, 167, 69)  # Green
            prefix = "[OK] "
        elif finding_type == "warning" or "FOUND" in finding or '🔴' in finding:
            self.set_text_color(220, 53, 69)  # Red
            prefix = "[!] "
        else:
            self.set_text_color(51, 51, 51)
            prefix = "- "

        # Emoji and non-latin-1 characters are sanitised by normalize_text()
        clean_finding = finding.strip()
        self.multi_cell(0, 5, prefix + clean_finding, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_text_color(51, 51, 51)


def _classify_finding(finding: str) -> str:
    """Classify a raw finding line for PDF styling (evaluated before sanitising)"""
    if '✅' in finding or 'no pathogenic' in finding.lower():
        return "positive"
    if '🔴' in finding or 'FOUND' in finding:
        return "warning"
    return "info"


def generate_pdf_report(
    vcf_filename: str,
    results: Dict[str, Any],
    analyses_metadata: Dict[str, Dict],
    vcf_metrics: Dict[str, Any] = None
) -> bytes:
    """
    Generate a PDF report from analysis results

    Args:
        vcf_filename: Name of the analyzed VCF file
        results: Dictionary of analysis_id -> AnalysisResult
        analyses_metadata: Metadata about analyses (icons, names, etc.)
        vcf_metrics: Optional VCF quality metrics

    Returns:
        PDF content as bytes
    """
    pdf = HelixightReport()
    pdf.alias_nb_pages()
    pdf.add_page()

    # Title page content
    pdf.set_font('Helvetica', 'B', 24)
    pdf.set_text_color(102, 126, 234)
    pdf.ln(20)
    pdf.cell(0, 15, 'HELIXIGHT', 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')

    pdf.set_font('Helvetica', '', 14)
    pdf.set_text_color(100, 100, 100)
    pdf.cell(0, 10, 'Genetic Analysis Report', 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')

    pdf.ln(20)

    # Report metadata
    pdf.set_font('Helvetica', '', 11)
    pdf.set_text_color(51, 51, 51)
    pdf.cell(0, 8, f"File: {vcf_filename}", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')
    pdf.cell(0, 8, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", 0,
             new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')
    pdf.cell(0, 8, f"Analyses Run: {len(results)}", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT, align='C')

    # VCF Quality Metrics section
    if vcf_metrics:
        pdf.ln(15)
        pdf.chapter_title("VCF Quality Metrics")
        pdf.set_font('Helvetica', '', 10)

        metrics_table = [
            ("Total Variants", f"{vcf_metrics.get('total_variants', 'N/A'):,}"),
            ("SNPs", f"{vcf_metrics.get('snps', 'N/A'):,}"),
            ("Indels", f"{vcf_metrics.get('indels', 'N/A'):,}"),
            ("Heterozygous", f"{vcf_metrics.get('heterozygous', 'N/A'):,}"),
            ("Homozygous Alt", f"{vcf_metrics.get('homozygous_alt', 'N/A'):,}"),
            ("Ti/Tv Ratio", f"{vcf_metrics.get('titv_ratio', 'N/A'):.2f}" if vcf_metrics.get('titv_ratio') else "N/A"),
        ]

        for label, value in metrics_table:
            pdf.cell(80, 7, label, 0, new_x=XPos.RIGHT, new_y=YPos.TOP)
            pdf.cell(0, 7, str(value), 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    # Summary section
    pdf.add_page()
    pdf.chapter_title("Summary")

    total_analyses = len(results)
    completed = sum(1 for r in results.values() if r.status.value == "completed")
    total_findings = sum(len(r.findings) for r in results.values())
    total_warnings = sum(len(r.warnings) for r in results.values())

    pdf.set_font('Helvetica', '', 11)
    summary_data = [
        ("Analyses Completed", f"{completed}/{total_analyses}"),
        ("Key Findings", str(total_findings)),
        ("Warnings", str(total_warnings)),
    ]

    for label, value in summary_data:
        pdf.cell(80, 8, label, 0, new_x=XPos.RIGHT, new_y=YPos.TOP)
        pdf.set_font('Helvetica', 'B', 11)
        pdf.cell(0, 8, value, 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_font('Helvetica', '', 11)

    # Individual analysis results
    for analysis_id, result in results.items():
        metadata = analyses_metadata.get(analysis_id, {})
        icon = metadata.get('icon', '')
        name = result.name

        pdf.add_page()
        # Clean icon for PDF (remove emoji)
        pdf.chapter_title(name)

        # Status
        pdf.set_font('Helvetica', 'B', 10)
        if result.status.value == "completed":
            pdf.set_text_color(40, 167, 69)
            pdf.cell(0, 8, "Status: COMPLETED", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        else:
            pdf.set_text_color(220, 53, 69)
            pdf.cell(0, 8, f"Status: {result.status.value.upper()}", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        pdf.set_text_color(51, 51, 51)

        # Scores
        if result.scores:
            pdf.ln(5)
            pdf.set_font('Helvetica', 'B', 11)
            pdf.cell(0, 8, "Scores", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

            for score_name, score_data in result.scores.items():
                if 'percentage' in score_data:
                    pdf.add_score_bar(score_name, score_data['percentage'])

        # Findings
        if result.findings:
            pdf.ln(5)
            pdf.set_font('Helvetica', 'B', 11)
            pdf.cell(0, 8, "Key Findings", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)

            for finding in result.findings[:15]:  # Limit to prevent overflow
                pdf.add_finding(finding, _classify_finding(finding))

        # Warnings
        if result.warnings:
            pdf.ln(5)
            pdf.set_font('Helvetica', 'B', 11)
            pdf.set_text_color(255, 193, 7)
            pdf.cell(0, 8, "Warnings", 0, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.set_text_color(51, 51, 51)

            for warning in result.warnings[:10]:
                pdf.add_finding(warning, "warning")

    # Disclaimer page
    pdf.add_page()
    pdf.chapter_title("Important Disclaimer")

    pdf.set_font('Helvetica', '', 10)
    pdf.set_text_color(220, 53, 69)

    disclaimer_text = """
This report is generated by Helixight, an open-source genetic analysis toolkit provided for EDUCATIONAL AND RESEARCH PURPOSES ONLY.

THIS IS NOT A MEDICAL DIAGNOSTIC TOOL.

Key Points:
- This software is NOT clinically validated
- Results may contain errors, false positives, or false negatives
- Do NOT make any health decisions based on this report
- Do NOT use this as a substitute for professional genetic counseling
- This software is NOT approved by FDA, EMA, or any regulatory body

If you discover concerning results:
- Do NOT panic - these are not clinical-grade results
- Consult a certified genetic counselor or medical geneticist
- Get proper clinical testing through accredited laboratories

Genetic predisposition does not equal destiny. Most traits and diseases are influenced by lifestyle, environment, and complex interactions between hundreds of genes.

By reading this report, you acknowledge that this is for educational purposes only and that the authors bear no responsibility for any actions taken based on these results.
"""

    pdf.multi_cell(0, 5, disclaimer_text.strip(), new_x=XPos.LMARGIN, new_y=YPos.NEXT)

    # Get PDF as bytes (fpdf2 returns a bytearray; st.download_button needs bytes)
    return bytes(pdf.output())
