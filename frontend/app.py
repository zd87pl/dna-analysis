"""
Helixight - DNA Analysis Toolkit
Streamlit Web Frontend
"""

import streamlit as st
import os
import tempfile
import shutil
from pathlib import Path
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as px

from analysis import (
    AnalysisRunner, AnalysisStatus, get_system_info,
    get_vcf_metrics, VCFMetrics, _result_cache
)
from pdf_report import generate_pdf_report

# Page configuration
st.set_page_config(
    page_title="Helixight - DNA Analysis",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        text-align: center;
        padding: 1rem 0;
    }
    .subtitle {
        text-align: center;
        color: #666;
        margin-bottom: 2rem;
    }
    .metric-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1rem;
        border-radius: 10px;
        color: white;
        text-align: center;
    }
    .warning-box {
        background-color: #fff3cd;
        border: 1px solid #ffc107;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .success-box {
        background-color: #d4edda;
        border: 1px solid #28a745;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .stProgress > div > div > div > div {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
    }
</style>
""", unsafe_allow_html=True)

# Constants
DATA_DIR = os.environ.get("DATA_DIR", "/data")
RESULTS_DIR = os.environ.get("RESULTS_DIR", "/app/results")

# Ensure directories exist
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


def init_session_state():
    """Initialize session state variables"""
    if 'vcf_path' not in st.session_state:
        st.session_state.vcf_path = None
    if 'vcf_info' not in st.session_state:
        st.session_state.vcf_info = None
    if 'vcf_metrics' not in st.session_state:
        st.session_state.vcf_metrics = None
    if 'results' not in st.session_state:
        st.session_state.results = {}
    if 'language' not in st.session_state:
        st.session_state.language = 'en'
    if 'use_cache' not in st.session_state:
        st.session_state.use_cache = True


def render_header():
    """Render the main header"""
    st.markdown('<h1 class="main-header">🧬 HELIXIGHT</h1>', unsafe_allow_html=True)
    st.markdown('<p class="subtitle">Open Source Genetic Analysis Toolkit</p>', unsafe_allow_html=True)


def render_sidebar():
    """Render the sidebar with file management and settings"""
    with st.sidebar:
        st.markdown("### 📁 Data Management")

        # File upload - accept vcf and gz extensions (covers .vcf.gz files)
        uploaded_file = st.file_uploader(
            "Upload VCF File",
            type=['vcf', 'gz'],
            help="Upload a VCF (.vcf) or compressed VCF (.vcf.gz) file"
        )

        if uploaded_file:
            # Validate filename pattern (must be .vcf or .vcf.gz)
            filename = uploaded_file.name
            if not (filename.endswith('.vcf') or filename.endswith('.vcf.gz')):
                st.error("❌ File must be a VCF file (.vcf or .vcf.gz)")
            else:
                # Save uploaded file
                file_path = os.path.join(DATA_DIR, filename)
                with open(file_path, 'wb') as f:
                    f.write(uploaded_file.getbuffer())

                # Validate
                is_valid, message = AnalysisRunner.validate_vcf(file_path)
                if is_valid:
                    st.session_state.vcf_path = file_path
                    st.session_state.vcf_info = message
                    # Calculate VCF metrics
                    with st.spinner("Calculating VCF metrics..."):
                        st.session_state.vcf_metrics = get_vcf_metrics(file_path)
                    st.success(f"✅ {message}")
                else:
                    st.error(f"❌ {message}")
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass  # Ignore removal errors

        # Or select existing file
        st.markdown("---")
        st.markdown("**Or select existing file:**")

        # Find VCF files (both .vcf and .vcf.gz)
        vcf_files = sorted(
            list(Path(DATA_DIR).glob("*.vcf")) + list(Path(DATA_DIR).glob("*.vcf.gz")),
            key=lambda f: f.stat().st_mtime,
            reverse=True
        )
        if vcf_files:
            file_options = ["Select a file..."] + [f.name for f in vcf_files]
            selected = st.selectbox("Available VCF files", file_options)

            if selected != "Select a file...":
                file_path = os.path.join(DATA_DIR, selected)
                if st.button("Load Selected File"):
                    is_valid, message = AnalysisRunner.validate_vcf(file_path)
                    if is_valid:
                        st.session_state.vcf_path = file_path
                        st.session_state.vcf_info = message
                        # Calculate VCF metrics
                        st.session_state.vcf_metrics = get_vcf_metrics(file_path)
                        st.success(f"✅ Loaded: {selected}")
                        st.rerun()
                    else:
                        st.error(f"❌ {message}")
        else:
            st.info("No VCF files found in /data directory")

        # Current file status
        st.markdown("---")
        st.markdown("### 📊 Current File")
        if st.session_state.vcf_path:
            st.success(f"**{os.path.basename(st.session_state.vcf_path)}**")
            st.caption(st.session_state.vcf_info)

            # Display VCF metrics if available
            if st.session_state.vcf_metrics:
                metrics = st.session_state.vcf_metrics
                with st.expander("📈 VCF Quality Metrics", expanded=False):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.metric("Total Variants", f"{metrics.total_variants:,}")
                        st.metric("SNPs", f"{metrics.snps:,}")
                        st.metric("Indels", f"{metrics.indels:,}")
                    with col2:
                        st.metric("Heterozygous", f"{metrics.heterozygous:,}")
                        st.metric("Homozygous Alt", f"{metrics.homozygous_alt:,}")
                        if metrics.titv_ratio > 0:
                            st.metric("Ti/Tv Ratio", f"{metrics.titv_ratio:.2f}")
                    st.caption(f"File size: {metrics.file_size_mb:.2f} MB")
        else:
            st.warning("No file loaded")

        # Settings
        st.markdown("---")
        st.markdown("### ⚙️ Settings")
        st.session_state.language = st.selectbox(
            "Language",
            ["en", "pl"],
            format_func=lambda x: "English" if x == "en" else "Polski"
        )

        st.session_state.use_cache = st.checkbox(
            "Use result caching",
            value=st.session_state.use_cache,
            help="Cache analysis results to avoid re-running on the same file"
        )

        if st.button("Clear Cache", help="Clear all cached analysis results"):
            _result_cache.clear()
            st.success("Cache cleared!")

        # About
        st.markdown("---")
        st.markdown("### ℹ️ About")
        st.caption("Version 1.1.0")
        st.caption("[GitHub](https://github.com/helixight/helixight-oss)")

        # System info
        with st.expander("System Info"):
            info = get_system_info()
            for tool, version in info.items():
                st.text(f"{tool}: {version}")


def render_analysis_selection():
    """Render analysis selection interface"""
    st.markdown("## 🔬 Select Analyses")

    if not st.session_state.vcf_path:
        st.warning("⚠️ Please upload or select a VCF file first")
        return []

    selected_analyses = []
    categories = AnalysisRunner.get_analyses_by_category()

    # Category display names
    category_names = {
        "fitness": "🏃 Fitness & Athletics",
        "health": "🏥 Health & Clinical",
        "comprehensive": "🧬 Comprehensive",
        "fun": "🎲 Fun & Ancestry",
        "advanced": "🔬 Advanced"
    }

    cols = st.columns(2)

    for idx, (category, analyses) in enumerate(categories.items()):
        col = cols[idx % 2]
        with col:
            st.markdown(f"### {category_names.get(category, category.title())}")

            for analysis in analyses:
                key = f"select_{analysis['id']}"
                if st.checkbox(
                    f"{analysis['icon']} {analysis['name']}",
                    key=key,
                    help=analysis['description']
                ):
                    selected_analyses.append(analysis['id'])

    return selected_analyses


def render_progress(analysis_id: str, progress: float, message: str):
    """Update progress display"""
    if 'progress_bar' in st.session_state:
        st.session_state.progress_bar.progress(progress, text=message)


def run_analyses(selected: list):
    """Run selected analyses and display results"""
    if not selected:
        st.warning("Please select at least one analysis")
        return

    st.markdown("## 📈 Running Analyses")

    # Progress tracking
    progress_bar = st.progress(0, text="Initializing...")
    status_text = st.empty()

    # Use caching based on user preference
    runner = AnalysisRunner(
        st.session_state.vcf_path,
        use_cache=st.session_state.use_cache
    )

    results = {}
    total = len(selected)
    cached_count = 0

    for i, analysis_id in enumerate(selected):
        analysis_info = AnalysisRunner.ANALYSES.get(analysis_id, {})
        analysis_name = analysis_info.get('name', analysis_id)

        progress = i / total
        progress_bar.progress(progress, text=f"Running {analysis_name}...")
        status_text.info(f"🔄 Processing: {analysis_name}")

        result = runner.run_analysis(analysis_id)
        results[analysis_id] = result

        if result.status == AnalysisStatus.COMPLETED:
            # Check if this was loaded from cache
            if st.session_state.use_cache and analysis_id in runner.results:
                cached_count += 1
            status_text.success(f"✅ Completed: {analysis_name}")
        else:
            status_text.error(f"❌ Failed: {analysis_name}")

    completion_text = "All analyses complete!"
    if cached_count > 0:
        completion_text += f" ({cached_count} loaded from cache)"
    progress_bar.progress(1.0, text=completion_text)
    st.session_state.results = results

    return results


def render_results():
    """Display analysis results"""
    if not st.session_state.results:
        return

    st.markdown("## 📊 Results")

    # Summary metrics
    total = len(st.session_state.results)
    completed = sum(1 for r in st.session_state.results.values()
                   if r.status == AnalysisStatus.COMPLETED)
    findings = sum(len(r.findings) for r in st.session_state.results.values())
    warnings = sum(len(r.warnings) for r in st.session_state.results.values())

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Analyses Run", total)
    with col2:
        st.metric("Completed", completed, delta=None if completed == total else f"-{total-completed}")
    with col3:
        st.metric("Findings", findings)
    with col4:
        st.metric("Warnings", warnings, delta_color="inverse" if warnings > 0 else "off")

    st.markdown("---")

    # Detailed results in tabs - filter to valid analysis IDs
    valid_results = {
        aid: result for aid, result in st.session_state.results.items()
        if aid in AnalysisRunner.ANALYSES
    }

    if valid_results:
        tabs = st.tabs([
            AnalysisRunner.ANALYSES[aid]['icon'] + " " + AnalysisRunner.ANALYSES[aid]['name']
            for aid in valid_results.keys()
        ])

        for tab, (analysis_id, result) in zip(tabs, valid_results.items()):
            with tab:
                render_single_result(analysis_id, result)


def render_single_result(analysis_id: str, result):
    """Render results for a single analysis"""
    analysis_info = AnalysisRunner.ANALYSES.get(analysis_id, {})

    # Status indicator
    if result.status == AnalysisStatus.COMPLETED:
        st.success(f"✅ Analysis completed successfully")
    else:
        st.error(f"❌ Analysis failed: {result.error}")
        return

    col1, col2 = st.columns([2, 1])

    with col1:
        # Scores visualization
        if result.scores:
            st.markdown("### 📊 Scores")

            # Create bar chart for scores
            score_names = []
            score_values = []
            score_colors = []

            for name, data in result.scores.items():
                if 'percentage' in data:
                    score_names.append(name)
                    pct = data['percentage']
                    score_values.append(pct)
                    # Color based on value
                    if pct >= 70:
                        score_colors.append('#28a745')
                    elif pct >= 40:
                        score_colors.append('#ffc107')
                    else:
                        score_colors.append('#dc3545')

            if score_names:
                fig = go.Figure(go.Bar(
                    x=score_values,
                    y=score_names,
                    orientation='h',
                    marker_color=score_colors,
                    text=[f"{v}%" for v in score_values],
                    textposition='auto',
                ))
                fig.update_layout(
                    height=max(200, len(score_names) * 50),
                    margin=dict(l=0, r=0, t=0, b=0),
                    xaxis_range=[0, 100],
                    xaxis_title="Score (%)",
                    showlegend=False
                )
                st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Findings
        if result.findings:
            st.markdown("### 🔍 Key Findings")
            for finding in result.findings[:10]:  # Limit display
                # Clean up the finding text
                clean_finding = finding.strip()
                if '✅' in clean_finding:
                    st.success(clean_finding)
                elif '🔴' in clean_finding or 'FOUND' in clean_finding:
                    st.error(clean_finding)
                elif '⚠️' in clean_finding:
                    st.warning(clean_finding)
                else:
                    st.info(clean_finding)

        # Warnings
        if result.warnings:
            st.markdown("### ⚠️ Warnings")
            for warning in result.warnings[:5]:
                st.warning(warning)

    # Raw output (collapsible)
    with st.expander("📝 Raw Output"):
        st.code(result.output, language="text")

    if result.error:
        with st.expander("❌ Errors"):
            st.code(result.error, language="text")


def render_report_download():
    """Render report download section"""
    if not st.session_state.results:
        return

    st.markdown("## 📄 Generate Report")

    col1, col2 = st.columns(2)

    with col1:
        # Generate report and provide direct download button (not nested)
        report = generate_text_report()
        st.download_button(
            "📥 Download Text Report",
            report,
            file_name=f"helixight_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
            mime="text/plain",
            use_container_width=True
        )

    with col2:
        # Generate PDF report
        vcf_name = os.path.basename(st.session_state.vcf_path) if st.session_state.vcf_path else "Unknown"

        # Convert VCFMetrics to dict for PDF generation
        vcf_metrics_dict = None
        if st.session_state.vcf_metrics:
            m = st.session_state.vcf_metrics
            vcf_metrics_dict = {
                'total_variants': m.total_variants,
                'snps': m.snps,
                'indels': m.indels,
                'heterozygous': m.heterozygous,
                'homozygous_alt': m.homozygous_alt,
                'titv_ratio': m.titv_ratio,
                'file_size_mb': m.file_size_mb,
            }

        pdf_bytes = generate_pdf_report(
            vcf_filename=vcf_name,
            results=st.session_state.results,
            analyses_metadata=AnalysisRunner.ANALYSES,
            vcf_metrics=vcf_metrics_dict
        )

        st.download_button(
            "📥 Download PDF Report",
            pdf_bytes,
            file_name=f"helixight_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
            mime="application/pdf",
            use_container_width=True
        )


def generate_text_report() -> str:
    """Generate a text report from results"""
    vcf_name = os.path.basename(st.session_state.vcf_path) if st.session_state.vcf_path else "Unknown"
    lines = [
        "=" * 80,
        "                    HELIXIGHT GENETIC ANALYSIS REPORT",
        f"                    Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 80,
        "",
        f"VCF File: {vcf_name}",
        f"Analyses Run: {len(st.session_state.results)}",
        "",
        "=" * 80,
    ]

    for analysis_id, result in st.session_state.results.items():
        analysis_info = AnalysisRunner.ANALYSES.get(analysis_id, {})
        lines.extend([
            "",
            f"{'─' * 80}",
            f"  {analysis_info.get('icon', '')} {result.name}",
            f"  Status: {result.status.value}",
            f"{'─' * 80}",
        ])

        if result.scores:
            lines.append("\n  SCORES:")
            for name, data in result.scores.items():
                if 'percentage' in data:
                    pct = data['percentage']
                    bar = '█' * (pct // 10) + '░' * (10 - pct // 10)
                    lines.append(f"    {name}: {bar} {pct}%")

        if result.findings:
            lines.append("\n  KEY FINDINGS:")
            for finding in result.findings:
                lines.append(f"    • {finding}")

        if result.warnings:
            lines.append("\n  WARNINGS:")
            for warning in result.warnings:
                lines.append(f"    ⚠ {warning}")

    lines.extend([
        "",
        "=" * 80,
        "                              END OF REPORT",
        "=" * 80,
        "",
        "DISCLAIMER: This report is for educational and research purposes only.",
        "Results should be interpreted by a qualified genetic counselor or physician.",
        "This is NOT a medical diagnostic tool.",
    ])

    return "\n".join(lines)


def render_disclaimer():
    """Render medical disclaimer"""
    st.markdown("---")
    st.markdown("""
    <div class="warning-box">
    <strong>⚠️ Important Disclaimer</strong><br>
    Helixight is an educational and research tool. Results are NOT intended for medical diagnosis,
    treatment decisions, or clinical use. Always consult a qualified healthcare professional
    or genetic counselor for interpretation of genetic data. This software is NOT FDA/EMA approved.
    </div>
    """, unsafe_allow_html=True)


def main():
    """Main application"""
    init_session_state()
    render_header()
    render_sidebar()

    # Main content
    tab1, tab2, tab3 = st.tabs(["🔬 Analysis", "📊 Results", "📚 Help"])

    with tab1:
        selected = render_analysis_selection()

        st.markdown("---")

        if st.button("▶️ Run Selected Analyses", type="primary",
                    disabled=not st.session_state.vcf_path or not selected,
                    use_container_width=True):
            with st.spinner("Running analyses..."):
                run_analyses(selected)
            st.rerun()

    with tab2:
        render_results()
        render_report_download()

    with tab3:
        render_help()

    render_disclaimer()


def render_help():
    """Render help section"""
    st.markdown("## 📚 Help & Documentation")

    with st.expander("🚀 Quick Start", expanded=True):
        st.markdown("""
        1. **Upload your VCF file** using the sidebar
        2. **Select analyses** you want to run
        3. **Click "Run Selected Analyses"**
        4. **View results** in the Results tab
        5. **Download report** for your records
        """)

    with st.expander("📁 Supported File Formats"):
        st.markdown("""
        - **VCF** (.vcf) - Variant Call Format
        - **Compressed VCF** (.vcf.gz) - gzip-compressed VCF

        VCF files should be aligned to **GRCh38** (hg38) reference genome.
        """)

    with st.expander("🔬 Analysis Descriptions"):
        for aid, info in AnalysisRunner.ANALYSES.items():
            st.markdown(f"**{info['icon']} {info['name']}**")
            st.caption(info['description'])
            st.markdown("---")

    with st.expander("🔒 Privacy & Security"):
        st.markdown("""
        - **100% Local Processing** - Your genetic data never leaves this container
        - **No Cloud Upload** - All analysis runs on your local machine
        - **No Data Collection** - We don't collect any usage data
        - **Open Source** - Audit the code yourself
        """)

    with st.expander("❓ FAQ"):
        st.markdown("""
        **Q: Where can I get a VCF file?**
        A: VCF files can be obtained from consumer genetic testing services
        (23andMe, AncestryDNA with raw data export), or from whole genome
        sequencing services.

        **Q: Is this medically accurate?**
        A: This tool is for educational purposes only. Results should be
        interpreted by healthcare professionals.

        **Q: Can I run this offline?**
        A: Yes! Once the Docker container is downloaded, no internet connection
        is required for analysis.
        """)


if __name__ == "__main__":
    main()
