"""
Shared pytest configuration for the Helixight python tests.

Evaluated before the test modules import `analysis`, so the environment the
frontend reads at import time points at a throwaway results directory and the
repository's scripts, and `frontend/` is importable.
"""

import os
import sys
import tempfile

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
FRONTEND_DIR = os.path.join(REPO_ROOT, 'frontend')

os.environ.setdefault('RESULTS_DIR', tempfile.mkdtemp(prefix='helixight-results-'))
os.environ.setdefault('SCRIPTS_DIR', os.path.join(REPO_ROOT, 'scripts'))
os.environ.setdefault('HELIXIGHT_WORK_DIR', os.path.join(os.environ['RESULTS_DIR'], '.work'))

if FRONTEND_DIR not in sys.path:
    sys.path.insert(0, FRONTEND_DIR)
