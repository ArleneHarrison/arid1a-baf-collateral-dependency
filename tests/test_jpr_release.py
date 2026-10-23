"""The v1.1.0 primary estimand is distinct from historical co-dependency."""
from pathlib import Path
import runpy

def test_current_release():
    root=Path(__file__).resolve().parents[1]
    runpy.run_path(str(root/'current/validate_release.py'))['main']()
