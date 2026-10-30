"""Open this file in Spyder and press F5. Edit settings.py for input/settings.

Every run creates a separate results folder. Source TIFFs are read-only. No
chat caches, prior masks, internet access, or bundled Codex runtime are needed.
RESULT is left in Spyder's variable explorer after a successful run.
"""
from __future__ import annotations

from pathlib import Path
import runpy
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))


def main(settings=None):
    if settings is None:
        # Reading this afresh avoids stale settings after repeated Spyder runs.
        settings = runpy.run_path(str(HERE / "settings.py"))["SETTINGS"]
    try:
        from phase_analysis.pipeline import run
    except ImportError as exc:
        raise ImportError(
            "A required Python package is unavailable. Install requirements.txt "
            "in the same environment as the Spyder kernel, then restart it. "
            f"Kernel interpreter: {sys.executable}. Original error: {exc}"
        ) from exc
    return run(settings)


if __name__ == "__main__":
    RESULT = main()
