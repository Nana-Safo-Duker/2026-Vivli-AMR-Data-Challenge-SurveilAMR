"""Run the full SurveilAMR pipeline end to end.

Executes, in order:
  1. run_analysis.py           — ATLAS full-cohort summary (chunked, ~1M rows)
  2. analyze_supplementary.py  — KEYSTONE / DREAM / GASAR summaries + figures
  3. analyze_spidaar.py        — SPIDAAR isolate + patient summaries + figures
  4. generate_figures.py       — ATLAS publication figures (fig1-fig6, fig10)
  5. eda_surveilamr.py         — fast sample-based exploratory figures

Run: python scripts/run_all.py
"""
import subprocess
import sys
import time
from pathlib import Path

SCRIPTS = [
    "run_analysis.py",
    "analyze_supplementary.py",
    "analyze_spidaar.py",
    "generate_figures.py",
    "eda_surveilamr.py",
]


def main() -> None:
    root = Path(__file__).resolve().parent
    for name in SCRIPTS:
        path = root / name
        print(f"\n{'=' * 70}\nRunning {name}\n{'=' * 70}")
        start = time.time()
        result = subprocess.run([sys.executable, str(path)])
        elapsed = time.time() - start
        if result.returncode != 0:
            print(f"{name} failed after {elapsed:.1f}s (exit code {result.returncode})")
            sys.exit(result.returncode)
        print(f"{name} completed in {elapsed:.1f}s")
    print("\nAll SurveilAMR analysis steps completed successfully.")


if __name__ == "__main__":
    main()
