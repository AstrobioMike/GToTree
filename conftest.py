import os
from pathlib import Path

os.environ["COVERAGE_PROCESS_START"] = str(Path(__file__).resolve().parent / "pyproject.toml")
