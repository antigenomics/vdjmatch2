from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def _binary_path() -> Path:
    package_dir = Path(__file__).resolve().parent
    candidates = [package_dir / "vdjmatch2"]
    if os.name == "nt":
        candidates.insert(0, package_dir / "vdjmatch2.exe")
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Bundled vdjmatch2 binary was not found inside the installed package."
    )


def main() -> int:
    binary = _binary_path()
    completed = subprocess.run([str(binary), *sys.argv[1:]])
    return int(completed.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
