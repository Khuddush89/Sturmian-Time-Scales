"""Exercise installed entry points and useful failure messages."""

import subprocess
import sys

from sturmian import __version__


def test_module_reports_package_version_from_another_directory(tmp_path):
    completed = subprocess.run(
        [sys.executable, "-m", "sturmian", "--version"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )
    assert completed.stdout.strip() == __version__


def test_figures_require_a_verification_table(tmp_path):
    completed = subprocess.run(
        [sys.executable, "-m", "sturmian", "figures", "--output", str(tmp_path / "run")],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 1
    assert "verify" in completed.stderr
    assert not list((tmp_path / "run").rglob("*.png"))
