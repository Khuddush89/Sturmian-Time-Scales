"""Keep the original independent 535-assertion audit executable in CI."""

import json

import pytest

from sturmian.verification import run


@pytest.mark.slow
def test_complete_reference_verification(tmp_path):
    report = run(tmp_path / "run")
    assert report["passed"] is True
    assert report["check_count"] == 535
    assert all(record["passed"] for record in report["checks"])
    written = json.loads((tmp_path / "run/reports/verification_results.json").read_text())
    assert written == report
    assert len(list((tmp_path / "run/data").glob("*.csv"))) == 13
