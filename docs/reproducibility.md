# 🧪 Reproducibility

## Reference environment

The prepared snapshot was checked with Python 3.12 and the versions in
`requirements.txt`. CI also runs Python 3.11 and 3.13. Only the local environment
has been executed before publishing; GitHub checks run after the first push.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
ruff check .
ruff format --check .
pytest -q
sturmian verify --output outputs/reference
sturmian figures --output outputs/reference
```

## Output layout

| Path relative to the output directory | Contents |
| :--- | :--- |
| `data/*.csv` | 13 verification tables and 7 figure datasets |
| `reports/verification_results.json` | Every assertion, residual, tolerance, and run environment |
| `figures/*.png` | Seven recreated research figures |

The committed `reports/verification_results.json`, `data/`, and PNG figures
are the reference snapshot. New output directories are ignored by git.
Verification must precede plotting in the same directory.

## What the checks cover

- Exact Fraction arithmetic for gauge states, Green matrices, and polynomial identities.
- Sixteen seeded variable-gain grids with independently evaluated identities.
- Spectral residuals, mass orthogonality, nodal counts, and order comparisons.
- Pulse transfer roots checked against independent 80-digit Decimal calculations.
- Continuum and hybrid refinement, and direct/reduced variable-gain state agreement.
- Generalized-zero counts, near-threshold Euler cases, and Riccati sign loss.

The random seed is `20261001`. The longest streaming examples contain one
million nodes. Last floating-point digits and image metadata may differ across
platforms; compare using recorded numerical tolerances.

## CI artifacts

The `CI` workflow lints, tests, builds a wheel, and runs the installed wheel's CLI.
Download `sturmian-reproduction` from the workflow run for fresh CSVs, PNGs,
logs, JSON evidence, and built distributions.

The complete verification test uses a temporary directory, so a normal `pytest`
run does not overwrite the committed reference snapshot.
