<div align="center">

<img src="assets/hero.svg" alt="Sturmian Time Scales — exact grids, clear spectra, reproducible research" width="100%" />

<a href="https://github.com/Khuddush89/Sturmian-Time-Scales">
  <img src="https://readme-typing-svg.demolab.com?font=Fira+Code&amp;weight=500&amp;size=20&amp;duration=3400&amp;pause=1300&amp;color=06B6D4&amp;center=true&amp;vCenter=true&amp;width=800&amp;height=55&amp;lines=Gauge+reduction+on+time+scales;Exact+identities+%E2%80%A2+spectra+%E2%80%A2+zero+counts;Reproduce+the+computations+with+Python" alt="Animated heading: gauge reduction, spectra, and zero counts" width="800" />
</a>

[![CI](https://img.shields.io/github/actions/workflow/status/Khuddush89/Sturmian-Time-Scales/ci.yml?branch=main&style=for-the-badge&logo=githubactions&logoColor=white&label=CI)](https://github.com/Khuddush89/Sturmian-Time-Scales/actions/workflows/ci.yml)
[![Version](https://img.shields.io/badge/version-0.1.0-8B5CF6?style=for-the-badge&logo=python&logoColor=white)](CHANGELOG.md)
[![License](https://img.shields.io/badge/license-MIT-06B6D4?style=for-the-badge&logo=opensourceinitiative&logoColor=white)](LICENSE)
[![PRs welcome](https://img.shields.io/badge/PRs-welcome-EC4899?style=for-the-badge&logo=git&logoColor=white)](CONTRIBUTING.md)
[![Stars](https://img.shields.io/github/stars/Khuddush89/Sturmian-Time-Scales?style=for-the-badge&logo=github&logoColor=white&color=8B5CF6)](https://github.com/Khuddush89/Sturmian-Time-Scales/stargazers)
[![Forks](https://img.shields.io/github/forks/Khuddush89/Sturmian-Time-Scales?style=for-the-badge&logo=github&logoColor=white&color=06B6D4)](https://github.com/Khuddush89/Sturmian-Time-Scales/forks)

**A reproducible numerical supplement for proportional conformable dynamic equations.**

[Quick start](#-quick-start) · [Figures](#-figures-and-demos) · [Mathematics](docs/mathematics.md) · [Verification](docs/reproducibility.md)

</div>

<img src="assets/divider.svg" alt="" width="100%" />

## 🚀 Overview

Explore Sturmian spectra and generalized zeros on harmonic, square-root,
alternating, pulse, and hybrid time scales. The code pairs exact finite identities
with numerical experiments and seven reproducible figures.

The computations accompany *Gauge Reduction and Sturmian Consequences for
Proportional Conformable Dynamic Equations on Time Scales*.
See [provenance](docs/provenance.md) for the scope of this supplement.

### 🧭 Contents

- [✨ Features](#-features)
- [🛠️ Tech stack](#-tech-stack)
- [📦 Quick start](#-quick-start)
- [📸 Figures and demos](#-figures-and-demos)
- [🧪 Verification](#-verification)
- [🗂️ Repository map](#-repository-map)
- [🤝 Contributing](#-contributing)
- [📈 Community and star history](#-community-and-star-history)
- [📜 License and citation](#-license-and-citation)

## ✨ Features

| | Capability | Evidence / output |
| :---: | :--- | :--- |
| 🔮 | Positive gauge reduction with variable gains | Exact rational examples and seeded identity checks |
| 🎼 | Finite Dirichlet spectra and eigenfunctions | Symmetric tridiagonal pencils and mass orthogonality |
| 🔢 | Generalized-zero counts | Nodal zeros and strict crossings across contained gaps |
| 🌉 | Dense-block and jump transfer | Pulse thresholds and hybrid state comparisons |
| 📐 | Independent high-precision checks | 80-digit Decimal calculations |
| 📊 | Reproducible research figures | Seven PNGs with accompanying CSV data |
| 🧪 | Explicit verification evidence | 535 assertions with recorded residuals and tolerances |

## 🛠️ Tech stack

[![Stack: Python, Git, GitHub, GitHub Actions](https://skillicons.dev/icons?i=py,git,github,githubactions&theme=dark)](pyproject.toml)

**Python 3.11+ · NumPy · SciPy · Matplotlib · pytest · Ruff**

Exact arithmetic uses Python's `fractions` and `decimal` modules.
The [reference dependency versions](requirements.txt) reproduce the bundled run.

## 📦 Quick start

```bash
git clone https://github.com/Khuddush89/Sturmian-Time-Scales.git
cd Sturmian-Time-Scales
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .
sturmian verify --output outputs/reference
sturmian figures --output outputs/reference
```

On Windows PowerShell, activate with `.venv\Scripts\Activate.ps1`.
Both commands also work as `python -m sturmian verify` and `python -m sturmian figures`.

```python
from sturmian.model import harmonic, spectrum

grid = harmonic(15)
eigenvalues = spectrum(grid, alpha=0.75)
print(eigenvalues[:3])
```

Runs write CSVs, PNGs, and a JSON verification report under the chosen output
directory. The committed `data/` and `assets/screenshots/` provide reference results.

## 📸 Figures and demos

| Eigenfunctions on two domains | Spectral ratios |
| :---: | :---: |
| ![Harmonic nodes and exact pulse-block eigenfunctions](assets/screenshots/fig1_eigenfunctions.png) | ![Spectral ratios and their variation across modes](assets/screenshots/fig2_spectrum.png) |

<details>
<summary>🖼️ Open the complete seven-figure gallery</summary>

| Figure | Preview |
| :--- | :--- |
| 3 · Euler zero counts | ![Euler graininess and generalized-zero counts](assets/screenshots/fig3_euler.png) |
| 4 · Hybrid dynamics | ![Hybrid dynamics on dense blocks with dotted gap guides](assets/screenshots/fig4_hybrid.png) |
| 5 · Square-root scale | ![Square-root profiles and cumulative generalized-zero counts](assets/screenshots/fig5_sqrt.png) |
| 6 · Riccati sign loss | ![Riccati denominator sign loss on an alternating grid](assets/screenshots/fig6_riccati.png) |
| 7 · Disconjugacy | ![Pulse threshold and independent solutions sharing a gap](assets/screenshots/fig7_disconjugacy.png) |

</details>

<details>
<summary>🎬 Recording and screenshot placeholders</summary>

![Placeholder for a terminal GIF recorded from a verification run](assets/demo-placeholder.svg)

Add your recording as `assets/demos/verification.gif` and replace the placeholder
link. Additional screenshot slots are listed in [assets/README.md](assets/README.md).

</details>

Solid arcs belong to dense blocks. Dotted gap connectors are guides;
they assign no solution values inside missing intervals.

## 🧪 Verification

```bash
python -m pip install -r requirements-dev.txt
ruff check .
ruff format --check .
pytest -q
```

`pytest` runs analytical-limit and boundary-convention tests, then the full
535-assertion audit. For a shorter development loop, use `pytest -q -m "not slow"`.

CI checks Python 3.11, 3.12, and 3.13, builds an installable wheel,
and regenerates all seven figures. Its downloadable artifact contains fresh results.

The suite checks finite identities and numerical assertions. General statements
and infinite-horizon results require mathematical proofs; see [the conventions](docs/mathematics.md).

## 🗂️ Repository map

| Path | Purpose |
| :--- | :--- |
| [`src/sturmian/`](src/sturmian/) | Model, verification engine, figure generator, and CLI |
| [`tests/`](tests/) | Analytical, domain, CLI, and complete verification tests |
| [`data/`](data/) | Reference CSV tables |
| [`reports/`](reports/) | Re-run verification evidence and environment details |
| [`assets/`](assets/) | Local SVG branding, PNG figures, and demo placeholders |
| [`docs/`](docs/) | Mathematics, API, reproducibility, and publishing instructions |
| [`.github/`](.github/) | CI, labels, maintenance, releases, and community templates |

## 🤝 Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow.
Bug reports should include the time scale, gains, boundary conditions, and a
minimal example. Mathematical changes need an independent check.

[Code of conduct](CODE_OF_CONDUCT.md) · [Security policy](SECURITY.md) · [Changelog](CHANGELOG.md)

<img src="assets/divider.svg" alt="" width="100%" />

## 📈 Community and star history

<details>
<summary>📊 Maintainer activity and repository growth</summary>

The stats and contribution graph describe the maintainer's public GitHub account.
The star-history chart tracks this repository after it is published and receives stars.

<a href="https://github.com/Khuddush89">
  <img src="https://github-readme-stats.vercel.app/api?username=Khuddush89&amp;show_icons=true&amp;hide_border=true&amp;bg_color=0D1020&amp;title_color=8B5CF6&amp;icon_color=EC4899&amp;text_color=CBD5E1" alt="Public account statistics for Khuddush89" />
</a>

<a href="https://github.com/Khuddush89">
  <img src="https://github-readme-activity-graph.vercel.app/graph?username=Khuddush89&amp;bg_color=0d1020&amp;color=cbd5e1&amp;line=8b5cf6&amp;point=ec4899&amp;area=true&amp;hide_border=true" alt="Recent public contribution activity for Khuddush89" width="100%" />
</a>

<a href="https://www.star-history.com/#Khuddush89/Sturmian-Time-Scales&amp;Date">
  <img src="https://api.star-history.com/svg?repos=Khuddush89/Sturmian-Time-Scales&amp;type=Date" alt="Star history for Sturmian-Time-Scales" width="100%" />
</a>

Live cards use third-party services. Local branding and research figures remain
available if those services are temporarily unavailable.

</details>

## 📜 License and citation

Released under the [MIT License](LICENSE).
Use [CITATION.cff](CITATION.cff) to cite the software, data, or figures,
and include the commit or release used in your experiment.

<div align="center">

**Explore the mathematics. Reproduce the evidence.**

[Documentation](docs/README.md) · [Report a bug](https://github.com/Khuddush89/Sturmian-Time-Scales/issues/new?template=bug_report.md) · [Suggest a feature](https://github.com/Khuddush89/Sturmian-Time-Scales/issues/new?template=feature_request.md)

<img src="assets/footer-wave.svg" alt="Decorative purple, pink, and cyan footer wave" width="100%" />

</div>
