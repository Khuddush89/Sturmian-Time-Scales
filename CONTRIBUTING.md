# 🤝 Contributing

Contributions to the computations, documentation, and reproducibility are welcome.
For a new numerical method or mathematical claim, open an issue before a large change.

## 🛠️ Set up

```bash
git clone https://github.com/Khuddush89/Sturmian-Time-Scales.git
cd Sturmian-Time-Scales
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
git switch -c fix/describe-your-change
```

## 🧪 Validate

```bash
ruff check .
ruff format --check .
pytest -q
sturmian verify --output outputs/contribution
sturmian figures --output outputs/contribution
```

For quick edits, use `pytest -q -m "not slow"` first.
Run the complete suite for changes to the model, verification, or figures.

## 📐 Scientific changes

- State the time scale and distinguish dense blocks from gaps.
- Preserve the generalized-zero convention in [mathematics.md](docs/mathematics.md).
- Explain any new tolerance and use an analytical or independent reference.
- Do not loosen a tolerance just to make a failing check pass.
- Keep seeds explicit. Use positive rescaling for sign counts, not amplitudes.
- Update reference data and PNGs only when the scientific change requires it.

## 📦 Pull requests

Keep each pull request focused. Describe the trigger, resulting behavior, and
validation. Include before/after figures when a plot changes and add an entry
under `Unreleased` in [CHANGELOG.md](CHANGELOG.md).

Do not commit `outputs/`, environments, caches, typesetting files, or credentials.
Contributions are licensed under the project's MIT License.
Follow the [code of conduct](CODE_OF_CONDUCT.md).
