# 🚀 Publish on GitHub

## Repository metadata

| Field | Value |
| :--- | :--- |
| Name | `Sturmian-Time-Scales` |
| Owner used in links | `Khuddush89` |
| Description | Reproducible Sturmian spectra, gauge identities, and zero counts on time scales with Python. |
| Default branch | `main` |
| Software version | `0.1.0` |
| License | MIT |

Topics: `time-scales`, `sturmian-theory`, `conformable-derivative`, `sturm-liouville`, `oscillation-theory`, `spectral-theory`, `numerical-analysis`, `reproducible-research`, `scientific-computing`, `python`, `numpy`, `scipy`, `matplotlib`.

All remote links use the target namespace `Khuddush89/Sturmian-Time-Scales`.
If publishing under another name, update that string across the tree and the
account name in the README's profile cards before committing.

## Exact git init → push sequence

Extract the ZIP and run these commands inside its `Sturmian-Time-Scales` folder.
Install GitHub CLI and sign in with your GitHub account if needed.

```bash
cd Sturmian-Time-Scales
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
ruff check .
ruff format --check .
pytest -q

git init -b main
git add .
git commit -m "Initial release: reproducible Sturmian time-scale computations"
gh auth login
gh repo create Khuddush89/Sturmian-Time-Scales --public --source=. --remote=origin --description "Reproducible Sturmian spectra, gauge identities, and zero counts on time scales with Python."
git push -u origin main
gh repo edit Khuddush89/Sturmian-Time-Scales --add-topic time-scales --add-topic sturmian-theory --add-topic conformable-derivative --add-topic sturm-liouville --add-topic oscillation-theory --add-topic spectral-theory --add-topic numerical-analysis --add-topic reproducible-research --add-topic scientific-computing --add-topic python --add-topic numpy --add-topic scipy --add-topic matplotlib
```

If the empty repository already exists, replace `gh repo create ...` with:

```bash
git remote add origin https://github.com/Khuddush89/Sturmian-Time-Scales.git
git push -u origin main
```

Create an empty remote without GitHub-generated README or license files when
using this local tree. No remote repository has been created by preparing this ZIP.

## Labels

Run after creating the repository. These commands create or update the labels
used by issue templates, release notes, and maintenance workflows.

```bash
gh label create "bug" --repo Khuddush89/Sturmian-Time-Scales --color "EC4899" --description "Incorrect behavior or numerical result" --force
gh label create "enhancement" --repo Khuddush89/Sturmian-Time-Scales --color "8B5CF6" --description "New capability or improvement" --force
gh label create "documentation" --repo Khuddush89/Sturmian-Time-Scales --color "06B6D4" --description "Documentation and examples" --force
gh label create "mathematics" --repo Khuddush89/Sturmian-Time-Scales --color "8B5CF6" --description "Scientific identities or correctness" --force
gh label create "tests" --repo Khuddush89/Sturmian-Time-Scales --color "22C55E" --description "Regression tests and verification" --force
gh label create "visuals" --repo Khuddush89/Sturmian-Time-Scales --color "EC4899" --description "Figures, branding, and visual assets" --force
gh label create "data" --repo Khuddush89/Sturmian-Time-Scales --color "06B6D4" --description "Reference tables and verification evidence" --force
gh label create "ci" --repo Khuddush89/Sturmian-Time-Scales --color "6366F1" --description "GitHub Actions and automation" --force
gh label create "dependencies" --repo Khuddush89/Sturmian-Time-Scales --color "0EA5E9" --description "Dependency updates" --force
gh label create "security" --repo Khuddush89/Sturmian-Time-Scales --color "DC2626" --description "Security-sensitive work" --force
gh label create "inactive" --repo Khuddush89/Sturmian-Time-Scales --color "94A3B8" --description "No recent activity; remains open" --force
gh label create "pinned" --repo Khuddush89/Sturmian-Time-Scales --color "F59E0B" --description "Keep active regardless of age" --force
gh label create "breaking" --repo Khuddush89/Sturmian-Time-Scales --color "DC2626" --description "Breaking interface or behavior change" --force
gh label create "help wanted" --repo Khuddush89/Sturmian-Time-Scales --color "14B8A6" --description "Contributions welcome" --force
gh label create "good first issue" --repo Khuddush89/Sturmian-Time-Scales --color "A78BFA" --description "Suitable for a first contribution" --force
```

## Suggested GitHub settings

| Setting | Recommendation |
| :--- | :--- |
| About | Use the description and topics above |
| Features | Enable Issues and Discussions; disable unused Wiki and Projects |
| Default branch | `main` |
| Branch rule | Require pull requests, resolved conversations, and passing CI |
| Required checks | `Lint`, `Tests (Python 3.11)`, `Tests (Python 3.12)`, `Tests (Python 3.13)`, `Reproduce figures` |
| Review policy | Require one approval when an independent reviewer is available |
| History | Block force pushes and branch deletion; use squash merging and delete merged branches |
| Actions token | Keep repository defaults read-only; workflows request their own scoped permissions |
| Dependencies | Enable dependency graph, Dependabot alerts, and security updates |
| Security | Enable private vulnerability reporting and available secret scanning |
| Social preview | Upload `assets/social-preview.png` (1280 × 640) |
| Funding | Configure a verified account before activating `FUNDING.yml` |

Wait for the first successful CI run before selecting required status checks.
Pinned action commits are refreshed by Dependabot. Labeling reads base-branch
configuration only; it never executes pull-request code.

## First release

Release Drafter creates drafts only. Review its version, changelog, and CI
artifacts before manually publishing `v0.1.0`. Update `pyproject.toml`,
`src/sturmian/__init__.py`, `CITATION.cff`, and the README version badge together.

The stale workflow marks inactive items but never closes them automatically.
Scientific correctness and security labels are exempt from stale marking.

## Official references

- [Branch protection](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/about-protected-branches)
- [Dependabot options](https://docs.github.com/en/code-security/reference/supply-chain-security/dependabot-options-reference)
- [Pull request labeler](https://github.com/actions/labeler)
- [Stale action](https://github.com/actions/stale)
- [Release Drafter](https://github.com/release-drafter/release-drafter)
- [readme-typing-svg](https://github.com/DenverCoder1/readme-typing-svg)
- [Skill Icons](https://github.com/tandpfun/skill-icons)
- [GitHub Readme Stats](https://github.com/anuraghazra/github-readme-stats)
- [Activity graph](https://github.com/Ashutosh00710/github-readme-activity-graph)
- [Star History](https://github.com/star-history/star-history)
