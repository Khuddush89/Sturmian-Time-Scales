# 📦 Changelog

## [Unreleased]

Add changes here as pull requests are merged.

## [0.1.0] — 2026-10-02

### Added

- Installable `src/sturmian` package and `sturmian verify` / `sturmian figures` commands.
- Analytical-limit, boundary-convention, CLI, and full verification tests.
- Gradient SVG branding, seven reference figures, and reproducibility documentation.
- CI, pull-request labeling, stale marking, Dependabot, and release drafting.
- Community templates, MIT License, and software citation metadata.

### Fixed

- Missing dependency instructions and obsolete filenames in the supplied README.
- Output paths tied to source location; runs now use an explicit output directory.
- Strict sign crossings hidden by floating-point underflow in `propagate`.
- Silent invalid grids and nonpositive spectral mass coefficients.
- Hard-coded verification dates; fresh reports now record the run environment.

This version identifies the prepared software snapshot. A GitHub release is
created only when a maintainer publishes it.
