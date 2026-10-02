# 📦 Provenance and scope

The input was the uploaded `Sturmian_Reviewed.zip`. Its numerical supplement
contained three Python modules, 20 CSVs, seven PNG figures, and verification
evidence. The manuscript and typesetting outputs are kept outside this repository.

The numerical material accompanies *Gauge Reduction and Sturmian Consequences
for Proportional Conformable Dynamic Equations on Time Scales*.
No publication venue, DOI, or acceptance status is asserted here.

## Source mapping

| Supplied path | Repository path |
| :--- | :--- |
| `code/sturmian_model.py` | `src/sturmian/model.py` |
| `code/verify.py` | `src/sturmian/verification.py` |
| `code/figures.py` | `src/sturmian/figures.py` |
| `figures/*.png` | `assets/screenshots/*.png` |
| `data/*.csv` | `data/*.csv`, regenerated with the packaged code |
| Verification report and log | Fresh evidence in `reports/` |

Packaging changes replace wildcard imports, add explicit output paths, format
the code, record the current environment, and expose CLI commands. PNG generation
retains the original plotting conventions.

Correctness changes are limited to sign comparisons that avoid multiplication
underflow and clearer validation of grids, gains, spectral mass, and gap
admissibility. The original 535-assertion audit remains intact.

Apple resource forks, operating-system metadata, LaTeX files, PDFs, temporary
build files, and manuscript editorial diffs are excluded from the ready-to-push tree.
The input archive remains the source for those separate materials.

The supplied claim of 535 explicit checks is re-run and recorded. Earlier
unavailable claims about other verification counts are not used as repository badges.
