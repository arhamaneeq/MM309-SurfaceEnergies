# Fe–Ni–Al surface-energy report draft (through convergence)

This archive contains the report draft using the requested eight-section outline. Sections 1–5 are drafted. Sections 6–8 remain placeholders because production results, interpretation, and final conclusions are not yet available.

**Page limit:** 7 pages total, including title, abstract, tables, figures, and references. The compact two-column layout is intended to resemble the supplied sample lab report.

## Compile on Windows
Open PowerShell in this directory and run:

```powershell
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Important scientific notes
- The production scripts define three compositions: equiatomic FeNiAl (33.33 at.% each), Fe20Ni50Al30, and Fe50Ni30Al20.
- The saved convergence test is a separate Fe25Ni25Al50 FCC (111) model. Do not treat it as a production composition or use its bulk reference for production surface energies.
- A mapping discrepancy must be resolved: convergence documentation says type 1 = Fe, type 2 = Ni, type 3 = Al; production structure scripts document type 1 = Al, type 2 = Fe, type 3 = Ni. Check the LAMMPS data-file Masses section and pair_coeff element order before trusting energies.
- In the saved test, the 12- and 14-layer surface-energy estimates differ by about 0.0052 J/m² (about 0.36%). The result is preliminary and does not prove convergence for the production (110) surface or other compositions.
- Vacuum padding is per side. The 2 Å case differs from the 25 Å reference by about 0.75037 eV; 4 Å and larger match at the reported precision for this test. A 10 Å per-side gap is a conservative provisional choice, not a universal threshold.
- The Sharifi–Wick paper is cited as the source publication associated with the supplied MEAM files. The model's stated limitations for temperature-dependent properties are noted.
- The production plan contains 3 compositions × 3 orientations × 3 temperatures × 3 configurations = 81 surface calculations. No production results, configuration uncertainty, or final conclusions are claimed yet.

## Repository hygiene
Keep simulation inputs, potential files, logs, and machine-readable CSV files tracked in the project repository. Do not commit generated PDF/build artifacts unless your team specifically wants them versioned.
