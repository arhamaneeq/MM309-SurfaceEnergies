# Fe-Ni-Al Surface Slab Convergence Tests

## 1. Purpose

This directory contains initial convergence tests for the Fe-Ni-Al alloy
surface-energy project.

The tests investigate:

1. Slab thickness for the FCC (111) surface.
2. Separation between periodically repeated slabs along the surface normal.

Only the selected controlled slab-thickness dataset is retained for analysis.

## 2. Alloy composition and potential

Target composition:

| Element | Atomic percentage | LAMMPS atom type |
|---|---:|---:|
| Fe | 25% | 1 |
| Ni | 25% | 2 |
| Al | 50% | 3 |

The calculations use the MEAM potential files:

- `potential/library.meam`
- `potential/FeNiAl.meam`

The potential source and publication are documented below. Its accuracy
for the specific surface-energy calculations and temperatures in this
assignment still requires critical assessment.

### Potential provenance and limitations

The local `FeNiAl.meam` and `library.meam` parameter values were checked
against the Fe-Ni-Al MEAM files listed in the NIST Interatomic Potentials
Repository.

Reference:

H. Sharifi and C. D. Wick, "Developing interatomic potentials for complex
concentrated alloys of Cu, Ti, Ni, Cr, Co, Al, Fe, and Mn,"
Computational Materials Science 248 (2025), 113595.
DOI: 10.1016/j.commatsci.2024.113595.

Official record:
https://www.ctcms.nist.gov/potentials/system/Al-Fe-Ni/

The NIST entry notes that this potential was not optimized for
temperature-dependent properties and was not fitted to density, thermal
expansion, or thermal conductivity data. Its performance for the specific
surface-energy calculations in this project must therefore be assessed rather
than assumed. This limitation is particularly relevant when interpreting
temperature-dependent results.

## 3. Bulk reference

The controlled slab surface energies use a bulk reference energy of
approximately

`E_bulk = -4.174619518535 eV/atom`

This value comes from the tight bulk relaxation. The initial bulk input is
`FeNiAl_bulk.data`; the corresponding LAMMPS input is
`bulk/in.bulk_relax_tight`.

The bulk input structure is retained so that the calculation can be repeated.
The relaxed bulk structure is a generated output and is not required to be
committed when it can be regenerated using the input and LAMMPS script.

## 4. Slab-thickness convergence

The FCC (111) slab series contains 4, 6, 8, 10, 12, and 14 atomic planes.
The in-plane dimensions are held constant, with a surface area of approximately
96.711966 Å² per surface.

Each slab has two surfaces. The surface energy is calculated as

`gamma = (E_slab - N * E_bulk) / (2 * A)`

where:

- `E_slab` is the relaxed total slab energy in eV.
- `N` is the number of atoms in the slab.
- `E_bulk` is the bulk reference energy in eV/atom.
- `A` is the area of one surface in Å².
- The factor of 2 accounts for the two slab surfaces.

The conversion used is

`1 eV/Å² = 16.02176634 J/m²`.

### Controlled results

| Layers | Atoms | Relaxed energy (eV) | Surface energy (J/m²) |
|---:|---:|---:|---:|
| 4 | 64 | -252.043494 | 1.2534 |
| 6 | 96 | -386.264657 | 1.2010 |
| 8 | 128 | -517.120893 | 1.4272 |
| 10 | 160 | -650.773999 | 1.4218 |
| 12 | 192 | -784.234158 | 1.4324 |
| 14 | 224 | -917.759036 | 1.4376 |

Machine-readable results are in `slab/results/111_thickness_convergence.csv`.
The LAMMPS inputs are in `slab/inputs/`.

The 12- and 14-layer surface energies differ by approximately
0.0052 J/m² (about 0.36% relative to the 14-layer result). However, the
change between 6 and 8 layers is substantial. The complete series should
therefore be considered when selecting a production slab thickness.

These results apply to the tested alloy arrangement, surface termination,
potential, and relaxation settings. Additional configurations may be needed
to assess sensitivity to chemical disorder.

## 5. Vacuum-separation convergence

The vacuum test uses the relaxed controlled 10-layer FCC (111) slab. The
listed vacuum value is the padding on **each side** of the slab. The vacuum
input files use periodic boundary conditions along all three directions.

| Vacuum per side (Å) | Cell height (Å) | Final energy (eV) | Energy difference from 25 Å case (eV) |
|---:|---:|---:|---:|
| 2 | 23.952639 | -651.524369 | 0.750370 |
| 4 | 27.952639 | -650.773999 | 0.000000 |
| 6 | 31.952639 | -650.773999 | 0.000000 |
| 8 | 35.952639 | -650.773999 | 0.000000 |
| 10 | 39.952639 | -650.773999 | 0.000000 |
| 15 | 49.952639 | -650.773999 | 0.000000 |
| 20 | 59.952639 | -650.773999 | 0.000000 |
| 25 | 69.952639 | -650.773999 | 0.000000 |

The 4 Å-per-side case matches the 25 Å reference to the reported energy
precision, while the 2 Å-per-side case differs by approximately 0.75037 eV.
This conclusion is specific to this slab, potential, and numerical setup; it
should not be treated as a universal vacuum requirement.

The data are in `vacuum_convergence/results/vacuum_convergence.csv`.
The scripts used to generate the vacuum test structures are in
`vacuum_convergence/`.

The vacuum generators depend on the relaxed controlled 10-layer slab produced
by the corresponding slab LAMMPS input. The starting slab input structure and
the scripts needed to generate the vacuum cases should be retained.

## 6. Reproduction files

- `FeNiAl_bulk.data`: starting bulk structure.
- `bulk/in.bulk_relax`: standard bulk relaxation input.
- `bulk/in.bulk_relax_tight`: tighter bulk relaxation input used for the
  reference energy above.
- `bulk/in.meam_test`: potential smoke-test input.
- `potential/`: MEAM potential files.
- `slab/structures/FeNiAl_111_L*_2x2_controlled.data`: starting structures
  for the controlled thickness series.
- `slab/inputs/`: LAMMPS inputs for the controlled thickness series.
- `slab/results/111_thickness_convergence.csv`: summarized thickness results.
- `vacuum_convergence/`: vacuum-generation scripts, LAMMPS inputs, and results.

Run the LAMMPS commands from the repository root because the input files use
repository-relative paths.

Generated relaxed structures and simulation logs are outputs and are excluded
from version control. The essential starting structures, input scripts,
generation scripts, potential files, and summarized results should be retained
so the work can be reproduced.

## 7. Limitations and next steps

- These are initial convergence tests, not the complete surface-energy study.
- Thickness results may depend on the particular alloy arrangement and surface
  termination.
- The vacuum test assesses total-energy behavior for one representative slab.
- The potential source, parameter provenance, and validation must be verified
  and documented.
- The selected model should be justified against an explicit convergence
  tolerance before production calculations.
- Further assignment work includes the required surfaces, temperatures,
  configurations, uncertainties, and comparison of the results.
