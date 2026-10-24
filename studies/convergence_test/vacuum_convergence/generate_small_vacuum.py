from pathlib import Path
from ase.io import read, write

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
TEST = ROOT / "convergence_test"

SOURCE = (
    TEST / "slab" / "structures"
    / "FeNiAl_111_L10_2x2_controlled_relaxed.data"
)
STRUCTURES = BASE / "structures"
INPUTS = BASE / "inputs"

STRUCTURES.mkdir(parents=True, exist_ok=True)
INPUTS.mkdir(parents=True, exist_ok=True)

if not SOURCE.is_file():
    raise FileNotFoundError(f"Source slab not found: {SOURCE}")

for vacuum in (2, 4, 6, 8):
    slab = read(SOURCE, format="lammps-data", atom_style="atomic")

    # Add the specified vacuum on both sides of the slab.
    slab.center(vacuum=vacuum, axis=2)
    slab.pbc = (True, True, True)
    slab.wrap()

    datafile = STRUCTURES / f"FeNiAl_111_L10_vacuum{vacuum}.data"

    write(
        datafile,
        slab,
        format="lammps-data",
        atom_style="atomic",
        specorder=["Fe", "Ni", "Al"],
        masses=True,
        units="metal",
    )

    infile = INPUTS / f"in.vacuum{vacuum}"
    infile.write_text(
        f"""units metal
atom_style atomic
boundary p p p

read_data convergence_test/vacuum_convergence/structures/FeNiAl_111_L10_vacuum{vacuum}.data

pair_style meam
pair_coeff * * convergence_test/potential/library.meam Fe Ni Al convergence_test/potential/FeNiAl.meam Fe Ni Al

neighbor 2.0 bin
neigh_modify delay 10 check yes

thermo 100
thermo_style custom step atoms pe press pxx pyy pzz lx ly lz

min_style cg
minimize 0.0 1.0e-6 50000 500000

run 0

write_data convergence_test/vacuum_convergence/structures/FeNiAl_111_L10_vacuum{vacuum}_relaxed.data
"""
    )

    print(
        f"Vacuum: {vacuum:2d} Å per side | "
        f"Cell height: {slab.cell.lengths()[2]:.6f} Å | "
        f"Atoms: {len(slab)} | "
        f"Structure: {datafile.name} | Input: {infile.name}"
    )

print("Small-vacuum structure and input generation complete.")
