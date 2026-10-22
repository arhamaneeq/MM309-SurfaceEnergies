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

STRUCTURES.mkdir(parents=True, exist_ok=True)

if not SOURCE.is_file():
    raise FileNotFoundError(f"Source slab not found: {SOURCE}")

for vacuum in (10, 15, 20, 25):
    slab = read(SOURCE, format="lammps-data", atom_style="atomic")

    # Add the specified vacuum on both sides of the slab.
    slab.center(vacuum=vacuum, axis=2)
    slab.pbc = (True, True, True)
    slab.wrap()

    output = STRUCTURES / f"FeNiAl_111_L10_vacuum{vacuum}.data"

    write(
        output,
        slab,
        format="lammps-data",
        atom_style="atomic",
        specorder=["Fe", "Ni", "Al"],
        masses=True,
        units="metal",
    )

    print(
        f"Vacuum: {vacuum:2d} Å per side | "
        f"Cell height: {slab.cell.lengths()[2]:.6f} Å | "
        f"Atoms: {len(slab)} | Output: {output.name}"
    )

print("Large-vacuum structure generation complete.")
