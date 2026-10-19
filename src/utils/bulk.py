import argparse
import sys
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import write


ELEMENTS = ("Fe", "Ni", "Al")
FRACTION_TOL = 1.0e-8

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a random substitutional FCC Fe-Ni-Al alloy and write "
            "a LAMMPS-readable data file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--Fe", type=float, help="Fe atomic fraction")
    parser.add_argument("--Ni", type=float, help="Ni atomic fraction")
    parser.add_argument("--Al", type=float, help="Al atomic fraction")
    parser.add_argument(
        "--size",
        type=int,
        nargs="+",
        default=[4, 4, 4],
        metavar="N",
        help="FCC conventional-cell repeats: one integer (cubic) or nx ny nz",
    )
    parser.add_argument(
        "--lattice-constant",
        type=float,
        default=3.60,
        help="FCC lattice parameter in Angstrom; choose/justify for your model",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("fcc_random_alloy.data"),
        help="Output LAMMPS data filename",
    )
    return parser.parse_args()

def resolve_fractions(args):
    supplied = {
        "Fe": args.Fe,
        "Ni": args.Ni,
        "Al": args.Al,
    }
    provided = [element for element, value in supplied.items() if value is not None]

    if len(provided) < 2:
        raise ValueError("Specify at least two of --Fe, --Ni and --Al.")

    for element in provided:
        value = supplied[element]
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"--{element} must be between 0 and 1.")

    if len(provided) == 2:
        missing = next(element for element in ELEMENTS if supplied[element] is None)
        inferred = 1.0 - sum(supplied[element] for element in provided)
        if inferred < -FRACTION_TOL:
            raise ValueError(
                f"The supplied fractions sum to more than 1; inferred {missing} "
                f"fraction would be {inferred:.8g}."
            )
        supplied[missing] = max(0.0, inferred)
        print(f"Inferred --{missing} = {supplied[missing]:.8g}")
    else:
        total = sum(supplied[element] for element in ELEMENTS)
        if not np.isclose(total, 1.0, rtol=0.0, atol=FRACTION_TOL):
            raise ValueError(
                f"When all three fractions are supplied, they must sum to 1 "
                f"(received {total:.12g})."
            )

    # Normalize tiny floating-point sum discrepancies after validation.
    total = sum(supplied[element] for element in ELEMENTS)
    return {element: supplied[element] / total for element in ELEMENTS}

def allocate_integer_counts(fractions, n_atoms):
    """Largest-remainder allocation: integer counts sum exactly to n_atoms."""
    raw = {element: fractions[element] * n_atoms for element in ELEMENTS}
    counts = {element: int(np.floor(raw[element])) for element in ELEMENTS}
    remainder = n_atoms - sum(counts.values())

    order = sorted(
        ELEMENTS,
        key=lambda element: raw[element] - counts[element],
        reverse=True,
    )
    for element in order[:remainder]:
        counts[element] += 1
    return counts

def main():
    args = parse_args()

    try:
        fractions = resolve_fractions(args)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    if len(args.size) == 1:
        nx = ny = nz = args.size[0]
    elif len(args.size) == 3:
        nx, ny, nz = args.size
    else:
        print("Error: --size takes either one integer or three integers.", file=sys.stderr)
        return 2

    if min(nx, ny, nz) < 1:
        print("Error: all --size values must be positive integers.", file=sys.stderr)
        return 2
    if args.lattice_constant <= 0:
        print("Error: --lattice-constant must be positive.", file=sys.stderr)
        return 2

    # A conventional FCC cell contains four sites.
    atoms = bulk("Fe", crystalstructure="fcc",
                 a=args.lattice_constant, cubic=True)
    atoms = atoms.repeat((nx, ny, nz))
    n_atoms = len(atoms)

    counts = allocate_integer_counts(fractions, n_atoms)

    # Create the exact species list, then randomly permute site occupations.
    symbols = []
    for element in ELEMENTS:
        symbols.extend([element] * counts[element])

    rng = np.random.default_rng(args.seed)
    symbols = np.asarray(symbols, dtype=object)
    rng.shuffle(symbols)
    atoms.set_chemical_symbols(symbols.tolist())

    # specorder fixes LAMMPS atom-type IDs: 1=Fe, 2=Ni, 3=Al.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write(
        str(args.output),
        atoms,
        format="lammps-data",
        atom_style="atomic",
        specorder=list(ELEMENTS),
    )

    print("\nGenerated random FCC alloy")
    print("--------------------------")
    print(f"Output file       : {args.output}")
    print(f"Atom count        : {n_atoms}")
    print(f"Supercell repeats : {nx} x {ny} x {nz}")
    print(f"Lattice parameter : {args.lattice_constant:.6f} Angstrom")
    print(f"Cell lengths      : {nx * args.lattice_constant:.6f}, "
          f"{ny * args.lattice_constant:.6f}, "
          f"{nz * args.lattice_constant:.6f} Angstrom")
    print(f"Random seed       : {args.seed}")
    print("LAMMPS atom types : 1=Fe, 2=Ni, 3=Al")
    print("\nElement  Target fraction  Atom count  Actual fraction")
    for element in ELEMENTS:
        actual = counts[element] / n_atoms
        print(f"{element:>6}  {fractions[element]:>15.8f}  "
              f"{counts[element]:>10d}  {actual:>14.8f}")
    print(
        "\nNote: fractions are approximated by integer atom counts. "
        "The lattice parameter is an input assumption, not an equilibrium "
        "lattice parameter calculated by this script."
    )

if __name__ == "__main__":
    main()