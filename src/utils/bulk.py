from pathlib import Path
from typing import Mapping

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write


ELEMENT_ORDER = ("Fe", "Ni", "Al")
FRACTION_TOL = 1.0e-8

def validate_composition(
    composition: Mapping[str, float],
    *,
    tolerance: float = FRACTION_TOL,
) -> dict[str, float]:
    """Validate and normalize a complete Fe-Ni-Al atomic-fraction mapping.

    All three elements must be present. Values must be finite and lie in [0, 1],
    and their sum must equal one within ``tolerance``.

    The returned mapping follows ``ELEMENT_ORDER``.
    """

    if set(composition) != set(ELEMENT_ORDER):
        missing = sorted(set(ELEMENT_ORDER) - set(composition))
        extra = sorted(set(composition) - set(ELEMENT_ORDER))

        details = []

        if missing:
            details.append(f"missing elements: {', '.join(missing)}")
        if extra:
            details.append(f"unsupported elements: {', '.join(extra)}")

        raise ValueError("Invalid composition mapping (" + "; ".join(details) + ").")

    fractions = {}

    for element in ELEMENT_ORDER:
        try:
            value = float(composition[element])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Fraction for {element} must be numeric.") from exc
        if not np.isfinite(value):
            raise ValueError(f"Fraction for {element} must be finite.")
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"Fraction for {element} must lie between 0 and 1.")
        fractions[element] = value

    total = sum(fractions.values())

    if not np.isclose(total, 1.0, rtol=0.0, atol=tolerance):
        raise ValueError(
            f"Atomic fractions must sum to 1; received {total:.12g}."
        )

    # Correct only the tiny floating-point discrepancy permitted above.
    return {element: fractions[element] / total for element in ELEMENT_ORDER}

def allocate_integer_counts(
    composition: Mapping[str, float],
    n_atoms: int,
) -> dict[str, int]:
    
    """Allocate integer species counts using the largest-remainder method.

    The counts sum exactly to ``n_atoms`` and minimize rounding discrepancies
    in the usual largest-remainder sense.
    """

    fractions = validate_composition(composition)

    if isinstance(n_atoms, bool) or not isinstance(n_atoms, (int, np.integer)):
        raise TypeError("n_atoms must be an integer.")
    if n_atoms < 1:
        raise ValueError("n_atoms must be positive.")

    raw = {element: fractions[element] * n_atoms for element in ELEMENT_ORDER}
    counts = {element: int(np.floor(raw[element])) for element in ELEMENT_ORDER}
    remainder = n_atoms - sum(counts.values())

    ranking = sorted(
        ELEMENT_ORDER,
        key=lambda element: raw[element] - counts[element],
        reverse=True,
    )
    
    for element in ranking[:remainder]:
        counts[element] += 1

    return counts

def generate_bulk(
    composition: Mapping[str, float],
    size: tuple[int, int, int] = (4, 4, 4),
    lattice_constant: float = 3.60,
    seed: int | None = None,
) -> Atoms:
    
    """Return a random chemical occupation of an FCC conventional supercell.

    Parameters
    ----------
    composition
        Complete atomic-fraction mapping, e.g.
        ``{"Fe": 0.25, "Ni": 0.50, "Al": 0.25}``.
    size
        Positive repeats ``(nx, ny, nz)`` of the conventional cubic FCC cell.
    lattice_constant
        Conventional cubic lattice parameter in Angstrom.
    seed
        Seed for NumPy's random generator. Set it for reproducibility.

    Returns
    -------
    ase.Atoms
        Periodic FCC supercell with randomized element occupations. The cell
        geometry is ideal and is not energy-minimized by this function.
    """

    fractions = validate_composition(composition)

    if len(size) != 3:
        raise ValueError("size must contain exactly three integers (nx, ny, nz).")
    
    if any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) for n in size):
        raise TypeError("Every size entry must be an integer.")
    
    if any(n < 1 for n in size):
        raise ValueError("All size entries must be positive.")
    
    if not np.isfinite(lattice_constant) or lattice_constant <= 0:
        raise ValueError("lattice_constant must be finite and positive.")

    # Use Fe only to establish the FCC parent geometry; all site occupations
    # are assigned below.
    atoms = bulk(
        "Fe",
        crystalstructure="fcc",
        a=float(lattice_constant),
        cubic=True,
    ).repeat(tuple(int(n) for n in size))

    counts = allocate_integer_counts(fractions, len(atoms))
    symbols = [
        element
        for element in ELEMENT_ORDER
        for _ in range(counts[element])
    ]

    rng = np.random.default_rng(seed)
    rng.shuffle(symbols)
    atoms.set_chemical_symbols(symbols)

    # Useful provenance for in-memory workflows. ASE's LAMMPS data writer
    # does not preserve arbitrary Python metadata, so persist metadata separately
    # if it is needed after writing.
    atoms.info["composition_target"] = fractions
    atoms.info["atom_counts"] = counts
    atoms.info["element_order"] = list(ELEMENT_ORDER)
    atoms.info["random_seed"] = seed
    atoms.info["lattice_constant_angstrom"] = float(lattice_constant)
    atoms.info["supercell_repeats"] = tuple(int(n) for n in size)

    return atoms

def write_bulk(atoms: Atoms, path: str | Path) -> Path:
    """Write an ASE structure as a LAMMPS atomic-style data file.

    Atom types are assigned by ``ELEMENT_ORDER``: Fe=1, Ni=2, Al=3.
    Parent directories are created as needed. Returns the output path.
    """
    path = Path(path)
    if path.suffix.lower() != ".data":
        raise ValueError("LAMMPS data output must use the .data extension.")

    symbols = set(atoms.get_chemical_symbols())
    unexpected = symbols - set(ELEMENT_ORDER)
    if unexpected:
        raise ValueError(f"Structure contains unsupported species: {sorted(unexpected)}")
    if not len(atoms):
        raise ValueError("Cannot write an empty structure.")

    path.parent.mkdir(parents=True, exist_ok=True)
    write(
        str(path),
        atoms,
        format="lammps-data",
        atom_style="atomic",
        specorder=list(ELEMENT_ORDER),
    )
    return path

def composition_report(atoms: Atoms) -> dict[str, dict[str, float | int]]:
    """Return atom counts and actual atomic fractions in canonical order."""
    
    if not len(atoms):
        raise ValueError("Cannot report composition of an empty structure.")

    symbols = atoms.get_chemical_symbols()
    result = {}
    for element in ELEMENT_ORDER:
        count = symbols.count(element)
        result[element] = {
            "count": count,
            "fraction": count / len(atoms),
        }
    return result