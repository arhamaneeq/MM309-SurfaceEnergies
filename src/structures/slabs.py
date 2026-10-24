
from math import gcd
from functools import reduce
from pathlib import Path
from typing import Sequence
from dataclasses import dataclass

import numpy as np
from ase import Atoms
from ase.build import bulk as ase_bulk, surface
from ase.io import write

from .bulk import (
    ELEMENT_ORDER,
    CompositionReport,
    allocate_integer_counts,
    composition_report,
    validate_composition,
)


@dataclass(frozen=True)
class SlabReport:
    n_atoms: int
    composition: CompositionReport
    surface_area_angstrom2: float
    atomic_thickness_angstrom: float
    cell_height_normal_angstrom: float
    cell_lengths_angstrom: tuple[float, float, float]
    pbc: tuple[bool, bool, bool]
    miller_indices: tuple[int, int, int] | None
    parent_random_seed: int | None


def _normalize_miller_indices(
    miller_indices: Sequence[int],
) -> tuple[int, int, int]:
    """Validate Miller indices and reduce them to their smallest integer form."""
    if len(miller_indices) != 3:
        raise ValueError("miller_indices must contain exactly three integers.")

    values = []
    for value in miller_indices:
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError("Miller indices must be integers.")
        values.append(int(value))

    if values == [0, 0, 0]:
        raise ValueError("Miller indices cannot all be zero.")

    divisor = reduce(gcd, (abs(value) for value in values))
    return tuple(value // divisor for value in values)


def _get_parent_composition(bulk_atoms: Atoms) -> dict[str, float]:
    """Obtain target composition from bulk metadata, falling back to actual composition."""
    target = bulk_atoms.info.get("composition_target")

    if target is None:
        report = composition_report(bulk_atoms)
        target = {
            element: getattr(report, element).fraction
            for element in ELEMENT_ORDER
        }

    return validate_composition(target)


def _get_lattice_constant(
    bulk_atoms: Atoms,
    lattice_constant: float | None,
) -> float:
    """Get the FCC conventional lattice parameter in Angstrom."""
    if lattice_constant is None:
        lattice_constant = bulk_atoms.info.get("lattice_constant_angstrom")

    if lattice_constant is None:
        raise ValueError(
            "Could not determine lattice constant. Pass lattice_constant "
            "explicitly or provide it in bulk_atoms.info."
        )

    lattice_constant = float(lattice_constant)

    if not np.isfinite(lattice_constant) or lattice_constant <= 0:
        raise ValueError("lattice_constant must be finite and positive.")

    return lattice_constant


def generate_slab(
    bulk_atoms: Atoms,
    miller_indices: Sequence[int],
    layers: int = 1,
    vacuum: float = 15.0,
    lattice_constant: float | None = None,
    inplane_repeats: tuple[int, int] = (4, 4),
) -> Atoms:
    """Construct an FCC slab, then randomly assign its chemical occupations.

    Parameters
    ----------
    bulk_atoms
        Existing periodic structure used as the source of target composition,
        lattice parameter, and random-seed metadata. Its atomic arrangement
        is NOT copied into the slab.
    miller_indices
        Conventional cubic Miller indices, e.g. (1, 0, 0), (1, 1, 0),
        or (1, 1, 1).
    layers
        Number of atomic layers requested from ASE's surface builder.
    vacuum
        Vacuum thickness in Angstrom on each side of the slab.
    lattice_constant
        Conventional cubic FCC lattice parameter in Angstrom. If omitted,
        it is read from bulk_atoms.info.
    inplane_repeats
        Repeats of the surface cell along its two in-plane directions.

    Returns
    -------
    ase.Atoms
        Randomly decorated FCC slab with periodic in-plane boundaries and
        a non-periodic surface-normal direction.

    Notes
    -----
    Element counts are allocated using the largest-remainder method, so
    the requested fractions are matched as closely as integer counts allow.
    Chemical occupations are shuffled deterministically when a parent seed
    is available.
    """

    if not isinstance(bulk_atoms, Atoms):
        raise TypeError("bulk_atoms must be an ase.Atoms object.")
    if len(bulk_atoms) == 0:
        raise ValueError("bulk_atoms cannot be empty.")
    if not np.all(bulk_atoms.pbc):
        raise ValueError("bulk_atoms must be periodic in all three directions.")

    if isinstance(layers, bool) or not isinstance(layers, (int, np.integer)):
        raise TypeError("layers must be an integer.")
    if layers < 1:
        raise ValueError("layers must be >= 1.")

    if not np.isfinite(vacuum) or vacuum <= 0:
        raise ValueError("vacuum must be finite and positive.")

    if (
        len(inplane_repeats) != 2
        or any(
            isinstance(n, bool) or not isinstance(n, (int, np.integer))
            for n in inplane_repeats
        )
        or any(n < 1 for n in inplane_repeats)
    ):
        raise ValueError("inplane_repeats must contain two positive integers.")

    hkl = _normalize_miller_indices(miller_indices)
    composition = _get_parent_composition(bulk_atoms)
    a = _get_lattice_constant(bulk_atoms, lattice_constant)

    unexpected = set(bulk_atoms.get_chemical_symbols()) - set(ELEMENT_ORDER)
    if unexpected:
        raise ValueError(
            f"Parent contains unsupported species: {sorted(unexpected)}"
        )

    # Build a pure conventional cubic FCC cell. The conventional cubic basis
    # ensures that the requested Miller indices refer to cubic crystallographic
    # axes, not to the non-orthogonal primitive FCC basis.
    fcc = ase_bulk(
        "Fe",
        crystalstructure="fcc",
        a=a,
        cubic=True,
    )

    slab = surface(
        fcc,
        indices=hkl,
        layers=int(layers),
        vacuum=None,
        periodic=False,
    )

    # Expand the surface area without increasing the slab thickness.
    slab = slab.repeat(
        (int(inplane_repeats[0]), int(inplane_repeats[1]), 1)
    )

    # Center the slab and add vacuum on both sides.
    slab.center(vacuum=float(vacuum), axis=2)
    slab.pbc = (True, True, False)

    # Allocate exact integer species counts for this slab's atom count.
    counts = allocate_integer_counts(composition, len(slab))
    symbols = [
        element
        for element in ELEMENT_ORDER
        for _ in range(counts[element])
    ]

    parent_seed = bulk_atoms.info.get("random_seed")

    # Derive a reproducible, orientation-specific decoration seed.
    decoration_seed = None
    if parent_seed is not None:
        decoration_seed = (
            int(parent_seed)
            + 10000 * hkl[0]
            + 100 * hkl[1]
            + hkl[2]
        )

    rng = np.random.default_rng(decoration_seed)
    rng.shuffle(symbols)
    slab.set_chemical_symbols(symbols)

    # Preserve useful provenance in memory.
    slab.info["parent_n_atoms"] = len(bulk_atoms)
    slab.info["miller_indices"] = hkl
    slab.info["layers_parameter"] = int(layers)
    slab.info["vacuum_angstrom"] = float(vacuum)
    slab.info["inplane_repeats"] = tuple(map(int, inplane_repeats))
    slab.info["lattice_constant_angstrom"] = a
    slab.info["element_order"] = list(ELEMENT_ORDER)
    slab.info["composition_target"] = composition
    slab.info["atom_counts"] = counts
    slab.info["random_seed"] = decoration_seed
    slab.info["parent_random_seed"] = parent_seed

    return slab


def write_slab(slab: Atoms, path: str | Path) -> Path:
    """Write a slab to a LAMMPS atomic-style .data file.

    LAMMPS type IDs follow ELEMENT_ORDER: Fe=1, Ni=2, Al=3.
    """
    path = Path(path)

    if path.suffix.lower() != ".data":
        raise ValueError("LAMMPS data output must use the .data extension.")
    if not len(slab):
        raise ValueError("Cannot write an empty slab.")

    unexpected = set(slab.get_chemical_symbols()) - set(ELEMENT_ORDER)
    if unexpected:
        raise ValueError(
            f"Slab contains unsupported species: {sorted(unexpected)}"
        )

    path.parent.mkdir(parents=True, exist_ok=True)

    write(
        str(path),
        slab,
        format="lammps-data",
        atom_style="atomic",
        specorder=list(ELEMENT_ORDER),
    )

    return path


def slab_report(slab: Atoms) -> SlabReport:
    """Summarize atom counts, composition, surface area, and geometry."""
    if not len(slab):
        raise ValueError("Cannot report an empty slab.")

    a_vec, b_vec, c_vec = np.asarray(slab.cell)
    normal_vector = np.cross(a_vec, b_vec)
    area = float(np.linalg.norm(normal_vector))

    if area <= 0:
        raise ValueError("Slab has a degenerate in-plane cell.")

    normal = normal_vector / area
    projected_positions = np.asarray(slab.positions) @ normal

    thickness = (
        float(np.ptp(projected_positions))
        if len(slab) > 1
        else 0.0
    )

    cell_height = abs(float(np.dot(c_vec, normal)))

    return SlabReport(
        n_atoms=len(slab),
        composition=composition_report(slab),
        surface_area_angstrom2=area,
        atomic_thickness_angstrom=thickness,
        cell_height_normal_angstrom=cell_height,
        cell_lengths_angstrom=tuple(float(x) for x in slab.cell.lengths()),
        pbc=tuple(bool(x) for x in slab.pbc),
        miller_indices=slab.info.get("miller_indices"),
        parent_random_seed=slab.info.get("parent_random_seed"),
    )
