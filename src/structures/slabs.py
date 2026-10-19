from math import gcd
from functools import reduce
from pathlib import Path
from typing import Sequence
from dataclasses import dataclass

import numpy as np
from ase import Atoms
from ase.build import surface
from ase.io import write

from .bulk import ELEMENT_ORDER, composition_report, CompositionReport

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


def generate_slab(
    bulk_atoms: Atoms,
    miller_indices: Sequence[int],
    layers: int = 1,
    vacuum: float = 15.0,
) -> Atoms:
    """Construct an oriented slab from an existing periodic bulk structure.

    Parameters
    ----------
    bulk_atoms
        Parent periodic ASE structure, normally returned by ``generate_bulk``.
        Its chemical occupations are retained in the slab's repeated parent
        basis; this function does not randomize or reassign species.
    miller_indices
        Miller indices such as ``(1, 0, 0)``, ``(1, 1, 0)`` or ``(1, 1, 1)``.
    layers
        Number of repeats of the supplied parent cell in the surface-normal
        direction. Since the parent is a supercell, this is NOT generally the
        number of individual atomic planes. Inspect the resulting thickness.
    vacuum
        Vacuum thickness in Angstrom added along the surface-normal cell axis.

    Returns
    -------
    ase.Atoms
        Oriented slab with in-plane periodicity and non-periodic boundary
        condition along the surface normal. The original bulk object is not
        modified.

    Notes
    -----
    ASE's surface builder treats the supplied Atoms cell as the parent lattice.
    If that cell is already a large disordered supercell, the resulting slab
    may be substantially larger than a slab built from a primitive cell.
    Always inspect the atom count, geometry, layer composition, and termination.
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

    hkl = _normalize_miller_indices(miller_indices)

    symbols = set(bulk_atoms.get_chemical_symbols())
    unexpected = symbols - set(ELEMENT_ORDER)

    if unexpected:
        raise ValueError(f"Parent contains unsupported species: {sorted(unexpected)}")


    # ASE constructs the oriented surface cell from the parent periodic cell.
    # Keep the parent untouched; surface() creates a new Atoms object.
    slab = surface(
        bulk_atoms,
        indices=hkl,
        layers=int(layers),
        vacuum=None,
        periodic=False,
    )


    # Add vacuum along the surface-normal cell axis and set slab periodicity.
    slab.center(vacuum=float(vacuum), axis=2)
    slab.pbc = (True, True, False)


    # Record provenance for downstream scripts. This metadata remains in memory;
    # save a sidecar JSON if provenance must survive writing to LAMMPS data.
    slab.info["parent_n_atoms"] = len(bulk_atoms)
    slab.info["miller_indices"] = hkl
    slab.info["layers_parameter"] = int(layers)
    slab.info["vacuum_angstrom"] = float(vacuum)
    slab.info["element_order"] = list(ELEMENT_ORDER)
    if "random_seed" in bulk_atoms.info:
        slab.info["parent_random_seed"] = bulk_atoms.info["random_seed"]
    if "composition_target" in bulk_atoms.info:
        slab.info["parent_composition_target"] = bulk_atoms.info[
            "composition_target"
        ]

    return slab


def write_slab(slab: Atoms, path: str | Path) -> Path:
    """Write a slab to a LAMMPS atomic-style .data file.

    LAMMPS type IDs are fixed by ELEMENT_ORDER: Fe=1, Ni=2, Al=3.
    """
    path = Path(path)
    if path.suffix.lower() != ".data":
        raise ValueError("LAMMPS data output must use the .data extension.")
    if not len(slab):
        raise ValueError("Cannot write an empty slab.")

    symbols = set(slab.get_chemical_symbols())
    unexpected = symbols - set(ELEMENT_ORDER)
    if unexpected:
        raise ValueError(f"Slab contains unsupported species: {sorted(unexpected)}")

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