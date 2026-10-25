import numpy as np
from ase.build import bulk as ase_bulk, surface

LATTICE_CONSTANT = 3.60  # Illustrative only; replace with your actual a(T)
LAYERS = 12
INPLANE_REPEATS = (4, 4)
MILLER_INDICES = ((1, 1, 1), (1, 1, 0), (1, 0, 0))


def count_surface_atoms(slab, tolerance=1e-5):
    """Count atoms in each outermost atomic plane."""
    a_vec, b_vec, _ = np.asarray(slab.cell)
    normal = np.cross(a_vec, b_vec)
    normal /= np.linalg.norm(normal)

    heights = np.asarray(slab.positions) @ normal
    bottom = heights.min()
    top = heights.max()

    n_bottom = np.count_nonzero(
        np.isclose(heights, bottom, atol=tolerance, rtol=0)
    )
    n_top = np.count_nonzero(
        np.isclose(heights, top, atol=tolerance, rtol=0)
    )

    return n_bottom, n_top


fcc = ase_bulk(
    "Fe",
    crystalstructure="fcc",
    a=LATTICE_CONSTANT,
    cubic=True,
)

for hkl in MILLER_INDICES:
    slab = surface(
        fcc,
        indices=hkl,
        layers=LAYERS,
        vacuum=None,
        periodic=False,
    )
    slab = slab.repeat((*INPLANE_REPEATS, 1))
    slab.center(vacuum=10.0, axis=2)
    slab.pbc = (True, True, False)

    bottom, top = count_surface_atoms(slab)

    print(
        f"{hkl}: total slab atoms={len(slab)}, "
        f"bottom-face atoms={bottom}, top-face atoms={top}"
    )