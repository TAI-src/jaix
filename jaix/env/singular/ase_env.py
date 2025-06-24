from ase import Atoms, Atom
from ase.optimize import BFGS
from ase.calculators.lj import LennardJones
import numpy as np
from typing import Optional, Tuple, Any, Literal
from scipy.spatial.distance import pdist


def quick_atom_test():
    rng= np.random.default_rng()

    atom_str = "C13"
    covalent_radius: float = 1.0
    num_atoms = len(Atoms(atom_str))
    box_length: float = (2 *covalent_radius * (0.5 + ((3.0 * num_atoms) / (4 * np.pi * np.sqrt(2))) ** (1 / 3)) )
    num_clusters = 4

    positions = (rng.random((num_atoms, 3)) - 0.5) * box_length * 1.5
    while not configuration_validity(positions):
        positions = (
            (rng.random((num_atoms, 3)) - 0.5) * box_length * 1.5
        )
    atoms = Atoms(
        atom_str,
        positions=positions,
        cell=np.array(
            [
                [box_length, 0, 0],
                [0, box_length, 0],
                [0, 0, box_length],
            ]
        ),
        calculator=LennardJones(),
    )
    e = atoms.get_potential_energy()
    print('Energy', e)
    f = atoms.get_forces()
    print('Forces')
    print(f)

    loc_opt = BFGS(atoms)
    e = atoms.get_potential_energy()
    print('Energy', e)
    f = atoms.get_forces()
    print('Forces')
    print(f)



def configuration_validity(
    positions: np.ndarray[Tuple[Any, Literal[3]], np.dtype[np.float64]]
) -> bool:
    """
    Checks if a potential configuration doesn't invalidate the physical laws.
    :param positions: Numpy array of the potential atomic configuration.
    :return: Boolean indicating stability of configuration.
    """
    if positions.shape[0] == 0:
        return True
    distances = pdist(positions)
    return bool(float(np.min(distances)) >= 0.15)
