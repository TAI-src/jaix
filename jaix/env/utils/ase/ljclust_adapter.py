import requests
from ase.calculators.lj import LennardJones
from ase.optimize.optimize import Optimizer
from ase.optimize import BFGS
from ase import Atoms
import os
import numpy as np
from typing import Optional
from scipy.spatial.distance import pdist

import logging
from jaix import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class LJClustAdapter:
    @staticmethod
    def _download_tar(target_dir) -> str:
        """
        Makes a request to the Cambridge database and creates a .xyz file from it.
        Taken from https://github.com/Martin092/CSE_Project/blob/main/auxiliary/cambridge_database.py
        :param atoms: Number of atoms in cluster.
        :param root: Directory root folder.
        :return: The file path.
        """
        target_path = os.path.join(target_dir, "LJ.tar")
        data_link = "https://doye.chem.ox.ac.uk/jon/structures/LJ/LJ.tar"

        # Download data if not already exists
        if not os.path.exists(target_path):
            try:
                response = requests.get(
                    data_link,
                    timeout=10,
                )
            except requests.exceptions.ConnectionError:
                raise RuntimeError(
                    "ERROR: Web request failed, please check your internet connection."
                )
            logger.debug("GET request sent to the database")
            if response.status_code != 200:
                raise RuntimeError(
                    f"ERROR: Web request failed with {response.status_code}"
                )
            with open(target_path, "wb") as file:
                file.write(response.content)
        else:
            logger.debug("Data already downloaded, skipping download step.")
        return target_path

    @staticmethod
    def _unpack_tar(tar_path: str, target_dir: str) -> str:
        """
        Unpacks the tar file to the target directory.
        :param tar_path: Path to the tar file.
        :param target_dir: Directory to unpack the tar file into.
        """
        import tarfile

        with tarfile.open(tar_path, "r") as tar:
            tar.extractall(path=target_dir)
            logger.debug(f"Unpacked {tar_path} to {target_dir}")
        # Quick and dirty workaround because file for 115 is not formatted correctly
        # Remove first line from the 115 file
        cluster_file = os.path.join(target_dir, "115")
        if os.path.exists(cluster_file):
            with open(cluster_file, "r") as file:
                lines = file.readlines()
            with open(cluster_file, "w") as file:
                file.writelines(lines[1:])
        return target_dir

    @staticmethod
    def _retrieve_cluster_data(num_atoms, target_dir: str = ".") -> np.ndarray:
        """
        Downloads and unpacks the Cambridge database of Lennard-Jones clusters.
        :return: Path to the directory containing the .xyz files.
        """
        lj_data_dir = os.path.join(target_dir, "LJ_data")
        # Check if the directory already exists
        if not os.path.exists(lj_data_dir):
            tar_path = LJClustAdapter._download_tar(target_dir)
            LJClustAdapter._unpack_tar(tar_path, target_dir=lj_data_dir)
            logger.debug("Lennard-Jones clusters database downloaded and unpacked.")
        # TODO: Identify what the "i" versions are for
        # Lowest energy icosahedral minima at sizes with non-icosahedral global minima.

        cluster_file = os.path.join(lj_data_dir, str(num_atoms))
        if not os.path.exists(cluster_file):
            logger.error(
                f"Cluster file for {num_atoms} atoms does not exist in {lj_data_dir}."
            )
            raise FileNotFoundError(
                f"Cluster file for {num_atoms} atoms does not exist in {lj_data_dir}."
            )
        logger.debug(f"Cluster file for {num_atoms} atoms found at {cluster_file}.")
        positions = np.loadtxt(cluster_file)
        assert (
            positions.shape[1] == 3
        ), "Positions should have three columns for x, y, z coordinates."
        assert (
            positions.shape[0] == num_atoms
        ), f"Expected {num_atoms} atoms, but got {positions.shape[0]}."
        return positions

    @staticmethod
    def retrieve_known_min(
        num_atoms: int, target_dir: str = ".", lj_params: dict = {}
    ) -> tuple[float, Atoms]:
        """
        Retrieves the known minimum configuration for a given number of atoms.
        :param num_atoms: Number of atoms in the cluster.
        :param target_dir: Directory to store the data.
        :return: Path to the .xyz file containing the known minimum configuration.
        """
        def_lj_params = LennardJones.default_parameters.copy()
        def_lj_params.update(lj_params)
        positions = LJClustAdapter._retrieve_cluster_data(num_atoms, target_dir)
        atoms = Atoms(
            "C" * num_atoms,  # Assuming all atoms are carbon for LJ clusters
            positions=positions,
            calculator=LennardJones(**def_lj_params),
        )
        # TODO: Do we need to specify the cell?
        # TODO: Are they always carbon atoms?
        e = atoms.get_potential_energy()
        return e, atoms

    def __init__(
        self,
        num_atoms: int,
        opt_alg: Optimizer = BFGS,
        opt_params: dict = {},
        lj_params: dict = {},
        target_dir: str = ".",
        fmax: float = 0.5,
        local_steps: int = 1000,
        covalent_radius: float = 1.0,
    ):
        """
        Initializes the LJClustAdapter with the number of atoms and optional Lennard-Jones parameters.
        :param num_atoms: Number of atoms in the cluster.
        :param lj_params: Optional dictionary of Lennard-Jones parameters.
        :param target_dir: Directory to store the data.
        """
        self.num_atoms = num_atoms
        self.lj_params = lj_params
        self.target_dir = target_dir
        self.positions = LJClustAdapter._retrieve_cluster_data(num_atoms, target_dir)
        self.min_val, self.min_pos = LJClustAdapter.retrieve_known_min(
            num_atoms, target_dir, lj_params
        )
        self.atom_str = f"C{num_atoms}"  # Assuming carbon atoms
        self.opt_alg = opt_alg
        self.opt_params = opt_params
        self.fmax = fmax
        self.local_steps = local_steps
        self.covalent_radius = covalent_radius
        self.box_length = (
            2
            * covalent_radius
            * (0.5 + ((3.0 * num_atoms) / (4 * np.pi * np.sqrt(2))) ** (1 / 3))
        )
        # TODO: Importance of box_length?

    def evaluate(self, positions: np.ndarray) -> tuple[float, np.ndarray]:
        """
        Evaluates the potential energy of a given configuration of atoms.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Tuple containing the potential energy and forces.
        """
        atoms = Atoms(
            self.atom_str,
            positions=positions,
            calculator=LennardJones(**self.lj_params),
        )
        energy = atoms.get_potential_energy()
        return energy

    def local_opt(self, positions: np.ndarray) -> tuple[float, np.ndarray]:
        """
        Performs local optimization on the given atomic positions.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :param steps: Number of optimization steps.
        :return: Optimized atomic positions.
        """

        atoms = Atoms(
            self.atom_str,
            positions=positions,
            calculator=LennardJones(**self.lj_params),
        )
        opt = self.opt_alg(atoms, **self.opt_params)
        opt.run(fmax=self.fmax, steps=self.local_steps)
        return atoms.get_potential_energy(), atoms.get_positions()

    def validate(self, positions: np.ndarray) -> bool:
        """
        Validates the atomic configuration based on physical laws.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Boolean indicating whether the configuration is valid.
        """

        if positions.shape[0] == 0:
            return True

        distances = pdist(positions)
        return float(np.min(distances)) >= 0.15

    def random_generate(self, seed: Optional[int] = None) -> np.ndarray:
        """
        Generates random atomic configurations within a specified box length.
        :param num_samples: Number of random configurations to generate.
        :param box_length: Length of the cubic box in which atoms are placed.
        :return: Numpy array of shape (num_samples, num_atoms, 3) representing the positions of the atoms.
        """
        rng = np.random.default_rng(seed)

        valid = False
        positions = None
        while not valid:
            positions = (rng.random((self.num_atoms, 3)) - 0.5) * self.box_length * 1.5
            valid = self.validate(positions)
        return positions
