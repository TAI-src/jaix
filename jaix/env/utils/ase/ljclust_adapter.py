from csv import DictReader
import requests
from ase.calculators.lj import LennardJones
from ase.optimize.optimize import Optimizer
from ase.optimize import BFGS, FIRE
from ase import Atoms
import os
import numpy as np
from typing import Optional, Type
from scipy.spatial.distance import pdist
from ttex.config import ConfigurableObject
from ase.calculators.kim import KIM
from os import path

import logging
from jaix import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


class LJClustAdapterConfig:
    def __init__(
        self,
        target_dir: str = "./ljclust_data",
        opt_alg: Type[Optimizer] = BFGS,
        opt_params: dict = {},
        fmax: float = 0.5,
        local_steps: int = 1000,
        covalent_radius: float = 1.0,
    ):
        self.target_dir = target_dir
        self.opt_alg = opt_alg
        self.opt_params = opt_params
        self.fmax = fmax
        self.local_steps = local_steps
        self.covalent_radius = covalent_radius


class LJClustAdapter(ConfigurableObject):
    config_class = LJClustAdapterConfig

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
    def retrieve_lj_params(atom_str: str) -> dict:
        """
        Retrieves the Lennard-Jones parameters for a given atom string.
        :param atom_str: Atom string in the format "C{num_atoms}" or "X{num_atoms}".
        :return: Tuple containing sigma, epsilon, and cutoff.
        """
        atom_numbers = Atoms(atom_str).get_atomic_numbers()
        assert (
            len(set(atom_numbers)) == 1
        ), "This part only supports single-element clusters."
        num_atoms = len(atom_numbers)
        if atom_str == f"X{num_atoms}":
            # Default settings from KIM model
            return {"sigma": 1.0, "epsilon": 1.0, "cutoff": 4}

        # TODO: there is probably a nicer way to do this
        material = atom_str.replace(str(num_atoms), "")  # Remove the number of atoms
        params = {}
        file_name = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "lj_params.csv"
        )
        with open(file_name, newline="") as csvfile:
            reader = DictReader(csvfile, delimiter=",")
            for row in reader:
                if row["Species_i"] == material:
                    assert (
                        row["Species_j"] == material
                    ), "This part only supports single-element clusters."
                    params = {
                        "sigma": float(row["sigma"]),
                        "epsilon": float(row["epsilon"]),
                        "cutoff": float(row["cutoff"]),
                    }
                    break
        if not params:
            raise ValueError(
                f"No Lennard-Jones parameters found for species {material}."
            )
        return params

    @staticmethod
    def finst2species(function: int, instance: int) -> str:
        """
        Converts function and instance to a species string.
        :param function: Function index.
        :param instance: Instance index.
        :return: Species string.
        """
        assert function == 0, "LJClustAdapter only supports one function (0)."
        assert instance >= 0, "Instance must be a non-negative integer."
        available_numbers = list(range(3, 151))  # Valid number of atoms
        if instance >= len(available_numbers):
            raise ValueError(
                f"Instance {instance} is out of range for LJClustAdapter. "
                f"Valid instances are from 0 to {len(available_numbers) - 1}."
            )
        num_atoms = available_numbers[instance]
        return f"C{num_atoms}"  # Species string for LJ cluster

    @staticmethod
    def _construct_atoms(positions: np.ndarray, atom_str: str) -> Atoms:
        """
        Constructs an Atoms object from the given positions.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Atoms object.
        """
        num_atoms = len(Atoms(atom_str).get_atomic_numbers())

        # TODO: figure out the cell
        if atom_str == f"X{num_atoms}":
            # No specific Material set, assuming theoretical LJ setting
            # Easiest to compute with in-built LJ and specific (near) infinite cutoff
            calc = LennardJones(
                sigma=1.0,
                epsilon=1.0,
                rc=1e6,
                smooth=False,
            )
            atoms = Atoms(
                positions=positions,
                calculator=calc,
            )
        else:
            # If atom_str is set, we assume a specific LJ setting
            # Best to compute with KIM, which knows about the LJ parameters (sigma, epsilon, cutoff)
            calc = KIM("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
            atoms = Atoms(
                atom_str,
                positions=positions,
                calculator=calc,
            )

        atoms.set_pbc(False)  # No periodic boundary conditions for LJ clusters

        return atoms

    @staticmethod
    def retrieve_known_min(
        atom_str: str,
        target_dir: str = ".",
        local_opt: bool = True,
    ) -> tuple[float, Atoms]:
        """
        Retrieves the known minimum configuration for a given number of atoms.
        :param num_atoms: Number of atoms in the cluster.
        :param target_dir: Directory to store the data.
        :return: Path to the .xyz file containing the known minimum configuration.
        """
        num_atoms = len(Atoms(atom_str).get_atomic_numbers())

        positions = LJClustAdapter._retrieve_cluster_data(num_atoms, target_dir)
        if atom_str != f"X{num_atoms}":
            # If the atom_str is not a theoretical LJ setting, we scale the positions
            # to match the covalent radius of the atoms
            try:
                lj_params = LJClustAdapter.retrieve_lj_params(atom_str)
            except AssertionError as e:
                logger.warning(
                    f"Failed to retrieve LJ parameters for {atom_str}: {e}. "
                    "Returning None"
                )
                return np.nan, None
            except ValueError as e:
                logger.warning(
                    f"Failed to retrieve LJ parameters for {atom_str}: {e}. "
                    "Returning None"
                )
                return np.nan, None

            # Scale positions based on the sigma value from the Lennard-Jones parameters
            positions *= lj_params["sigma"]

        atoms = LJClustAdapter._construct_atoms(positions, atom_str)

        if local_opt:
            # Perform local optimization to finetune optimum
            opt_alg = FIRE(atoms)
            opt_alg.run(fmax=1e-10, steps=1000)
        e = atoms.get_potential_energy()
        return e, atoms

    def __init__(self, config: LJClustAdapterConfig):
        """
        Initializes the LJClustAdapter with the number of atoms and optional Lennard-Jones parameters.
        :param num_atoms: Number of atoms in the cluster.
        :param lj_params: Optional dictionary of Lennard-Jones parameters.
        :param target_dir: Directory to store the data.
        """
        ConfigurableObject.__init__(self, config)

        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)
            logger.debug(f"Created target directory: {self.target_dir}")

    def set_species(self, species_str: str) -> None:
        """
        Sets the number of atoms for the adapter.
        :param num_atoms: Number of atoms in the cluster.
        """
        tmp_atoms = Atoms(species_str)
        # need to figure out what makes sense here
        self.num_atoms = len(tmp_atoms.get_atomic_numbers())
        # TODO: We only do carbon atm
        assert (
            species_str == f"C{self.num_atoms}"
        ), f"Species string {species_str} does not match expected format C{self.num_atoms}."
        self.min_val, self.min_atoms = LJClustAdapter.retrieve_known_min(
            species_str, self.target_dir
        )
        self.min_pos = self.min_atoms.get_positions()
        self.atom_str = species_str
        self.box_length = (
            2
            * self.covalent_radius
            * (0.5 + ((3.0 * self.num_atoms) / (4 * np.pi * np.sqrt(2))) ** (1 / 3))
        )
        # TODO: Importance of box_length?

    def evaluate(self, positions: np.ndarray) -> tuple[float, np.ndarray]:
        """
        Evaluates the potential energy of a given configuration of atoms.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Tuple containing the potential energy and forces.
        """
        atoms = LJClustAdapter._construct_atoms(positions, self.atom_str)
        energy = atoms.get_potential_energy()
        return energy, self.info(atoms)

    def info(self, atoms: Atoms) -> dict:
        """
        Returns information about the atom wrt the distance to the known minimum.
        :param atom: Atom object.
        :return: Dictionary containing information about the atom.
        """
        # TODO: Add additional info, such as distance in atom space, isomerism, etc.
        return {
            "energy": atoms.get_potential_energy() - self.min_val,
        }

    def local_opt(self, positions: np.ndarray) -> tuple[float, np.ndarray]:
        """
        Performs local optimization on the given atomic positions.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :param steps: Number of optimization steps.
        :return: Optimized atomic positions.
        """

        atoms = LJClustAdapter._construct_atoms(positions, self.atom_str)
        opt = self.opt_alg(atoms, **self.opt_params)
        opt.run(fmax=self.fmax, steps=self.local_steps)
        return atoms.get_potential_energy(), atoms.get_positions()

    def validate(self, positions: np.ndarray) -> bool:
        """
        Validates the atomic configuration based on physical laws.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Boolean indicating whether the configuration is valid.
        """
        # TODO: Why is this a good validation?

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
        # TODO: Figure out proper seeding across jaix

        valid = False
        positions: np.ndarray
        while not valid:
            positions = (rng.random((self.num_atoms, 3)) - 0.5) * self.box_length * 1.5
            valid = self.validate(positions)
        return positions
