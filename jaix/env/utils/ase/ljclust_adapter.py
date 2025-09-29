from csv import DictReader
import requests
from ase.calculators.lj import LennardJones
from ase.optimize.optimize import Optimizer
from ase.optimize import BFGS, FIRE
from ase import Atoms
import os
import numpy as np
from typing import Optional, Type, Union, List, Tuple
from scipy.spatial.distance import pdist
from ttex.config import ConfigurableObject, Config
from ase.calculators.kim import KIM
from os import path

import logging
from jaix.utils.globals import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)

LJ_MAX_RC = 1e6  # near infinite cutoff to get as close to the optimum values known in literature as possible
MIN_ATOM_DISTANCE = 0.15  # Minimum distance between atoms in the Lennard-Jones clusters


class LJClustAdapterConfig(Config):
    def __init__(
        self,
        target_dir: str = "./ljclust_data",
        opt_alg: Type[Optimizer] = BFGS,
        opt_alg_params: Optional[
            dict
        ] = None,  # Parameters for the optimizer initialization
        opt_run_params: Optional[
            dict
        ] = None,  # Parameters for the optimizer run, including fmax and steps
    ):
        self.target_dir = target_dir
        self.opt_alg = opt_alg
        self.opt_alg_params = {} if opt_alg_params is None else opt_alg_params
        self.opt_run_params = {} if opt_run_params is None else opt_run_params


class LJClustAdapter(ConfigurableObject):
    config_class = LJClustAdapterConfig

    @staticmethod
    def _download_tar(target_dir) -> str:
        """
        Downloads the Lennard-Jones clusters data from the Cambridge database.
        :param target_dir: Directory to store the downloaded data.
        :return: Path to the downloaded tar file.
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
        Unpacks the downloaded tar file containing Lennard-Jones clusters data.
        :param tar_path: Path to the downloaded tar file.
        :param target_dir: Directory to unpack the data into.
        :return: Path to the directory containing the unpacked data.
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
        Retrieves the cluster data for a given number of atoms from the Lennard-Jones clusters database.
        :param num_atoms: Number of atoms in the cluster.
        :param target_dir: Directory to store the data.
        :return: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        """
        lj_data_dir = os.path.join(target_dir, "LJ_data")
        # Check if the directory already exists
        if not os.path.exists(lj_data_dir):
            tar_path = LJClustAdapter._download_tar(target_dir)
            LJClustAdapter._unpack_tar(tar_path, target_dir=lj_data_dir)
            logger.debug("Lennard-Jones clusters database downloaded and unpacked.")
        # i versions are for Lowest energy icosahedral minima at sizes with non-icosahedral global minima.

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
    def _retrieve_lj_params(material: Optional[str] = None) -> dict:
        """
        Retrieves the Lennard-Jones parameters for a given material from a CSV file.
        :param material: Material for which to retrieve the parameters. If None, retrieves all parameters.
        :return: Dictionary containing Lennard-Jones parameters (sigma, epsilon, cutoff).
        """
        if material == "X":
            # Special case for theoretical LJ clusters, return default parameters
            return {"sigma": 1.0, "epsilon": 1.0, "cutoff": 4.0}

        file_name = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "lj_params.csv"
        )
        if not os.path.exists(file_name):
            raise FileNotFoundError(
                f"Lennard-Jones parameters file {file_name} not found."
            )

        params = {}
        with open(file_name, newline="") as csvfile:
            reader = DictReader(csvfile, delimiter=",")
            if material is None:
                return {
                    row["Species_i"]: {
                        "sigma": float(row["sigma"]),
                        "epsilon": float(row["epsilon"]),
                        "cutoff": float(row["cutoff"]),
                    }
                    for row in reader
                }
            # If material is specified, retrieve parameters for that material
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
                f"Lennard-Jones parameters for {material} not found in {file_name}."
            )
        logger.debug(f"Retrieved Lennard-Jones parameters for {material}: {params}")
        return params

    @staticmethod
    def retrieve_lj_params(atom_str: str) -> Union[dict, None]:
        """
        Retrieves the Lennard-Jones parameters for a given atom string.
        :param atom_str: String representing the atoms, e.g., "Cu10" for 10 copper atoms.
        :return: Dictionary containing Lennard-Jones parameters (sigma, epsilon, cutoff).
        """
        try:
            atom_numbers = Atoms(atom_str).get_atomic_numbers()
        except ValueError as e:
            logger.error(f"Failed to parse atom string {atom_str}: {e}")
            return None
        if len(set(atom_numbers)) > 1:
            logger.error(
                "LJClustAdapter only supports single-element clusters, "
                f"but found multiple elements in {atom_str}."
            )
            return None
        num_atoms = len(atom_numbers)
        material = atom_str.replace(str(num_atoms), "")  # Remove the number of atoms
        try:
            params = LJClustAdapter._retrieve_lj_params(material)
        except FileNotFoundError as e:
            logger.error(f"Failed to retrieve LJ parameters file: {e}")
            return None
        except ValueError as e:
            logger.error(f"Failed to retrieve LJ parameters for {material}: {e}")
            return None
        return params

    @staticmethod
    def get_info(by_species: bool):
        """
        Retrieves available species and numbers for Lennard-Jones clusters.
        :param by_species: If True, returns species as functions; otherwise, returns numbers as functions.
        :return: Dictionary containing available species, available numbers, number of functions, and number of instances.
        """
        params = LJClustAdapter._retrieve_lj_params()  # Retrieve all parameters
        available_species: List[str] = ["X"] + list(
            params.keys()
        )  # Add "X" for theoretical LJ clusters
        available_numbers: List[int] = list(range(3, 151))  # Valid number of atoms
        functions = available_species if by_species else available_numbers
        instances = available_numbers if by_species else available_species
        return {
            "avail_species": available_species,
            "avail_numbers": available_numbers,
            "num_funcs": len(functions),  # type: ignore[arg-type]
            "num_insts": len(instances),  # type: ignore[arg-type]
        }

    @staticmethod
    def finst2species(function: int, instance: int, by_species: bool = True) -> str:
        """
        Converts function and instance indices to a species string for LJ clusters.
        :param function: Index of the function (species).
        :param instance: Index of the instance (number of atoms).
        :param by_species: If True, function is species index, instance is number of atoms; otherwise, vice versa.
        :return: Species string, e.g., "Cu10" for 10 copper atoms.
        """
        available_info = LJClustAdapter.get_info(by_species)
        available_species = available_info["avail_species"]
        available_numbers = available_info["avail_numbers"]

        species_idx = function if by_species else instance
        instance_idx = instance if by_species else function

        if species_idx < 0 or species_idx >= len(available_species):
            raise ValueError(
                f"Function index {function} is out of range for LJClustAdapter."
            )
        species = available_species[species_idx]
        if instance_idx < 0 or instance_idx >= len(available_numbers):
            raise ValueError(
                f"Instance index {instance} is out of range for LJClustAdapter."
            )
        num_atoms = available_numbers[instance_idx]

        return f"{species}{num_atoms}"  # Species string for LJ cluster

    @staticmethod
    def validate(positions: np.ndarray) -> bool:
        """
        Validates the atomic configuration based on physical laws.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Boolean indicating whether the configuration is valid.
        """
        # TODO: Why is this a good validation?

        if positions.shape[0] == 0:
            return True

        distances = pdist(positions)
        return float(np.min(distances)) >= MIN_ATOM_DISTANCE

    @staticmethod
    def _construct_atoms(positions: np.ndarray, atom_str: str) -> Atoms:
        """
        Constructs an ASE Atoms object from positions and atom string.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :param atom_str: String representing the atoms, e.g., "Cu10" for 10 copper atoms.
        :return: ASE Atoms object with the specified positions and calculator.
        """
        assert LJClustAdapter.validate(positions), "Invalid atomic configuration."
        num_atoms = len(Atoms(atom_str).get_atomic_numbers())

        # TODO: figure out the cell
        if atom_str == f"X{num_atoms}":
            # No specific Material set, assuming theoretical LJ setting
            # Easiest to compute with in-built LJ and specific (near) infinite cutoff
            calc = LennardJones(
                sigma=1.0,
                epsilon=1.0,
                rc=LJ_MAX_RC,  # near infinite cutoff
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
    ) -> Tuple[float, Union[Atoms, None]]:
        """
        Retrieves the known minimum energy and corresponding atomic positions for a given atom string.
        :param atom_str: String representing the atoms, e.g., "Cu10" for 10 copper atoms.
        :param target_dir: Directory to store the data.
        :param local_opt: Whether to perform local optimization on the retrieved positions.
        :return: Tuple containing the minimum energy and ASE Atoms object with the positions.
        """
        num_atoms = len(Atoms(atom_str).get_atomic_numbers())

        positions = LJClustAdapter._retrieve_cluster_data(num_atoms, target_dir)
        try:
            lj_params = LJClustAdapter.retrieve_lj_params(atom_str)
        except AssertionError as e:
            logger.error(f"Failed to retrieve LJ parameters for {atom_str}: {e}")
            return np.nan, None
        if lj_params is None:
            logger.warning(
                f"Lennard-Jones parameters for {atom_str} not found. " "Returning None."
            )
            return np.nan, None

        # Scale positions based on the sigma value from the Lennard-Jones parameters
        positions *= lj_params["sigma"]
        atoms = LJClustAdapter._construct_atoms(positions, atom_str)

        if local_opt:
            # Perform local optimization to finetune optimum
            opt_alg = FIRE(atoms, logfile=None)
            fmax = (
                1e-10 * lj_params["epsilon"] / lj_params["sigma"]
            )  # Force convergence threshold
            opt_alg.run(fmax=fmax, steps=1000)
        energy = atoms.get_potential_energy()
        return energy, atoms

    def __init__(self, config: LJClustAdapterConfig):
        """
        Initializes the LJClustAdapter with the given configuration.
        :param config: Configuration object containing parameters for the adapter.
        """
        ConfigurableObject.__init__(self, config)

        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)
            logger.debug(f"Created target directory: {self.target_dir}")

    def set_species(self, atom_str: str) -> None:
        """
        Sets the species for the LJClustAdapter, and corresponding values
        :param atom_str: String representing the atoms, e.g., "Cu10" for 10 copper atoms.
        """
        self.atom_str = atom_str
        self.min_val, self.min_atoms = LJClustAdapter.retrieve_known_min(
            self.atom_str, self.target_dir
        )
        self.num_atoms = len(Atoms(self.atom_str).get_atomic_numbers())
        self.min_pos = None
        if self.min_atoms is not None:
            self.min_pos = self.min_atoms.get_positions()
        # box length source? https://www.researchgate.net/publication/235583835_Local_structures_in_medium-sized_Lennard-Jones_clusters_Monte_Carlo_simulations
        self.lj_params = LJClustAdapter.retrieve_lj_params(self.atom_str)
        covalent_radius = (
            self.lj_params.get("sigma", 1.0) if self.lj_params is not None else 1.0
        )
        # TODO: figure out good values here
        self.box_length = (
            2
            * covalent_radius
            * (0.5 + ((3.0 * self.num_atoms) / (4 * np.pi * np.sqrt(2))) ** (1 / 3))
        )
        # TODO: Importance of box_length?

    def evaluate(self, positions: np.ndarray) -> Tuple[float, dict]:
        """
        Evaluates the potential energy of the given atomic positions.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Tuple containing the potential energy and a dictionary with information about the atom.
        """
        atoms = LJClustAdapter._construct_atoms(positions, self.atom_str)
        energy = atoms.get_potential_energy()
        return energy, self.info(atoms)

    def info(self, atoms: Atoms) -> dict:
        """
        Returns a dictionary with information about the atom configuration.
        :param atoms: ASE Atoms object containing the atomic configuration.
        :return: Dictionary containing information about the atom configuration.
        """
        # TODO: Add additional info, such as distance in atom space, isomerism, etc.
        return {
            "atom_str": self.atom_str,
            "num_atoms": self.num_atoms,
            "min_val": self.min_val,
            "min_pos": self.min_pos,
            "energy_diff": atoms.get_potential_energy() - self.min_val,
        }

    def local_opt(self, positions: np.ndarray) -> Tuple[float, np.ndarray]:
        """
        Performs local optimization on the given atomic positions.
        :param positions: Numpy array of shape (num_atoms, 3) representing the positions of the atoms.
        :return: Tuple containing the potential energy and the optimized positions.
        """

        atoms = LJClustAdapter._construct_atoms(positions, self.atom_str)
        opt = self.opt_alg(atoms, **self.opt_alg_params, logfile=None)
        opt.run(**self.opt_run_params)
        return atoms.get_potential_energy(), atoms.get_positions()

    def random_generate(self, seed: Optional[int] = None) -> np.ndarray:
        """
        Generates a random atomic configuration within the specified box length.
        :param seed: Optional seed for random number generation.
        :return: Numpy array of shape (num_atoms, 3) representing the random positions of the atoms.
        """
        rng = np.random.default_rng(seed)
        # TODO: Figure out proper seeding across jaix

        valid = False
        positions: np.ndarray
        while not valid:
            positions = (rng.random((self.num_atoms, 3)) - 0.5) * self.box_length * 1.5
            valid = self.validate(positions)
        return positions
