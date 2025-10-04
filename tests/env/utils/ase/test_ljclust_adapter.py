from jaix.env.utils.ase import LJClustAdapter, LJClustAdapterConfig

if LJClustAdapter is not None:
    from ase.optimize import LBFGS
    from ase import Atoms
    from ase.optimize.optimize import Optimizer
    from ase.calculators.lj import LennardJones


import os
import shutil
import csv
import pytest
import numpy as np
from ttex.config import ConfigFactory, ConfigurableObjectFactory as COF
from ...singular.test_ljclust_env import skip_remaining_tests


target_dir = "./tmp_data"


@pytest.fixture(scope="module", autouse=True)
def data_manager():
    # create a temporary directory for the tests
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    assert os.path.exists(target_dir), "Target directory could not be created."
    yield
    # cleanup after the tests
    assert os.path.exists(target_dir), "Target directory does not exist."
    shutil.rmtree(target_dir, ignore_errors=True)


def test_download_unpack():
    # Test a separate target dir to make sure downloading works properly
    tdir = "./tmp_data2"
    shutil.rmtree(tdir, ignore_errors=True)
    os.makedirs(tdir)
    tar_file = LJClustAdapter._download_tar(tdir)
    assert os.path.exists(tar_file), "Tar file was not downloaded."
    LJClustAdapter._unpack_tar(tar_file, target_dir=os.path.join(tdir, "LJ_data"))
    assert os.path.exists(
        os.path.join(tdir, "LJ_data")
    ), "LJ_data directory was not created."

    # make sure it still works if existing data is there
    tar_file = LJClustAdapter._download_tar(tdir)
    assert os.path.exists(tar_file), "Tar file not found"
    LJClustAdapter._unpack_tar(tar_file, target_dir=os.path.join(tdir, "LJ_data"))
    assert os.path.exists(
        os.path.join(tdir, "LJ_data")
    ), "LJ_data directory was not created."
    shutil.rmtree(tdir, ignore_errors=False)


def test_retrieve_cluster_data():
    num_atoms = 13
    positions = LJClustAdapter._retrieve_cluster_data(num_atoms, target_dir=target_dir)
    assert positions.shape == (num_atoms, 3)


@pytest.mark.parametrize("material", ["He", "C", "X"])
def test_retrieve_lj_params_file(material):
    params = LJClustAdapter._retrieve_lj_params(material)
    assert isinstance(params, dict), "LJ parameters should be a dictionary."
    assert "sigma" in params, "LJ parameters should contain 'sigma'."
    assert "epsilon" in params, "LJ parameters should contain 'epsilon'."
    assert isinstance(params["sigma"], float), "'sigma' should be a float."


def test_retrieve_lj_params_file_special():
    # Unknown material, throws error
    with pytest.raises(ValueError):
        LJClustAdapter._retrieve_lj_params("Abc")
    with pytest.raises(ValueError):
        LJClustAdapter._retrieve_lj_params("C3H4")
    with pytest.raises(ValueError):
        LJClustAdapter._retrieve_lj_params("C3")
    # Return all parameters
    params = LJClustAdapter._retrieve_lj_params()
    assert isinstance(params, dict), "LJ parameters should be a dictionary."
    assert "C" in params, "LJ parameters should contain 'C'."
    assert "sigma" in params["C"], "LJ parameters for 'C' should contain 'sigma'."
    assert "epsilon" in params["C"], "LJ parameters for 'C' should contain 'epsilon'."


@pytest.mark.parametrize("atom_str", ["C3", "X17", "He33", "Abc12", "C3H4"])
def test_retrieve_lj_params(atom_str):
    params = LJClustAdapter.retrieve_lj_params(atom_str)
    if atom_str == "Abc12" or atom_str == "C3H4":
        assert params is None, "Unknown material should return None."
    else:
        assert isinstance(params, dict), "LJ parameters should be a dictionary."
        assert "sigma" in params, "LJ parameters should contain 'sigma'."
        assert "epsilon" in params, "LJ parameters should contain 'epsilon'."
        assert isinstance(params["sigma"], float), "'sigma' should be a float."


def test_finst2species():
    with pytest.raises(ValueError):
        LJClustAdapter.finst2species(-1, 0)
    with pytest.raises(ValueError):
        LJClustAdapter.finst2species(0, -1)
    with pytest.raises(ValueError):
        LJClustAdapter.finst2species(0, 148)
    with pytest.raises(ValueError):
        LJClustAdapter.finst2species(0, 148, False)
    assert LJClustAdapter.finst2species(0, 0) == "X3"
    assert LJClustAdapter.finst2species(0, 147) == "X150"
    assert LJClustAdapter.finst2species(2, 147) == "He150"
    assert LJClustAdapter.finst2species(118, 147) == "Og150"
    with pytest.raises(ValueError):
        LJClustAdapter.finst2species(119, 1)
    with pytest.raises(ValueError):
        LJClustAdapter.finst2species(1, 119, False)


def create_random_positions(num_atoms):
    covalent_radius = 1
    box_length = (
        2
        * covalent_radius
        * (0.5 + ((3.0 * num_atoms) / (4 * np.pi * np.sqrt(2))) ** (1 / 3))
    )
    rng = np.random.default_rng(42)
    positions = (rng.random((num_atoms, 3)) - 0.5) * box_length * 1.5
    return positions


def test_validate():
    # TODO: once I know what we actually want to validate
    pass


@pytest.mark.parametrize("atom_str", ["C3", "X13", "He150", "C3H4"])
def test_construct_atoms(atom_str):
    num_atoms = len(Atoms(atom_str).get_atomic_numbers())
    positions = create_random_positions(num_atoms)
    atoms = LJClustAdapter._construct_atoms(positions, atom_str)
    assert isinstance(atoms, Atoms), "Constructed atoms is not an instance of Atoms."
    assert isinstance(
        atoms.get_potential_energy(), float
    ), "Potential energy is not a float."
    if atom_str.startswith("X"):
        assert isinstance(
            atoms.calc, LennardJones
        ), "Calculator is not an instance of LennardJones."
    else:
        assert not isinstance(
            atoms.calc, LennardJones
        ), "Calculator is not an instance of LennardJones."


def test_retrieve_known_min_X(request):
    with open(request.path.parent.joinpath("glob_min.csv"), newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        glob_min = {int(row["N"]): float(row["Energy"]) for row in reader}
    for num_atoms in range(3, 151):
        min_val, atoms = LJClustAdapter.retrieve_known_min(
            f"X{num_atoms}",
            target_dir=target_dir,
            local_opt=False,
        )
        # Check that the retrieved minimum roughly matches the known minimum

        assert (
            abs(min_val - glob_min[num_atoms]) < 1e-5
        ), f"Mismatch for {num_atoms} atoms: {min_val} != {glob_min[num_atoms]}"
        # Check that the positions are valid
        assert atoms.positions.shape == (
            num_atoms,
            3,
        ), f"Positions for {num_atoms} atoms do not match expected shape."
        assert LJClustAdapter.validate(
            atoms.positions
        ), f"Positions for {num_atoms} atoms are not valid."


@pytest.mark.parametrize("atom_str", ["C3", "X17", "He33", "Og150", "C150"])
def test_retrieve_known_min_material(atom_str):
    min_val, atoms = LJClustAdapter.retrieve_known_min(
        atom_str,
        target_dir=target_dir,
        local_opt=False,
    )
    from ase.optimize import FIRE

    opt = FIRE(atoms, logfile=None)
    opt.run(fmax=1e-6, steps=1000)
    assert (
        abs(min_val - atoms.get_potential_energy()) < 0.005
    ), f"Mismatch for {atom_str}: {min_val} != {atoms.get_potential_energy()}"
    assert atoms.get_potential_energy() <= min_val

    # TODO: What are good target accuracies for non-X?

    min_val_loc, atoms2 = LJClustAdapter.retrieve_known_min(
        atom_str,
        target_dir=target_dir,
        local_opt=True,
    )
    assert (
        min_val_loc <= min_val
    ), f"Local minimum for {atom_str} did not improve after minimization"
    assert (
        abs(min_val - atoms2.get_potential_energy()) < 0.005
    ), f"Mismatch for {atom_str}: {min_val} != {atoms2.get_potential_energy()}"


@pytest.mark.parametrize("atom_str", ["C3H4", "Abc17", "Dum3"])
def test_retrieve_known_min_invalid(atom_str):
    try:
        Atoms(atom_str)  # Check if the atom string is valid
    except ValueError:
        with pytest.raises(ValueError):
            LJClustAdapter.retrieve_known_min(
                atom_str,
                target_dir=target_dir,
                local_opt=False,
            )
    else:
        # If the atom string is valid, we should not raise an error
        min_val, atoms = LJClustAdapter.retrieve_known_min(
            atom_str,
            target_dir=target_dir,
            local_opt=False,
        )
        assert min_val is np.nan, f"Expected None for {atom_str}, got {min_val}."
        assert atoms is None, f"Expected None for {atom_str}, got {atoms}."


def get_config(def_vals: bool) -> LJClustAdapterConfig:
    adapter_params = {
        "target_dir": target_dir,
    }
    if not def_vals:
        spec_params = {
            "opt_alg": LBFGS,
            "opt_alg_params": {"alpha": 71},
            "opt_run_params": {"fmax": 0.1, "steps": 1000},
        }
        adapter_params.update(spec_params)
    config = LJClustAdapterConfig(**adapter_params)
    return config


@pytest.mark.parametrize("def_vals", [True, False])
def test_init(def_vals):
    # Test initialization of the adapter
    adapter_config = get_config(def_vals=def_vals)
    adapter = LJClustAdapter(adapter_config)
    assert isinstance(
        adapter, LJClustAdapter
    ), "Adapter is not an instance of LJClustAdapter."
    for key, value in adapter_config.__dict__.items():
        assert (
            adapter.__dict__[key] == value
        ), f"{key} is not set correctly in the adapter."


def test_init_advanced():
    # Test that we can pass optimizer as a string
    tdir = "./tmp_data2"
    config = get_config(def_vals=False)
    adapter_params = config.__dict__
    adapter_params["target_dir"] = tdir
    adapter_params["opt_alg"] = "ase.optimize.LBFGS"
    config_dict = {"jaix.env.utils.ase.LJClustAdapterConfig": adapter_params}
    nconfig = ConfigFactory.from_dict(config_dict)
    assert issubclass(
        nconfig.opt_alg, Optimizer
    ), "Optimizer is not an subclass of ase.optimize.Optimizer."

    # Test folder creation works
    adapter = LJClustAdapter(nconfig)
    assert isinstance(
        adapter, LJClustAdapter
    ), "Adapter is not an instance of LJClustAdapter."
    assert os.path.exists(tdir), "Target directory was not created."
    shutil.rmtree(tdir, ignore_errors=True)


def test_set_species():
    # Test setting species works
    adapter_config = get_config(def_vals=True)
    adapter = LJClustAdapter(adapter_config)
    adapter.set_species("Ar13C3")
    assert adapter.num_atoms == 16, "Number of atoms was not set to 16."
    assert adapter.atom_str == "Ar13C3", "Atom string was not set to 'Ar13C3'."
    assert adapter.min_pos is None

    adapter.set_species("C5")
    assert adapter.num_atoms == 5, "Number of atoms was not set to 5."
    assert adapter.atom_str == "C5", "Atom string was not set to '5C'."

    assert isinstance(adapter.min_val, float), "Minimum value is not a float."
    assert isinstance(adapter.box_length, float), "Box length is not a float."
    assert isinstance(
        adapter.min_pos, np.ndarray
    ), "Minimum positions are not a numpy array."


@pytest.mark.parametrize(
    "def_vals,atom_str", [(True, "C13"), (False, "He33"), (True, "C3H4"), (True, "X17")]
)
def test_generate(def_vals, atom_str):
    # Test the generate method of the adapter
    adapter_config = get_config(def_vals=def_vals)
    adapter = LJClustAdapter(adapter_config)
    adapter.set_species(atom_str)  # Set number of atoms to 13
    pos = adapter.random_generate()
    assert pos.shape == (
        adapter.num_atoms,
        3,
    ), "Generated positions do not match the number of atoms."
    # Check that the positions are valid
    assert adapter.validate(pos), "Generated positions are not valid."
    # Check we can create an atoms object
    atoms = LJClustAdapter._construct_atoms(pos, atom_str)
    assert isinstance(
        atoms, Atoms
    ), "Atoms object could not be created from generated positions."
    assert isinstance(
        atoms.get_potential_energy(), float
    ), "Potential energy could not be calculated from generated positions."


@pytest.mark.parametrize(
    "def_vals,atom_str", [(True, "C13"), (False, "He33"), (True, "C3H4"), (True, "X17")]
)
def test_evaluate(def_vals, atom_str):
    # Test the evaluate method of the adapter
    adapter_config = get_config(def_vals=def_vals)
    adapter = LJClustAdapter(adapter_config)
    adapter.set_species(atom_str)  # Set number of atoms to 13
    pos = adapter.random_generate()
    energy, info = adapter.evaluate(pos)
    assert not np.isnan(energy), "Energy should not be NaN."
    assert isinstance(energy, float), "Energy is not a float."
    assert "energy_diff" in info, "Info dictionary does not contain 'energy' key."
    assert isinstance(
        info["energy_diff"], float
    ), "Energy in info dictionary is not a float."
    if atom_str == "C3H4":
        assert np.isnan(
            info["energy_diff"]
        ), "Energy difference for C3H4 should be NaN."
    else:
        assert not np.isnan(
            info["energy_diff"]
        ), "Energy difference should not be NaN for other materials."


@pytest.mark.parametrize(
    "def_vals,atom_str", [(True, "C13"), (False, "He33"), (True, "C3H4"), (True, "X17")]
)
def test_local_opt(def_vals, atom_str):
    # Test the local optimization method of the adapter
    adapter_config = get_config(def_vals=def_vals)
    adapter = LJClustAdapter(adapter_config)
    adapter.set_species(atom_str)  # Set number of atoms to 13
    pos = adapter.random_generate()
    # Get current energy
    initial_energy, _ = adapter.evaluate(pos)
    energy, opt_pos = adapter.local_opt(pos)
    assert isinstance(energy, float), "Energy after optimization is not a float."
    assert opt_pos.shape == (
        adapter.num_atoms,
        3,
    ), "Optimized positions do not match the number of atoms."
    # Check that the optimized positions are valid
    assert adapter.validate(opt_pos), "Optimized positions are not valid."
    # Check that the energy is lower than the initial energy
    assert (
        energy <= initial_energy
    ), "Energy after optimization is not lower than initial energy."
