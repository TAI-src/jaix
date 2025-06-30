from jaix.env.utils.ase import LJClustAdapter
import os
import shutil
import csv
import pytest
from ase.optimize import LBFGS
from ase import Atoms
from ase.calculators.lj import LennardJones
import numpy as np

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
    assert os.path.exists(os.path.join(tdir, "LJ_data")), "LJ_data directory was not created."

    # make sure it still works if existing data is there
    tar_file = LJClustAdapter._download_tar(tdir)
    assert os.path.exists(tar_file), "Tar file not found"
    LJClustAdapter._unpack_tar(tar_file, target_dir=os.path.join(tdir, "LJ_data"))
    assert os.path.exists(os.path.join(tdir, "LJ_data")), "LJ_data directory was not created."
    shutil.rmtree(tdir, ignore_errors=False)


def test_retrieve_cluster_data():
    num_atoms = 13
    positions = LJClustAdapter._retrieve_cluster_data(num_atoms, target_dir=target_dir)
    assert positions.shape == (num_atoms, 3)


def test_retrieve_known_min(request):
    with open(request.path.parent.joinpath("glob_min.csv"), newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter=";")
        glob_min = {int(row["N"]): float(row["Energy"]) for row in reader}
    for num_atoms in range(3, 151):
        min_val, _ = LJClustAdapter.retrieve_known_min(num_atoms, target_dir=target_dir, lj_params = {"sigma": 1.0, "epsilon": 1.0, "rc":1e6, "smooth":False} )
        # Check that the retrieved minimum roughly matches the known minimum
        assert abs(min_val - glob_min[num_atoms]) < 1e-5, f"Mismatch for {num_atoms} atoms: {min_val} != {glob_min[num_atoms]}"


def get_parameters(def_vals:bool) -> dict:
    adapter_params = {
        "target_dir": target_dir,
    }
    if not def_vals:
        spec_params = {
            "opt_alg": LBFGS,
            "opt_params": {"alpha": 71},
            "lj_params": {"sigma": 1.1, "epsilon": 1.1},
            "fmax": 0.1,
            "local_steps": 1001,
            "covalent_radius": 1.1,
        }
        adapter_params.update(spec_params)
    return adapter_params

@pytest.mark.parametrize("def_vals", [True, False])
def test_init(def_vals):
    # Test initialization of the adapter
    adapter_params = get_parameters(def_vals=def_vals)
    adapter = LJClustAdapter(**adapter_params)
    assert isinstance(adapter, LJClustAdapter), "Adapter is not an instance of LJClustAdapter."
    for key, value in adapter_params.items():
        assert adapter.__dict__[key] == value, f"{key} is not set correctly in the adapter."

def test_init_advanced():
    # Test folder generation works
    # Test opt_alg extraction works
    tdir = "./tmp_data2"
    adapter_params = get_parameters(def_vals=False)
    adapter_params["target_dir"] = tdir
    adapter_params["opt_alg"] = "ase.optimize.LBFGS"
    adapter = LJClustAdapter(**adapter_params)
    assert isinstance(adapter, LJClustAdapter), "Adapter is not an instance of LJClustAdapter."
    assert os.path.exists(tdir), "Target directory was not created."
    shutil.rmtree(tdir, ignore_errors=True)

def test_set_species():
    # Test setting species works
    adapter_params = get_parameters(def_vals=True)
    adapter = LJClustAdapter(**adapter_params)
    with pytest.raises(AssertionError):
        # Should raise an error if species is not set
        adapter.set_species("Ar13C3")
    adapter.set_species("C5")
    assert adapter.num_atoms == 5, "Number of atoms was not set to 5."
    assert adapter.atom_str == "C5", "Atom string was not set to '5C'."

    assert isinstance(adapter.min_val, float), "Minimum value is not a float."
    assert isinstance(adapter.box_length, float), "Box length is not a float."
    assert isinstance(adapter.min_pos, np.ndarray), "Minimum positions are not a numpy array."
  
# TODO: Figure out a way to test validation

@pytest.mark.parametrize("def_vals", [True, False])
def test_generate(def_vals):
    # Test the generate method of the adapter
    adapter_params = get_parameters(def_vals=def_vals)
    adapter = LJClustAdapter(**adapter_params)
    adapter.set_species("C13")  # Set number of atoms to 13
    pos = adapter.random_generate()
    assert pos.shape == (adapter.num_atoms, 3), "Generated positions do not match the number of atoms."
    # Check that the positions are valid
    assert adapter.validate(pos), "Generated positions are not valid."
    # Check we can create an atoms object
    atoms = Atoms(
        adapter.atom_str,
        positions=pos,
        calculator=LennardJones(**adapter.lj_params),
    )
    assert isinstance(atoms, Atoms), "Atoms object could not be created from generated positions."
    assert isinstance(atoms.get_potential_energy(), float), "Potential energy could not be calculated from generated positions."

@pytest.mark.parametrize("def_vals", [True, False])
def test_evaluate(def_vals):
    # Test the evaluate method of the adapter
    adapter_params = get_parameters(def_vals=def_vals)
    adapter = LJClustAdapter(**adapter_params)
    adapter.set_species("C13")  # Set number of atoms to 13
    pos = adapter.random_generate()
    energy, info = adapter.evaluate(pos)
    assert isinstance(energy, float), "Energy is not a float."
    assert "energy" in info, "Info dictionary does not contain 'energy' key."
    assert isinstance(info["energy"], float), "Energy in info dictionary is not a float."
    assert info["energy"] >= 0

@pytest.mark.parametrize("def_vals", [True, False])
def test_local_opt(def_vals):
    # Test the local optimization method of the adapter
    adapter_params = get_parameters(def_vals=def_vals)
    adapter = LJClustAdapter(**adapter_params)
    adapter.set_species("C13")  # Set number of atoms to 13
    pos = adapter.random_generate()
    # Get current energy
    initial_energy, _ = adapter.evaluate(pos)
    energy, opt_pos = adapter.local_opt(pos)
    assert isinstance(energy, float), "Energy after optimization is not a float."
    assert opt_pos.shape == (adapter.num_atoms, 3), "Optimized positions do not match the number of atoms."
    # Check that the optimized positions are valid
    assert adapter.validate(opt_pos), "Optimized positions are not valid."
# Check that the energy is lower than the initial energy
    assert energy <= initial_energy, "Energy after optimization is not lower than initial energy."

def test_finst2species():
    with pytest.raises(AssertionError):
        LJClustAdapter.finst2species(2, 3)
    with pytest.raises(ValueError):
        LJClustAdapter.finst2species(0, 148)
    assert LJClustAdapter.finst2species(0, 0) == "C3"
    assert LJClustAdapter.finst2species(0, 147) == "C150"
