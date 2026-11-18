import shutil
import pickle

import pytest
from jaix.env.utils.problem.re_problem.reproblem_adapter import (
    REProblem,
    REProblemConfig,
)

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.problems.problem_info import ProblemInfo


def sim_run(batch: Batch, incomplete=False):
    # Simulate running the batch by creating a result file in the output directory
    batch.out_dir.mkdir(parents=True, exist_ok=True)
    with open(batch.out_dir / "result.txt", "w") as f:
        f.write("done")
    if incomplete:
        # Create an incomplete run by not creating the batch file
        return
    # Create the batch file to simulate a completed batch
    batch.write()
    assert batch.exists()  # Check if the batch now exists (completed run)


def test_batch_creation(tmp_path):
    problem = REProblem(REProblemConfig(), 0)
    pinfo = ProblemInfo.from_problem(problem)
    batch = Batch(
        problem=problem,
        pinfo=pinfo,
        pid=pinfo.uuid,
        sid=5,
        algorithm_name="NSGA2",
        selector="tournament",
        n_gen=100,
        parent_dir=tmp_path,
        seed=42,  # Set a seed for reproducibility
    )
    assert batch.name == f"p{pinfo.uuid}_s5"
    settings = batch.settings
    assert ["algorithm_name", "selector", "n_gen"] == list(settings.keys())
    assert batch.out_dir.exists()  # Check if the output directory is create
    assert batch.run_id.startswith(f"{pinfo.uuid}_")  # Check if the run ID
    assert batch.run_id.endswith(
        f"_{batch.seed}"
    )  # Check if the run ID ends with the seed
    # check that the folder contains the run_id as a subfolder
    assert (tmp_path / batch.run_id).exists()  # Check if the run
    exp_id = batch.experiment_id
    assert exp_id.startswith(f"{pinfo.uuid}_")  # Check if the experiment ID

    batch2 = Batch(
        problem=problem,
        pinfo=pinfo,
        pid=pinfo.uuid,
        sid=5,
        rep=1,
        algorithm_name="NSGA2",
        selector="tournament",
        n_gen=100,
        parent_dir=tmp_path,
        seed=43,  # Set a different seed for the second batch
    )
    assert batch2.name == f"p{pinfo.uuid}_s5_r1"
    assert (
        batch.experiment_id == batch2.experiment_id
    )  # Experiment ID should be the same for different reps
    assert (
        batch2.run_id != batch.run_id
    )  # Run ID should be different for different seeds
    assert (
        batch2.out_dir.exists()
    )  # Check if the output directory is created for the second batch

    _, seed = batch2.parse_name(batch.run_id)
    assert seed == batch.seed  # Check if the parsed seed matches the original seed

    assert batch.get_num_runs() == 2  # batch 1 and 2
    # create folder
    assert batch.out_dir.exists()  # Check if the output directory exists
    assert not batch.exists()  # Check if the batch does not exist yet (empty dir)

    for b in [batch, batch2]:
        # Simulate half-finished batches by creating a dummy result file in the output directory
        sim_run(b, incomplete=True)
    run_seeds = batch.get_run_seeds()
    assert len(run_seeds) == 0  # No completed runs yet
    assert (
        not batch.exists()
    )  # Check if the batch does not exist yet (no completed runs)

    for b in [batch, batch2]:
        # Create the batch file to simulate a completed batch
        sim_run(b, incomplete=False)

    run_seeds = batch.get_run_seeds()
    assert set(run_seeds) == {batch.seed, batch2.seed}  # Check if the run seeds
    assert batch.exists()
    shutil.rmtree(batch.out_dir)  # Remove the output directory for batch 1
    assert not batch.exists()  # Now the batch should not exist

    batch3 = Batch(
        problem=problem,
        pinfo=pinfo,
        pid=pinfo.uuid,
        sid=5,
        rep=2,
        algorithm_name="NSGA3",
        selector="tournament",
        n_gen=100,
        parent_dir=tmp_path,
        seed=44,  # Set a different seed for the third batch
    )
    assert (
        batch3.experiment_id != batch.experiment_id
    )  # Experiment ID should be the same for different reps


def test_batch_missing_values(tmp_path):
    problem = REProblem(REProblemConfig(), 0)
    pinfo = ProblemInfo.from_problem(problem)
    batch = Batch(
        problem=problem,
        pinfo=pinfo,
        pid=pinfo.uuid,
        sid=5,
        algorithm_name="NSGA2",
        selector="tournament",
        n_gen=100,
        parent_dir=tmp_path,
    )
    assert batch.seed is None
    with pytest.raises(ValueError):
        _ = batch.run_id  # Should raise ValueError because seed is None


@pytest.mark.parametrize("legacy", [True, False])
def test_load(tmp_path, legacy):
    problem = REProblem(REProblemConfig(), 0)
    pinfo = ProblemInfo.from_problem(problem)
    batch = Batch(
        problem=problem,
        pinfo=pinfo,
        pid=pinfo.uuid,
        sid=5,
        algorithm_name="NSGA2",
        selector="tournament",
        n_gen=100,
        parent_dir=tmp_path,
        seed=42,  # Set a seed for reproducibility
    )
    if legacy:
        batch_file = batch.batch_file.with_suffix(".pkl")
        # Previously saved as pickle, now saved as json. Simulate a legacy batch by saving it as a pickle file.
        with open(batch_file, "wb") as f:
            b_copy = batch.model_copy(update={"problem": str(batch.problem)})
            pickle.dump(b_copy, f)
    else:
        batch_file = batch.batch_file
        batch.write()  # Write the batch to create the batch file
    assert batch.exists()  # Check if the batch now exists (completed run)
    b2 = Batch.load(batch_file)  # Load the batch from the file
    assert b2.name == batch.name  # Check if the loaded batch has the same name
    assert (
        b2.settings == batch.settings
    )  # Check if the loaded batch has the same settings
    assert b2.pid == batch.pid  # Check if the loaded batch has the same pid
    assert b2.sid == batch.sid  # Check if the loaded batch has the same sid
    assert b2.pinfo == batch.pinfo  # Check if the loaded batch has the same pinfo
    assert (
        b2.problem is None
    )  # turned None for serialization, check if the loaded batch has the same problem
    assert str(b2.parent_dir) == str(
        batch.parent_dir
    )  # Check if the loaded batch has the same parent_dir
    assert b2.seed == batch.seed  # Check if the loaded batch has the same seed
    assert b2.rep == batch.rep  # Check if the loaded batch has the same rep
    assert b2.out_dir == batch.out_dir  # Check if the loaded batch has the same out_dir
