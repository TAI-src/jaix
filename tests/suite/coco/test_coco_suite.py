from jaix.suite.coco import COCOSuiteConfig, COCOSuite, COCOProblem
import pytest
from jaix.env.singular.ec_env import ECEnvironmentConfig
from ttex.config import (
    ConfigurableObjectFactory as COF,
)
import os
from jaix.suite.suite import AggType
from .test_coco_problem import skip_remaining_tests


@pytest.fixture(scope="function")
def env_config(skip_remaining_tests):
    config = ECEnvironmentConfig(budget_multiplier=1)
    return config


def test_init(env_config):
    config = COCOSuiteConfig(env_config, "bbob")
    suite = COF.create(COCOSuite, config)

    assert suite.suite_name == "bbob"
    assert suite.env_config.budget_multiplier == 1


def test_get_envs(env_config):
    config = COCOSuiteConfig(
        env_config,
        "bbob",
        suite_instance="instances: 1",
        suite_options="function_indices: 1,2 dimensions: 2,3",
        num_batches=2,
        current_batch=0,
        output_folder="test_run",
    )
    assert config.output_folder == "test_run"
    suite = COF.create(COCOSuite, config)
    assert suite.output_folder == "test_run_batch000of2"
    counter = 0
    for env in suite.get_envs():
        assert isinstance(env.func, COCOProblem)
        assert not env.stop()
        assert env.func.id_function == 1
        counter += 1

        # Test observation
        assert env.func.is_observed
        assert env.func.evaluations == 0
        rec_file = env.close()
        if rec_file is not None:
            os.remove(rec_file)
    assert counter == 2


def test_get_envs_agg(env_config):
    config = COCOSuiteConfig(
        env_config,
        "bbob",
        suite_instance="instances: 1,2,3",
        suite_options="function_indices: 1,2 dimensions: 2,3",
        num_batches=1,
        current_batch=0,
        output_folder="test_run",
    )
    suite = COF.create(COCOSuite, config)
    for envs in suite.get_agg_envs(agg_type=AggType.INST, seed=1):
        assert len(envs) == 3
        for env in envs:
            assert isinstance(env.func, COCOProblem)
            assert env.func.is_observed

            # env.reset()
            env.step(env.action_space.sample())
        rec_files = [env.close() for env in envs]
        for rec_file in rec_files:
            if rec_file is not None:
                os.remove(rec_file)


def test_get_problem_dict(env_config):
    import cocoex as ex

    suite = ex.Suite("bbob", "", "")
    # TODO: check how this works for bi-objective

    config = COCOSuiteConfig(env_config, "bbob", output_folder="test_run")
    coco_suite = COF.create(COCOSuite, config)
    with pytest.raises(NotImplementedError):
        problems = coco_suite._get_agg_problem_dict(AggType.NONE, seed=1337)
    problems = coco_suite._get_agg_problem_dict(AggType.INST, seed=1337)
    assert problems.keys() == set(suite.dimensions)
    # check all functions are in each dimensions
    assert all([len(func_per_dim.keys()) == 24 for func_per_dim in problems.values()])
    for dim, fun in problems.items():
        for fun_id, coco_problems in fun.items():
            assert len(coco_problems) == 15
            for prob in coco_problems:
                f, d, i = prob.id_triple
                assert f == int(fun_id)
                assert d == dim
    suite.free()
