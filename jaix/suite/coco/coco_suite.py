from jaix.env.singular.ec_env import (
    ECEnvironment,
    ECEnvironmentConfig,
)
from ttex.config import ConfigurableObjectFactory as COF, Config
from jaix.suite.suite import Suite, AggType
from jaix.suite.coco import COCOProblem
import cocoex as ex
import regex as re
import random as rnd
from typing import Optional, Dict


class COCOSuiteConfig(Config):
    def __init__(
        self,
        env_config: ECEnvironmentConfig,
        suite_name: str,
        suite_instance: str = "",
        suite_options: str = "",
        num_batches: int = 1,
        current_batch: int = 0,
        output_folder: str = "",
    ):
        self.suite_name = suite_name
        self.suite_instance = suite_instance
        self.suite_options = suite_options
        self.env_config = env_config
        self.num_batches = num_batches
        self.current_batch = current_batch
        self.output_folder = output_folder


class COCOSuite(Suite):
    config_class = COCOSuiteConfig  # type: ignore[assignment]

    def __init__(self, config: COCOSuiteConfig):
        super().__init__(config)
        if self.num_batches > 1:
            self.output_folder += "_batch%03dof%d" % (
                self.current_batch,
                self.num_batches,
            )

        self.suite = ex.Suite(self.suite_name, self.suite_instance, self.suite_options)

    def _get_agg_problem_dict(self, agg_type: AggType, seed: Optional[int] = None):
        if agg_type != AggType.INST:
            raise NotImplementedError()
        problems = {}  # type: Dict[int, Dict[int, ex.Problem]]
        for dim in self.suite.dimensions:
            problems[dim] = {}
            function_names = set(
                [
                    re.findall(r"_f[0-9]+_", name)[0]
                    for name in self.suite.ids("", f"d{dim:02d}", "")
                ]
            )
            functions = [re.findall(r"[0-9]+", name)[0] for name in function_names]
            for func in functions:
                instance_names = set(
                    [
                        re.findall(r"i[0-9]+", name)[0]
                        for name in self.suite.ids(f"f{func}", f"d{dim:02d}", "")
                    ]
                )
                instances = [
                    int(re.findall(r"[0-9]+", name)[0]) for name in instance_names
                ]
                # TODO: qd workaraound until observer moves to later
                shuff_inst = rnd.Random(seed).sample(instances, len(instances))
                problems[dim][func] = [
                    self.suite.get_problem_by_function_dimension_instance(
                        int(func), dim, inst
                    )
                    for inst in shuff_inst
                ]
        return problems

    def get_envs(self):
        observer = ex.Observer(self.suite_name, "result_folder: " + self.output_folder)
        for batch_counter, coco_func in enumerate(self.suite):
            # Only responsible for running part of the experiments
            if (
                batch_counter % self.num_batches
                != self.current_batch % self.num_batches
            ):
                continue
            func = COCOProblem(coco_func)
            func.observe_with(observer)
            env = COF.create(
                ECEnvironment,
                self.env_config,
                func=func,
                func_id=func.problem.id_function,
                inst=func.problem.id_instance,
            )
            yield env

    def get_agg_envs(
        self, agg_type: AggType = AggType.NONE, seed: Optional[int] = None
    ):
        # Currently, this only makes sense for single batches
        assert self.num_batches == 1
        assert self.current_batch == 0
        problems_dict = self._get_agg_problem_dict(agg_type, seed)
        for dim, funcs in problems_dict.items():
            for fun, coco_problems in funcs.items():
                observers = [
                    ex.Observer(
                        self.suite_name,
                        f"result_folder: {self.output_folder}_s{i}/{prob.id}",
                    )
                    for i, prob in zip(range(len(coco_problems)), coco_problems)
                ]
                funcs = [COCOProblem(prob) for prob in coco_problems]
                for obs, func in zip(observers, funcs):
                    func.observe_with(obs)
                envs = [
                    COF.create(ECEnvironment, self.env_config, func) for func in funcs
                ]
                yield envs
                # TODO: combine logs for post-processing (or add it to the post-processing)
