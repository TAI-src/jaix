from ttex.config import Config, ConfigurableObject

from jaix.env.utils.problem.static_problem import StaticProblem


class ExpDesignProblemConfig(Config):
    def __init__(
        self,
    ):
        Config.__init__(self)


class ExpDesignProblem(ConfigurableObject, StaticProblem):
    config_class = ExpDesignProblemConfig

    def __init__(self, config: ExpDesignProblemConfig, inst: int):
        ConfigurableObject.__init__(self, config)
        StaticProblem.__init__(self, dimension=0, num_objectives=0, precision=0.0)
        self.inst = inst

    def _eval(self, x):
        raise NotImplementedError("This is a placeholder for the actual problem.")
