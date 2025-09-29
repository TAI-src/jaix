from jaix.runner.ask_tell.at_strategy import ATStrategy
from ttex.config import Config, ConfigurableObject
from gymnasium import Env, spaces
import itertools


class EnumerateATStratConfig(Config):
    def __init__(self, ask_size: int):
        self.ask_size = ask_size


class EnumerateATStrat(ConfigurableObject, ATStrategy):
    config_class = EnumerateATStratConfig

    def __init__(self, config: EnumerateATStratConfig, env: Env, *args, **kwargs):
        ConfigurableObject.__init__(self, config)
        # For now only for specific spaces
        assert isinstance(env.action_space, spaces.MultiBinary)
        self.xstart = [env.action_space.sample() for _ in range(self.ask_size)]
        ATStrategy.__init__(self, self.xstart)
        # TODO: figure out No overload variant of "product" matches argument types "list[int]", "tuple[int, ...]"  [call-overload]
        self.options = list(itertools.product([0, 1], repeat=env.action_space.n))  # type: ignore

    def initialize(self):
        self.current = 0

    @property
    def name(self):
        return "Enumerate Search"

    def ask(self, env, **kwargs):
        # Check if we are done
        assert self.current < len(self.options)
        # Sample the next ask_size options
        end = min(self.current + self.ask_size, len(self.options))
        self.xcurrent = self.options[self.current : end]
        self.current += self.ask_size
        return self.xcurrent

    def tell(self, **kwargs):
        pass

    def stop(self):
        if self.current >= len(self.options):
            return {"enumeration_done": True}
        return {}
