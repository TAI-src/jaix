from abc import abstractmethod
from enum import Enum
from typing import Optional
from jaix.env.wrapper import WrappedEnvFactory as WEF, ClosingWrapper, OnlineWrapper


class AggType(Enum):
    NONE = 0
    INST = 1


class Suite:

    @abstractmethod
    def get_envs_generator(self,
                 agg_type: AggType = AggType.NONE,
                 seed: Optional[int] = None):
        raise NotImplementedError()
    
    def get_envs(self,
                 agg_type: AggType = AggType.NONE,
                 seed: Optional[int] = None):
        # TODO: Potentially add a non-generator version later
        for env in self.get_envs_generator(agg_type = agg_type, seed=seed):
            wrappers = [
                (ClosingWrapper, {}),
                (OnlineWrapper, {"online": True}),
            ]
            # TODO: And to force reset before use
            wrapped_env = WEF.wrap(env, wrappers)
            yield wrapped_env
            assert wrapped_env.closed
