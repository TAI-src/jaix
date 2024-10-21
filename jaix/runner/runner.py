from abc import abstractmethod
import logging

logger = logging.getLogger("DefaultLogger")


class Runner:
    @abstractmethod
    def run(self, env, opt_class, opt_config, *args, **kwargs):
        raise NotImplementedError()
