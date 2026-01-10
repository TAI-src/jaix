from jaix.utils.experiment_context import ExperimentContext
from ttex.config.config import Config
from ttex.log import initiate_logger, get_logging_config
from typing import Optional, Dict


class LoggingConfig(Config):
    def __init__(
        self,
        log_level: int = 30,
        logger_name: Optional[str] = None,
        disable_existing: Optional[bool] = False,
        dict_config: Optional[Dict] = None,
    ):
        self.log_level = log_level
        self.logger_name = logger_name
        self.disable_existing = disable_existing
        self.dict_config = dict_config

    def _setup(self, ctx: ExperimentContext):
        self.logger_name = (
            self.logger_name if self.logger_name is not None else ctx.get("logger_name")
        )
        self.dict_config = (
            self.dict_config
            if self.dict_config
            else get_logging_config(self.logger_name, False)
        )

        initiate_logger(
            log_level=self.log_level,
            logger_name=self.logger_name,
            disable_existing=self.disable_existing,
            logging_config=self.dict_config,
        )
        return True
