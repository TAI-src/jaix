from ttex.log import log_wandb_init, teardown_wandb_logger, log_wandb_artifact
from jaix.utils.experiment_context import ExperimentContext
from ttex.config import Config


class WandbSession:
    def __init__(self, config: Config):
        self.config = config

    def start(self, ctx: ExperimentContext):
        self.config_dict = {}
        try:
            # This is a workaround since not all
            # configs have to_dict method
            self.config_dict = self.config.to_dict()
        except AttributeError:
            pass
        except NotImplementedError:
            pass
        # actually init wandb
        self.run = log_wandb_init(
            run_config=self.config_dict, logger_name=ctx.get("wandb_logger_name")
        )
        ctx.set("exp_id", self.run.id)

    def end(self, ctx: ExperimentContext):
        if len(self.config_dict) == 0:
            return

        wandb_logger_name = ctx.get("wandb_logger_name")

        for artifact in ctx.get("artifact_paths", []):
            log_wandb_artifact(
                logger_name=wandb_logger_name,
                artifact_name=artifact.name,
                local_path=artifact.local_path,
                artifact_type=artifact.artifact_type,
                description=artifact.description,
            )

        teardown_wandb_logger(name=wandb_logger_name)
