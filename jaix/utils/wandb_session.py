from jaix.experiment import ExperimentConfig
from ttex.log import log_wandb_init, teardown_wandb_logger, log_wandb_artifact
from jaix.utils.experiment_context import ExperimentContext


class WandbSession:

    def __init__(self, config: ExperimentConfig):
        self.config = config

    def start(self, ctx: ExperimentContext):
        config_dict = self.config.to_dict()
        # actually init wandb
        self.run = log_wandb_init(
            run_config=config_dict, logger_name=ctx.get("wandb_logger_name")
        )
        ctx.set("exp_id", self.run.id)

    def end(self, ctx: ExperimentContext):
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
