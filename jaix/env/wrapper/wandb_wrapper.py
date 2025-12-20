from jaix.env.wrapper.value_track_wrapper import ValueTrackWrapper
import gymnasium as gym
from ttex.config import ConfigurableObject, Config
import logging
from typing import Optional, Dict, List
from ttex.log import setup_wandb_logger
import jaix.utils.globals as globals


class WandbWrapperConfig(Config):
    def __init__(
        self,
        logger_name: Optional[str] = None,
        custom_metrics: Optional[Dict] = None,
        snapshot: bool = True,
        snapshot_sensitive_keys: Optional[List[str]] = None,
        project: Optional[str] = None,
        group: Optional[str] = None,
        passthrough: bool = True,
        state_eval: str = "obs0",  # Which value should be logged
        is_min: bool = True,  # Whether lower is better for state_eval
    ):
        self.passthrough = passthrough
        self.custom_metrics = custom_metrics
        self.snapshot = snapshot
        self.snapshot_sensitive_keys = snapshot_sensitive_keys
        self.project = project
        self.group = group
        self.state_eval = state_eval
        self.is_min = is_min
        self.logger_name = (
            logger_name
            if (logger_name is not None)  # Avoid using root logger
            else globals.WANDB_LOGGER_NAME
        )
        if self.logger_name == globals.LOGGER_NAME:
            raise ValueError(
                "WandbWrapperConfig: logger_name cannot be the root logger name."
            )
        globals.WANDB_LOGGER_NAME = self.logger_name

    def _setup(self):  # Setup wandb logger
        logger = setup_wandb_logger(
            name=self.logger_name,
            custom_metrics=self.custom_metrics,
            snapshot=self.snapshot,
            snapshot_sensitive_keys=self.snapshot_sensitive_keys,
            project=self.project,
            group=self.group,
        )
        assert logger is not None
        return True


class WandbWrapper(ConfigurableObject, ValueTrackWrapper):
    """
    A wrapper that logs environment interactions to wandb.
    It logs rewards, resets, steps, and other relevant information.
    It can also log custom metrics and environment info on close.
    """

    config_class = WandbWrapperConfig

    def __init__(self, config: WandbWrapperConfig, env: gym.Env):
        ConfigurableObject.__init__(self, config)
        ValueTrackWrapper.__init__(
            self,
            env,
            passthrough=self.passthrough,
            state_eval=self.state_eval,
            is_min=self.is_min,
        )
        self.logger = logging.getLogger(self.logger_name)
        self.log_resets = 0
        self.log_env_steps = 0
        self.log_renv_steps = 0
        self.last_info = {}  # type: dict

    def reset(self, **kwargs):
        obs, info = self.env.reset(**kwargs)
        self.log_resets += 1
        self.log_renv_steps = 0
        self.last_info = info
        return obs, info

    def step(self, action):
        (
            obs,
            r,
            term,
            trunc,
            info,
        ) = self.env.step(action)
        self.last_info = info
        self.log_env_steps += 1
        self.log_renv_steps += 1

        # Get the value for tracking and update internal state
        val = self.get_val(obs, r, info, self.state_eval)
        self.update_vals(val)

        # TODO: Option for log frequency

        env_step = info["env_step"] if "env_step" in info else self.log_env_steps
        # Log per reset
        info_dict = {
            f"env/resets/{self.env.unwrapped}": self.log_resets,
            # f"restarts/r/{self.dim}/{self.env}/{self.log_resets}": r.item(),
            "env/step": env_step,
            "env/log_step": self.log_env_steps,
            # "restarts/step": self.log_renv_steps,
        }
        if r is not None:
            info_dict[f"env/r/{str(self.env.unwrapped)}"] = float(r.item())
        if f"raw_{self.state_eval}" in info:
            # The original value might have been overwritten, so get it from info if available
            raw_val = float(info[f"raw_{self.state_eval}"])
            best_raw_val = (
                float(info[f"best_raw_{self.state_eval}"])
                if f"best_raw_{self.state_eval}" in info
                else raw_val
            )
        else:
            # Get the raw value from the current step
            raw_val = float(val)
            best_raw_val = float(self.best_val if self.best_val is not None else val)

        info_dict[f"env/raw_{self.state_eval}/{str(self.env.unwrapped)}"] = raw_val
        info_dict[
            f"env/best_raw_{self.state_eval}/{str(self.env.unwrapped)}"
        ] = best_raw_val

        if term:
            info_dict[f"env/term/{str(self.env.unwrapped)}"] = float(
                self.log_renv_steps
            )

        self.logger.info(info_dict)
        # TODO: Figure out what info would be helpful from all the sub-wrappers etc
        return obs, r, term, trunc, info

    def close(self):
        env_step = (
            self.last_info["env_step"]
            if "env_step" in self.last_info
            else self.log_env_steps
        )
        closing_info = {
            "env/step": env_step,
            "env/log_step": self.log_env_steps,
        }
        for key, value in self.last_info.items():
            # TODO: for now not nested
            # TODO: Test for loggable values for wandb. if not, issue warning. Probably do this in ttex
            if isinstance(value, dict):
                continue
            closing_info[f"env/close/{str(self.env.unwrapped)}/{key}"] = value

        self.logger.info(closing_info)
        self.env.close()
