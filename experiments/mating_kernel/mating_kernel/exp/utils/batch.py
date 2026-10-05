import hashlib
import json
from pathlib import Path

from pydantic import BaseModel, ConfigDict

from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.problems.problem_info import ProblemInfo


class Batch(BaseModel):
    pid: str
    sid: int
    rep: int | None = None
    seed: int | None = None
    problem: StaticProblem
    pinfo: ProblemInfo
    parent_dir: str | Path = Path(".")
    model_config = ConfigDict(
        extra="allow", arbitrary_types_allowed=True
    )  # Allow extra fields in the model

    @property
    def name(self) -> str:
        """
        Generate a unique name for the batch based on its identifiers.
        """
        if self.rep is None:
            return f"p{self.pid}_s{self.sid}"
        return f"p{self.pid}_s{self.sid}_r{self.rep}"

    @property
    def out_dir(self) -> Path:
        """
        Generate the output directory path for the batch based on its identifiers.
        """
        out_path = Path(self.parent_dir) / self.run_id
        out_path.mkdir(parents=True, exist_ok=True)
        return out_path

    @property
    def settings(self) -> dict:
        return self.model_extra or {}

    @property
    def experiment_id(self) -> str:
        """
        Generate a unique name for the batch based on its identifiers.
        """
        hashed_settings = hash_dict({**self.settings})
        return f"{self.pid}_{hashed_settings}"

    @property
    def run_id(self) -> str:
        """
        Generate a unique name for the batch based on its identifiers.
        """
        if self.seed is None:
            raise ValueError("Seed must be set to generate run_id.")
        return f"{self.experiment_id}_{self.seed}"

    @staticmethod
    def parse_name(batch_name: str) -> tuple[str, int]:
        """
        Extract the batch identifiers from a batch name.
        """
        parts = batch_name.rsplit("_", 1)
        if len(parts) != 2:
            raise ValueError(f"Invalid batch name format: {batch_name}")
        settings_hash, seed_str = parts
        return settings_hash, int(seed_str)

    def get_num_runs(self, dir: str | Path | None = None) -> int:
        """
        Check in the given path how many folders with the experiment id exist
        """
        dir = Path(dir) if dir is not None else Path(self.parent_dir)
        if not dir.exists():
            return 0
        return len(
            [
                d
                for d in dir.iterdir()
                if d.is_dir() and d.name.startswith(self.experiment_id)
            ]
        )

    def get_run_seeds(self, dir: str | Path | None = None) -> list[int]:
        """
        Get the seeds of all runs that have been executed in the given path
        """
        dir = Path(dir) if dir is not None else Path(self.parent_dir)
        if not dir.exists():
            return []
        seeds = []
        for d in dir.iterdir():
            if (
                d.is_dir()
                and d.name.startswith(self.experiment_id)
                and any(f.name.startswith("batch_") for f in d.iterdir())
            ):
                # Check if the batch file is there, since it is the last thing to be written, if it is there, the run is complete
                # batch files contain "batch_" in the name
                # FIXME: This is a hacky way to check if the run is complete, but it works for now.
                # since it hardcodes the batch file name, it is not very robust.
                _, seed = self.parse_name(d.name)
                seeds.append(seed)
        return seeds

    def exists(self, dir: str | Path | None = None) -> bool:
        """
        Check if a folder with the same run id already exists in the given path
        """
        dir = Path(dir) if dir is not None else Path(self.parent_dir)
        if not dir.exists():
            return False
        # check if self.out_dir is empty (it will always exist as it is created in the out_dir property)
        run_dir = dir / self.run_id
        if not run_dir.exists():
            return False
        # Check for batch file since they are the last thing to be written
        # FIXME: This is a hacky way to check if the run is complete, but it works for now.
        if any(f.name.startswith("batch_") for f in run_dir.iterdir()):
            return True
        else:
            return False


def hash_dict(d: dict, length: int = 12) -> str:
    """
    Hash a dictionary to create a unique identifier.
    The dictionary is first converted to a JSON string with sorted keys to ensure consistent hashing.
    """
    dict_str = json.dumps(d, sort_keys=True)
    hash_str = hashlib.sha256(dict_str.encode()).hexdigest()
    if length < len(hash_str):
        return hash_str[:length]  # Return the first `length` characters of the hash
    else:
        return (
            hash_str  # Return the full hash if length is greater than the hash length
        )
