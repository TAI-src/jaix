from pydantic import BaseModel, ConfigDict
import hashlib
import json
from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.problems.problem_info import ProblemInfo
from pathlib import Path


class Batch(BaseModel):
    pid: str
    sid: int
    rep: int | None = None
    seed: int | None = None
    problem: StaticProblem
    pinfo: ProblemInfo
    parent_dir: str | Path = Path(".")
    model_config = ConfigDict(extra="allow")  # Allow extra fields in the model

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
        out_path = Path(self.parent_dir) / self.name
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
        hashed_settings = self.hash_dict({**self.settings})
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
            if d.is_dir() and d.name.startswith(self.experiment_id):
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
        return any(d.is_dir() and d.name == self.run_id for d in dir.iterdir())


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
