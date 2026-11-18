import hashlib
import json
import logging
import pickle
from pathlib import Path
from typing import Self

from pydantic import BaseModel, ConfigDict

from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.problems.problem_info import ProblemInfo

logger = logging.getLogger(__name__)


class Batch(BaseModel):
    pid: str
    sid: int
    rep: int | None = None
    seed: int | None = None
    pinfo: ProblemInfo
    problem: StaticProblem | None = None
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
    def batch_file(self) -> Path:
        """
        The file to save batch info to
        """
        return self.out_dir / f"batch_{self.name}.json"

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
        # TODO: update
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

    def exists(self) -> bool:
        """
        Check if a folder with the same run id already exists in the given path
        """
        # Check for batch file since they are the last thing to be written
        batch_exists = self.batch_file.exists()
        if not batch_exists:
            logger.debug(
                f"Batch file {self.batch_file} does not exist for batch {self.name}."
            )
            # For backward compatibility, also check for pickle file
            if self.out_dir.exists():
                pkl_file = self.batch_file.with_suffix(".pkl")
                logger.debug(
                    f"Checking if pickle file {pkl_file} exists for batch {self.name}."
                )
                batch_exists = Path(pkl_file).exists()
        logger.debug(f"Batch {self.name} exists: {batch_exists}")
        return batch_exists

    def write(self) -> Path:
        data = self.model_dump(
            mode="json",
            exclude={"problem"},
        )

        self.batch_file.write_text(json.dumps(data, indent=2))
        logger.debug(f"Batch {self.name} written to {self.batch_file}")

        return self.batch_file

    @classmethod
    def load(cls, batch_file: str | Path) -> Self:
        if Path(batch_file).suffix == ".pkl":
            with open(batch_file, "rb") as f:
                data = pickle.load(f)
            # Ensure problem is set to None when loading from pickle
            # This is to ensure consistency with the JSON loading, which does not include the problem object.
            setattr(data, "problem", None)  # noqa: B010
        elif Path(batch_file).suffix == ".json":
            data = json.loads(Path(batch_file).read_text())
        else:
            raise ValueError(f"Unsupported file format: {batch_file}")
        return cls.model_validate(data)


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
