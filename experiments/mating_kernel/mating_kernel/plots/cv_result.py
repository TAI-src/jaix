import io
import json
from pathlib import Path
from typing import Self

import polars as pl
from cocoviz import Indicator, ProblemDescription
from cocoviz import Result as CVResult
from cocoviz import ResultSet as CVResultSet
from cocoviz._typing import FilePath, FrameLike, NumericVector

from mating_kernel.exp.utils.batch import Batch
from mating_kernel.problems.problem_info import ProblemInfo


class Result(CVResult):
    def __init__(
        self,
        data: FrameLike,
        batch: Batch,
        fevals_column: str = "a_fevals",
        algorithm: str | None = None,
    ):
        if algorithm is None:
            algorithm = "_".join(f"{k}={v}" for k, v in batch.settings.items())
        problem = Result.get_description(batch.pinfo)
        self.batch = batch
        self.pq_file_name = f"{batch.run_id}.parquet"
        super().__init__(algorithm, problem, data, fevals_column)

    def to_parquet(self, path: FilePath):
        """
        Serialise the result to parquet.

        Writes both the polars DataFrame *and* metadata
        (algorithm name + problem JSON) into parquet schema metadata.

        Parameters
        ----------
        path : path-like
            Output path. Existing file will be overwritten.

        Raises
        ------
        ImportError
            If the optional `pyarrow` package is not installed.
        """
        # NOTE: had to override to add batch metadata
        try:
            import pyarrow as pa
            import pyarrow.parquet as pq
        except ImportError as e:
            raise ImportError(
                "to_parquet requires the pyarrow package to be installed. Install it with: `pip install coco-visualize[parquet]`."
            ) from e

        tbl = self._data.to_arrow()
        batch_dict = self.batch.model_dump(mode="json")
        metadata = {
            "algorithm": self.algorithm,
            "problem": self.problem.to_json(),
            "batch": json.dumps(batch_dict),
        }
        schema = tbl.schema.with_metadata(
            {k.encode(): v.encode() for k, v in metadata.items()}
        )
        tbl = pa.Table.from_arrays(list(tbl.itercolumns()), schema=schema)
        pq.write_table(tbl, str(path))

    def write(self, dir: FilePath) -> Path:
        path = Path(dir) / self.pq_file_name
        path.parent.mkdir(parents=True, exist_ok=True)
        self.to_parquet(path)
        return path

    def at_indicator(self, indicator: Indicator | str, targets: NumericVector):
        # Since parent has a different signature, we need to override it since this method creates a new instance
        cv_result = CVResult(
            self.algorithm,
            self.problem,
            self._data,
            "__fevals",  # This is what the parent class renames the column to
        )
        at_indicator_result = cv_result.at_indicator(indicator, targets)
        new_result = Result(
            at_indicator_result._data, self.batch, "__fevals", self.algorithm
        )
        return new_result

    @classmethod
    def parent_from_parquet(cls, path: FilePath) -> pl.DataFrame:
        """
        Load the parent `Batch` from a parquet file.

        Parameters
        ----------
        path : path-like
            Parquet file to load. May also be an ``http://`` or
            ``https://`` URL, in which case the optional `requests`
            dependency is required (``pip install coco-visualize[http]``).

        Returns
        -------
        Batch

        Raises
        ------
        ImportError
            If the optional `pyarrow` package is not installed, or if
            `path` is an HTTP(S) URL but `requests` is not installed.
        """
        try:
            import pyarrow.parquet as pq
        except ImportError as e:
            raise ImportError(
                "parent_from_parquet requires the pyarrow package to be installed. Install it with: `pip install coco-visualize[parquet]`."
            ) from e

        if isinstance(path, str) and path.startswith(("http://", "https://")):
            try:
                import requests
            except ImportError as e:
                raise ImportError(
                    "Loading a Result from an http(s) URL requires the 'requests' package. "
                    "Install it with `pip install coco-visualize[http]`."
                ) from e

            response = requests.get(path)
            response.raise_for_status()
            tbl = pq.read_table(io.BytesIO(response.content))
        else:
            tbl = pq.read_table(path)

        return tbl

    @classmethod
    def from_parquet(cls, path: FilePath) -> Self:
        """
        Load a `Result` from parquet.

        Expects parquet schema metadata entries:
        - algorithm
        - problem (in ProblemDescription.to_json() form)

        Parameters
        ----------
        path : path-like
            Parquet file to load. May also be an ``http://`` or
            ``https://`` URL, in which case the optional `requests`
            dependency is required (``pip install coco-visualize[http]``).

        Returns
        -------
        Result

        Raises
        ------
        ImportError
            If the optional `pyarrow` package is not installed, or if
            `path` is an HTTP(S) URL but `requests` is not installed.
        """
        # NOTE: had to reproduce to change metadata
        tbl = cls.parent_from_parquet(path)

        assert hasattr(tbl.schema, "metadata") and tbl.schema.metadata is not None
        batch_data = tbl.schema.metadata[b"batch"].decode("utf8")
        batch = Batch.model_validate_json(batch_data)
        data = pl.from_arrow(tbl)

        return cls(data, batch, "__fevals")

    @staticmethod
    def get_description(pinfo: ProblemInfo) -> ProblemDescription:
        cocoviz_description = ProblemDescription(
            name=pinfo.uuid,
            instance=0,
            number_of_variables=pinfo.num_variables,
            number_of_objectives=pinfo.num_objectives,
        )

        return cocoviz_description


class ResultSet(CVResultSet):

    @classmethod
    def from_directory(cls, dir: FilePath) -> Self:
        """
        Load a `ResultSet` from a directory containing parquet files.
        Also looks into subdirectories.

        Parameters
        ----------
        dir : path-like
            Directory containing parquet files.

        Returns
        -------
        ResultSet

        Raises
        ------
        ImportError
            If the optional `pyarrow` package is not installed.
        """
        dir = Path(dir)
        if not dir.is_dir():
            raise ValueError(f"{dir} is not a directory.")

        results = []
        for pq_file in dir.rglob("*.parquet"):
            result = Result.from_parquet(pq_file)
            results.append(result)

        return cls(results)
