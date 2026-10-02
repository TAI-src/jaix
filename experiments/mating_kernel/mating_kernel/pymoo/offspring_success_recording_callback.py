import copy

import pandas as pd

from jaix.env.utils.problem.static_problem import StaticProblem
from mating_kernel.pymoo.parser.population_parser import PopulationParser
from mating_kernel.pymoo.parser.reproduction_parser import ReproductionParser
from mating_kernel.pymoo.problem_wrapper import PymooProblemWrapper
from mating_kernel.pymoo.recording_callback import RecordingCallback


class OffspringSuccessRecordingCallback(RecordingCallback):
    def __init__(self, problem: StaticProblem):
        if not hasattr(problem, "static_problem"):
            # Wrap in pymoo problem wrapper if not already wrapped
            problem = PymooProblemWrapper(problem)
        if not hasattr(problem.static_problem, "get_archive_stats"):
            raise ValueError(
                "The static_problem must have a get_archive_stats method, i.e. be tracked using MOTrackingMixin."
            )
        assert hasattr(
            problem.static_problem, "ideal_point"
        ), "The static_problem must have an ideal_point attribute."
        parsers = [
            ReproductionParser(ideal=problem.static_problem.ideal_point),
            PopulationParser(ideal=problem.static_problem.ideal_point),
        ]
        super().__init__(recording_parsers=parsers)
        self._dirty = True  # Mark as dirty to indicate that new data has been added
        self._record_stats_df = (
            pd.DataFrame()
        )  # Initialize an empty DataFrame to store the full records

    def get_stat_gen(self, gen: int, attr: str) -> dict:
        if gen < 0 or gen >= len(self.data["record_stats"]):
            raise IndexError(f"Generation {gen} is out of bounds.")
        dat = self.data["record_stats"][gen].get(attr, [])
        assert (
            len(dat) == 1
        ), f"Expected a single record for {attr} at generation {gen}, but got {len(dat)}."
        return dat[0]

    def _get_full_record_dicts(self) -> list[dict]:
        full_records = []
        data = self.data["record_stats"]

        for gen in range(len(data)):
            if gen == 0:
                continue  # Skip the first generation as there are no offspring yet
            b_pop = self.get_stat_gen(gen - 1, "PopulationParser")
            a_pop = self.get_stat_gen(gen, "PopulationParser")
            b_archive = self.get_stat_gen(gen - 1, "archive_stats")
            a_archive = self.get_stat_gen(gen, "archive_stats")
            for fam in data[gen].get("ReproductionParser", []):
                fam_record = copy.deepcopy(fam)
                fam_record.update(
                    {f"b_{k}": v for k, v in b_pop.items()} if b_pop else {}
                )
                fam_record.update(
                    {f"a_{k}": v for k, v in a_pop.items()} if a_pop else {}
                )
                fam_record.update(
                    {f"b_{k}": v for k, v in b_archive.items()} if b_archive else {}
                )
                fam_record.update(
                    {f"a_{k}": v for k, v in a_archive.items()} if a_archive else {}
                )
                full_records.append(fam_record)
        return full_records

    def notify(self, algorithm):
        super().notify(algorithm)
        self._dirty = True  # Mark as dirty to indicate that new data has been added

    @property
    def record_stats(self) -> pd.DataFrame:
        if self._dirty:
            full_records = self._get_full_record_dicts()
            df = pd.DataFrame(full_records)
            # Remove duplicate child records based on 'o_X' and 'o_F' columns, keeping the first occurrence
            if df.empty:
                self._record_stats_df = df
                self._dirty = False
                return self._record_stats_df
            duplicate_keys = pd.DataFrame(
                {"o_X": df["o_X"].map(tuple), "o_F": df["o_F"].map(tuple)}
            )
            self._record_stats_df = df[
                ~duplicate_keys.duplicated(keep="first")
            ].reset_index(drop=True)
            self._dirty = False  # Reset the dirty flag after updating the DataFrame
        return self._record_stats_df
