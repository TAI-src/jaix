from mating_kernel.pymoo.parser.recording_parser import RecordingParser
import numpy as np
import pandas as pd


class PopulationParser(RecordingParser):
    record_args = ["survival"]
    record_attributes = []
    record_retrieval = ["survival"]

    def __init__(self, ideal: np.ndarray | None = None):
        self.ideal = ideal

    def parse(self, data: dict[str, list]) -> list[dict]:
        # Get population after final selection step per generation (should only be one entry per generation)
        population = data["survival"][-1]["output"]
        ind_data = [ind.data for ind in population]
        if self.ideal is not None:
            for data, ind in zip(ind_data, population):
                data["dist_to_ideal"] = np.linalg.norm(ind.F - self.ideal)
        # Aggregate ind data
        df = pd.DataFrame(ind_data)
        stats = df.replace([np.inf, -np.inf], np.nan).describe().to_dict()
        # flatten stats dictionary
        stats = {
            f"{k}_{stat}": v
            for k, stat_dict in stats.items()
            for stat, v in stat_dict.items()
        }
        return [stats]
