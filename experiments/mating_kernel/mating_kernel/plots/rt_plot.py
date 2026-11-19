from pathlib import Path

import matplotlib.pyplot as plt
from cocoviz import (
    Indicator,
    ResultSet,
    rtpplot,
)


def plot_rt(
    results: ResultSet,
    indicator_specs: list[tuple[str, bool]],
    output_dir: str | Path,
    file_prefix: str = "",
) -> list[str]:
    file_names = []
    for problem, by_prob in results.by_problem_name():
        for ind_col, lib in indicator_specs:
            ind = Indicator(ind_col, larger_is_better=lib)
ax = rtpplot(by_prob, ind)
            # save the plot to a file
            file_name = (
                f"{file_prefix}_rt_{problem}_{ind_col}.pdf"
                if file_prefix
                else f"rt_{problem}_{ind_col}.pdf"
            )
            fig = ax.figure
            file_path = Path(output_dir) / Path(file_name)
            fig.savefig(file_path, dpi=300)
            plt.close(fig)
            file_names.append(str(file_path))
    return file_names
