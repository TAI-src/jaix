from collections import defaultdict
from pathlib import Path


def find_data_files(
    folder: str | Path,
    file_type_pattern: str = "",
    problem_names: list[str] | None = None,
) -> dict[str, list[Path]]:
    path = Path(folder)
    files = list(path.rglob(file_type_pattern))

    # group files by problem names
    sorted_problem_names = sorted(problem_names or [])
    res_dict = defaultdict(list)
    for f in files:
        name_found = False
        for name in sorted_problem_names:
            if name in f.name:
                res_dict[name].append(f)
                name_found = True
                break
        if not name_found:
            res_dict["unknown"].append(f)
    return dict(res_dict)
