from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]

SOURCE = ROOT / "deps" / "reproblems"
DEST = ROOT / "jaix" / "env" / "utils" / "problem" / "re_problem"

ITEMS = [
    "reproblem_python_ver/reproblem.py",
    # "approximated_Pareto_fronts", # Don't need them for now
    "ideal_nadir_points",
]


def main():

    DEST.mkdir(parents=True, exist_ok=True)

    for item in ITEMS:
        source = SOURCE / item
        destination = DEST / source.name

        if source.is_dir():
            shutil.copytree(source, destination, dirs_exist_ok=True)
        elif source.is_file():
            shutil.copy2(source, destination)
        else:
            raise FileNotFoundError(source)


if __name__ == "__main__":
    main()
