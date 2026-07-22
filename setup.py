from distutils.core import setup
from setuptools import find_packages

__version__ = "0.1.1.69"

archive_packages = [
    "seaborn",
    "pandas",
    "matplotlib",
    "pymoo",
    "scikit-learn",
]

setup(
    name="tai_jaix",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        "numpy",
        "cma",
        "gymnasium",
        "tai-ttex",
        "cocopp",
    ]
    + archive_packages,
    extras_require={
        "ase": ["ase", "kimpy", "requests"],
        "tabrepo": ["tabrepo", "regex"],
        "coco": ["coco-experiment", "regex"],
        "test": ["pytest", "mypy", "pylint", "black"],
    },
    license="GPL3",
    long_description="jaix",
    long_description_content_type="text/x-rst",
)
