from distutils.core import setup
from setuptools import find_packages

__version__ = "0.1.1.19"

setup(
    name="tai_jaix",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        "cma",
        "gymnasium",
        "tai-ttex",
        "ase",
        "kimpy",
        "requests",
        "regex",
        "coco-experiment",
    ],
    # extras_require={
    #    'ase': ['ase', 'kimpy', 'requests'],
    #    'tabrepo': ['tabrepo', 'regex'],
    #    'coco': ['coco-experiment', 'regex'],
    # },
    license="GPL3",
    long_description="jaix",
    long_description_content_type="text/x-rst",
)
