from distutils.core import setup
from setuptools import find_packages

__version__ = "0.1.1.7"

setup(
    name="tai_jaix",
    version=__version__,
    packages=find_packages(),
    install_requires=["cma", "gymnasium", "coco-experiment", "regex", "tai-ttex"],
    license="GPL3",
    long_description="jaix",
    long_description_content_type="text/x-rst",
)
