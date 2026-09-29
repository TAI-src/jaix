from distutils.core import setup
from setuptools import find_packages

__version__ = "0.0.1"


setup(
    name="tai_mk",
    version=__version__,
    packages=find_packages(),
    install_requires=["tai-jaix", "cobi"],
    license="GPL3",
    long_description="mating kernel",
    long_description_content_type="text/x-rst",
)
