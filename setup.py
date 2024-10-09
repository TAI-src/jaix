from distutils.core import setup

setup(
    name="jaix",
    version="0.1",
    packages=[
        "jaix",
    ],
    install_requires=["cma", "gymnasium", "coco-experiment", "regex"],
    license="GPL3",
    long_description="jaix",
)
