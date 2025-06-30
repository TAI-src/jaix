FROM taisrc/jaix:tabrepo

RUN sudo apt install -y pkgconf
RUN sudo add-apt-repository -y ppa:openkim/latest && sudo apt-get update && sudo apt-get install -y libkim-api-dev openkim-models

# Preload python packages required
RUN pip install pytest black mypy pylint

WORKDIR /jaix
COPY setup.py setup.py
COPY jaix jaix

USER root
RUN pip install -e .

WORKDIR /tests
