FROM taisrc/jaix:tabrepo

# Preload python packages required
RUN pip install pytest black mypy pylint

WORKDIR /jaix
COPY setup.py setup.py
COPY jaix jaix

RUN pip install -e .

WORKDIR /tests
