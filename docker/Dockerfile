FROM python:3.12 AS jaix_base

WORKDIR /jaix
COPY setup.py setup.py
COPY jaix jaix
RUN pip install pip-tools && pip-compile --output-file requirements.txt setup.py

FROM python:3.12
COPY --from=jaix_base /jaix /jaix
WORKDIR /jaix

RUN pip install --no-cache-dir -r requirements.txt
ARG EXTRAS
RUN pip install --no-cache-dir -e .${EXTRAS}

WORKDIR /tests
