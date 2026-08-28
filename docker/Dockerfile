FROM python:3.12 AS jaix_base

WORKDIR /jaix
COPY setup.py setup.py
COPY jaix jaix
RUN pip install pip-tools && pip-compile --output-file requirements.txt setup.py

RUN pip install --prefix=/opt/jaix --no-cache-dir -r requirements.txt
RUN pip install --prefix=/opt/jaix --no-cache-dir -e .

FROM python:3.12
COPY --from=jaix_base /jaix /jaix
COPY --from=jaix_base /opt/jaix/ /usr/local
WORKDIR /jaix

ARG EXTRAS
RUN pip install --no-cache-dir -e .${EXTRAS}

WORKDIR /tests
