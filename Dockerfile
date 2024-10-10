# start with an image with conda installed
FROM condaforge/mambaforge AS compile-image

# set working directory
WORKDIR /data

# check for updates
RUN apt-get update -y && \
  apt-get upgrade -y && \
  apt install build-essential -y --no-install-recommends && \
  apt-get clean && apt-get autoclean

# copy in squirrel
COPY . /data/squirrel/
RUN cd /data/squirrel && \
  mamba install conda -n base -c conda-forge -c defaults -c bioconda && \
  mamba env create -f /data/squirrel/environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "squirrel", "/bin/bash", "-c"]
RUN mamba install -c conda-forge -n squirrel python=3.10 conda-pack && \
  cd /data/squirrel && \
  pip install .

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda list && \
  conda-pack -n squirrel -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar && /venv/bin/conda-unpack

SHELL ["/bin/bash", "-c"]

RUN conda clean --all &&\
  conda remove --name squirrel --all

# build squirrel
WORKDIR /data/squirrel
RUN source /venv/bin/activate && pip install --user --no-cache-dir . 

# build image
FROM debian:bullseye-slim AS runtime-image

COPY --from=compile-image /root/.local /root/.local
ENV PATH=/root/.local/bin:$PATH

# Copy /venv from the previous stage:
COPY --from=compile-image /venv /venv

# create directory to mount the basecalled directory and output directory
RUN mkdir -p /data/run_data/basecalled && mkdir -p /data/run_data/output

# check for updates
RUN apt-get update -y && \
  apt-get upgrade -y && \
  apt install build-essential -y --no-install-recommends && \
  apt install -y procps && \
  apt-get clean && apt-get autoclean

# to allow streamed log output
ENV PYTHONUNBUFFERED=1
ENV PATH=/venv/bin:$PATH
ENV MPLCONFIGDIR="."

CMD ["/bin/bash"]
