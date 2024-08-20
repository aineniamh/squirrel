FROM condaforge/mambaforge:latest AS conda

COPY environment.yml .

RUN /opt/conda/bin/mamba env create -f /environment.yml

ENV PATH=/opt/conda/envs/squirrel/bin:$PATH

RUN pip install git+https://github.com/aineniamh/squirrel.git

ENV MPLCONFIGDIR="."
ENV XDG_CACHE_HOME="./.cache"
ENV TMPDIR="./.tmp"

CMD ["/bin/bash"]
