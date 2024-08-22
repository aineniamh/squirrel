FROM condaforge/mambaforge:latest AS conda

COPY . ./squirrel/

RUN /opt/conda/bin/mamba env create -f ./squirrel/environment.yml

ENV PATH=/opt/conda/envs/squirrel/bin:$PATH

RUN cd squirrel && pip install .

ENV MPLCONFIGDIR="."

CMD ["/bin/bash"]
