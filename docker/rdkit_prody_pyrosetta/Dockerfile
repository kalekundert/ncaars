FROM condaforge/mambaforge
ARG PYROSETTA_CONDA_CHANNEL
COPY environment.yml /environment.yml

# Install our packages in the base conda environment.  This is a bad idea in 
# general (because it means that our dependencies need to be compatible with 
# condas dependencies), but activating conda environments in docker containers 
# is not really possible, so this is the only simple alternative.
RUN \
  conda config --prepend channels $PYROSETTA_CONDA_CHANNEL && \
  conda env update -n base /environment.yml

