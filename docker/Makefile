#PYROSETTA_CONDA_VERSION=2021.49+release
TAG = v2021.49-dev61812

all: molfile_to_params rdkit_prody_pyrosetta
sif: ncaars_molfile_to_params.sif ncaars_rdkit_prody_pyrosetta.sif

molfile_to_params: molfile_to_params/Dockerfile molfile_to_params/rosetta
	docker build -t ncaars-molfile-to-params:$(TAG) $@

molfile_to_params/rosetta: molfile_to_params/clone_rosetta.sh
	./$^ $@ $(TAG)

rdkit_prody_pyrosetta: rdkit_prody_pyrosetta/Dockerfile rdkit_prody_pyrosetta/environment.yml
	docker build \
		-t ncaars-rdkit-prody-pyrosetta \
		--build-arg PYROSETTA_CONDA_CHANNEL=$(PYROSETTA_CONDA_CHANNEL) \
		$@

ncaars_molfile_to_params.sif: molfile_to_params/Dockerfile molfile_to_params/rosetta
	singularity build $@ docker-daemon://ncaars-molfile-to-params:$(TAG)

ncaars_rdkit_prody_pyrosetta.sif: rdkit_prody_pyrosetta/Dockerfile rdkit_prody_pyrosetta/environment.yml
	singularity build $@ docker-daemon://ncaars-rdkit-prody-pyrosetta:latest

.PHONY: all sif molfile_to_params rdkit_prody_pyrosetta
