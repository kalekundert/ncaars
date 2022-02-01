TAG = v2021.49-dev61812
MODE = release

all: molfile_to_params rdkit_prody_pyrosetta

molfile_to_params: molfile_to_params/Dockerfile molfile_to_params/rosetta
	docker build -t ncaars-molfile-to-params:$(TAG) $@

molfile_to_params/rosetta: molfile_to_params/clone_rosetta.sh
	./$^ $@ $(TAG)

rdkit_prody_pyrosetta: rdkit_prody_pyrosetta/Dockerfile rdkit_prody_pyrosetta/environment.yml
	docker build -t ncaars-rdkit-prody-pyrosetta --build-arg PYROSETTA_CONDA_CHANNEL=$(PYROSETTA_CONDA_CHANNEL) $@

.PHONY: molfile_to_params rdkit_prody_pyrosetta
