TAG = v2021.49-dev61812
MODE = release

molfile_to_params: molfile_to_params/Dockerfile molfile_to_params/rosetta
	docker build -t ncaars-molfile-to-params:$(TAG) $@

molfile_to_params/rosetta: molfile_to_params/clone_rosetta.sh
	./$^ $@ $(TAG)

.PHONY: molfile_to_params
