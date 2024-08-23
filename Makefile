.PHONY: conda update_submodules docker singularity docker-upload

ORG := "logsdonlab"
PROJECT_NAME := "cenmap"
VERSION := "latest"
REPO := "${ORG}/${PROJECT_NAME}"
TAG_NAME := "${REPO}:${VERSION}"

conda:
	conda env create -f env.yaml --name "${TAG}"
	conda activate "${TAG}"

update_submodules:
	git submodule init
	git submodule update --remote --recursive

docker:
	snakemake -p -c 1 --sdm apptainer --configfile test/config/config.yaml --containerize | \
	sudo docker build -t "${TAG_NAME}" -f - .

singularity:
	$(MAKE) docker
	sudo apptainer build "${PROJECT_NAME}.sif" "docker-daemon://${TAG_NAME}"

docker-upload:
	$(MAKE) docker
	sudo docker image push "${TAG_NAME}"
