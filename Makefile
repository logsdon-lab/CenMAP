.PHONY: conda conda_all update_submodules docker_local singularity_local docker_local_upload dockerfile

ORG := "logsdonlab"
PROJECT_NAME := "cenmap"
VERSION := "latest"
REPO := "${ORG}/${PROJECT_NAME}"
TAG_NAME := "${REPO}:${VERSION}"

conda:
	conda env create -f env.yaml --name "${TAG}"
	conda activate "${TAG}"

# Make an environment with all dependencies.
# Need to manually clean.
conda_all:
	bash workflow/scripts/make_env_all.sh

update_submodules:
	git submodule init
	git submodule update --init --recursive --remote

dockerfile:
	snakemake -p -c 1 --sdm apptainer --configfile test/config/config.yaml --containerize > Dockerfile

docker_local:
	snakemake -p -c 1 --sdm apptainer --configfile test/config/config.yaml --containerize | \
	sudo docker build -t "${TAG_NAME}" -f - .

singularity_local:
	$(MAKE) docker_local
	sudo apptainer build "${PROJECT_NAME}.sif" "docker-daemon://${TAG_NAME}"

docker_local_upload:
	$(MAKE) docker_local
	sudo docker image push "${TAG_NAME}"
