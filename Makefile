.PHONY: conda dev update_submodules

bin=$(CURDIR)/venv/bin

conda:
	conda env create -f env.yaml --name cenmap
	conda activate cenmap

dev:
	python3 -m venv venv
	$(bin)/python3 -m pip install -r requirements-dev.txt

update_submodules:
	git submodule init
	git submodule update --remote --recursive
