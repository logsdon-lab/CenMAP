.PHONY: conda update_submodules

conda:
	conda env create -f env.yaml --name cenmap
	conda activate cenmap

update_submodules:
	git submodule init
	git submodule update --remote --recursive
