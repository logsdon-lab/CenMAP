.PHONY: venv dev update_submodules

bin=$(CURDIR)/venv/bin

venv:
	python3 -m venv venv
	$(bin)/python3 -m pip install -r requirements.txt

dev:
	$(MAKE) venv
	$(bin)/python3 -m pip install -r requirements-dev.txt

update_submodules:
	git submodule init
	git submodule update --remote --recursive
