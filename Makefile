.PHONY: dockerfile singularity dockerhub venv dev

tag_base="logsdonlab"
dir_name="$(shell basename $(CURDIR))"
local_tag=$(tag_base)/$(dir_name):latest
final_tag=$(account)/$(dir_name):latest
bin=$(CURDIR)/venv/bin

dockerfile:
	sudo docker build docker/ -t $(local_tag)
	sudo docker tag $(local_tag) $(final_tag)

singularity:
	$(MAKE) dockerfile account=$(account)
	sudo singularity build $(dir_name).sif docker-daemon://$(final_tag)

dockerhub:
	$(MAKE) dockerfile account=$(account)
	sudo docker image push $(final_tag)

venv:
	python3 -m venv venv
	$(bin)/python3 -m pip install -r requirements.txt

dev:
	$(MAKE) venv
	$(bin)/python3 -m pip install -r requirements-dev.txt
