# Docker/Singularity
Why?
* Modules not as reproducible.
* Certain tools cannot be installed via conda (ex. `dna-brnn`).

### Build
Build a local docker image and convert it to a singularity image.
* TODO: Should this be run from parent dir and entire workflow containerized?
```bash
tag="logsdonlab/hgsvc3:latest
sudo docker build -t $tag .
sudo singularity build hgsvc3.sif docker-daemon://$tag
```

Interactive shell.
```bash
singularity shell hgsvc3.sif
```

### Source
* https://micromamba-docker.readthedocs.io/en/latest/quick_start.html#running-commands-in-dockerfile-within-the-conda-environment
