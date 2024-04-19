# `CenMAP`
[![CI](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml)

A centromere mapping and annotation pipeline.

### Input
* `Verkko`/`hifiasm` genome assemblies
* PacBio HiFi reads used in the assemblies.
* A reference genome, ex. `CHM13`.

### Output
* A `bed` file with complete and correctly assembled centromeres.
* A `tsv` file with centromere alpha-satellite higher order repeat (HOR) array lengths.
* `RepeatMasker` repeat annotations and plots.
* `StainedGlass` sequence identity plots.
* `HumAS-HMMER` alpha-satellite HORs annotations and plots.

### Getting Started
`CenMAP` is implemented as a `Snakemake` pipeline and requires limited dependencies:
* `Python >= 3.7`
* `conda` / `mamba`
* `singularity`

#### Clone
```bash
git clone git@github.com:logsdon-lab/CenMAP.git --recurse-submodules

# Create a python virtualenv.
make venv && source venv/bin/activate
```

#### Data
By default, the following are expected.

Inputs can be modified in [`Configuration`](#configuration).

##### Assemblies
* Directory with subdirectories corresponding to sample names each with fasta files.
  * Compressed or uncompressed files with the `fa` or `fasta` extension are supported.
* `concat_asm.input_dir`

```yaml
concat_asm:
  input_dir: "data/assemblies"
```
```
data/assemblies/GM19129/
└── GM19129-hifiasm_v0.19.6.fa.gz
```

##### PacBio HiFi Reads
* Directory with subdirectories corresponding to sample names each with unaligned reads.
  * Compressed or uncompressed files with the `bam`, `fq`, or `fastq` extension are supported.
  * Defaults to reads with a `bam` extension.
  * Modify `nuc_freq.reads_ext` with a list of extensions.

    ```yaml
    nuc_freq:
      reads_ext:
        - bam
        - fq
        - fq.gz
        - fastq
        - fastq.gz
    ```


* `nuc_freq.hifi_reads_dir`

```yaml
nuc_freq:
  hifi_reads_dir: "data/raw_data"
```
```
data/raw_data/HG00171/
├── m54329U_220205_003428.hifi_reads.bam
├── m54329U_220205_003428.hifi_reads.bam.pbi
├── m54329U_220212_223600.hifi_reads.bam
├── m54329U_220212_223600.hifi_reads.bam.pbi
├── m54329U_220214_093304.hifi_reads.bam
└── m54329U_220214_093304.hifi_reads.bam.pbi
```

##### Reference Genome
* Fasta file reference genome.
* `align_asm_to_ref.reference`

```yaml
align_asm_to_ref:
    reference: "data/reference/T2T-CHM13v2.fasta"
```
```
data/reference/
└── T2T-CHM13v2.fasta
```

#### Configuration
Configuration is handled though `config/config.yaml`.

Samples and chromosomes can be specified via the `samples` and `chromosomes` sections, respectively.
```yaml
chromosomes:
  - "chr1"
  - "chr2"
  - "chr3"

samples:
  - HG00171
```

Each section specifies a step in the pipeline.
* Output directories, input reference annotations, threads, and other step-specific parameters can be modified.
* Parameter descriptions can be found in `config/schema/config.schema.yaml`.

### Usage

#### Local
```bash
snakemake -np -c 12 --use-conda --use-singularity
```

#### Cluster
The following commands shows cluster usage:

[`LSF`](https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=overview-lsf-introduction)
```bash
snakemake -np -j 100 \
--cluster "bsub -M {resources.mem_mb} -R 'rusage[mem={resources.mem_mb}]' -n {threads} -o /dev/null" \
--use-conda \
--use-singularity
```

### Steps
`CenMAP` is a centromere mapping and annotation pipeline that does the following.
* Extracts centromere regions from genome assemblies.
* Checks for centromeric contig misassemblies.
    * This is achieved via a [fork of `NucFreq`](https://github.com/logsdon-lab/CenMAP) which flags misassemblies based on per-base read coverage.
* Annotates repeats and alpha-satellite higher-order repeats in centromeres.
    * Using [`RepeatMasker`](https://github.com/rmhubley/RepeatMasker)
* Generates sequence identity plots across centromeres.
    * Using [`StainedGlass`](https://github.com/mrvollger/StainedGlass)

<p float="left">
    <img align="middle" src="docs/rulegraph.svg" width="30%">
</p>

### TODO
* Add test/move calculate HOR length script to new repo.
* Parameter validation.
