# `CenMAP`
[![CI](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml)
[![GitHub Release](https://img.shields.io/github/v/release/logsdon-lab/CenMAP)](https://github.com/logsdon-lab/CenMAP/releases)

A centromere mapping and annotation pipeline implemented in [`Snakemake`](https://snakemake.github.io/).

<table>
  <tr>
    <th>
      <figure float="left">
        <img align="middle" src="docs/all_cens_chr12_small.png" width="100%">
        <figcaption>chr12 centromere HOR structure.</figcaption>
      </figure>
    </th>
    <th>
      <figure float="left">
        <img align="middle" src="docs/all_AS-HOR_lengths.png" width="100%">
        <figcaption>Cumulative alpha-satellite HOR array lengths by chromosome.</figcaption>
      </figure>
    </th>
  </tr>
</table>

### [Input](https://github.com/logsdon-lab/CenMAP/wiki/2.-Getting-Started#data)
* `Verkko`/`hifiasm` genome assemblies
* PacBio HiFi reads used in the assemblies.
* A reference genome, ex. `CHM13`.

### [Output](https://github.com/logsdon-lab/CenMAP/wiki/5.-Output)
* Complete and correctly assembled centromere sequences and their regions.
* Centromere alpha-satellite higher order repeat (HOR) array lengths.
* `RepeatMasker` and `HumAS-HMMER` alpha-satellite HOR annotations and plots.
* `ModDotPlot` sequence identity plots.

### [Documentation]((https://github.com/logsdon-lab/CenMAP/wiki))
Read the docs on the `CenMAP` [wiki](https://github.com/logsdon-lab/CenMAP/wiki).
