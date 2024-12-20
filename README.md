# `CenMAP`
[![CI](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml/badge.svg)](https://github.com/logsdon-lab/hgsvc3/actions/workflows/main.yml)
[![GitHub Release](https://img.shields.io/github/v/release/logsdon-lab/CenMAP)](https://github.com/logsdon-lab/CenMAP/releases)

A centromere mapping and annotation pipeline for T2T human genome assemblies implemented in [`Snakemake`](https://snakemake.github.io/).

<table>
  <tr>
    <td>
      <figure float="center">
        <img align="middle" src="docs/HG02106_rc-chr1_haplotype2-0000103:121087497-126026281_large.tri.png" width="100%">
        <figcaption>HGSVC3 chr1 centromere HOR structure, centromere dip regions, and self-identity plot.</figcaption>
      </figure>
    </td>
    <td>
      <figure float="left">
        <img align="middle" src="docs/all_cens_chr12_small.png" width="100%">
        <figcaption>HGSVC3 chr12 centromere HOR structure.</figcaption>
      </figure>
      <figure float="left">
        <img align="middle" src="docs/all_AS-HOR_lengths.png" width="100%">
        <figcaption>HGSVC3 Cumulative alpha-satellite HOR array lengths by chromosome.</figcaption>
      </figure>
    </td>
  </tr>
</table>

### [Input](https://github.com/logsdon-lab/CenMAP/wiki/2.-Getting-Started#data)
* [`Verkko`](https://github.com/marbl/verkko) or [`hifiasm`](https://github.com/chhylp123/hifiasm) human genome assemblies
* PacBio HiFi reads used in the assemblies
* [`CHM13`](https://github.com/marbl/CHM13) reference genome assembly
* (Optional) Unaligned BAM files with 5mC modifications at CpG sites.

### [Output](https://github.com/logsdon-lab/CenMAP/wiki/5.-Output)
* Complete and correctly assembled centromere sequences and their regions.
* Centromere alpha-satellite higher order repeat (HOR) array lengths.
* [`RepeatMasker`](https://www.repeatmasker.org/) and [`HumAS-HMMER`](https://github.com/enigene/HumAS-HMMER) alpha-satellite HOR monomer annotations and plots.
* [`ModDotPlot`](https://github.com/marbl/ModDotPlot) sequence identity plots.
* Combined sequence identity and HOR array structure plots.
* (Optional) Centromere dip region (CDRs) with [`CDR-Finder`](https://github.com/koisland/CDR-Finder)

### [Documentation](https://github.com/logsdon-lab/CenMAP/wiki)
Read the docs on the `CenMAP` [wiki](https://github.com/logsdon-lab/CenMAP/wiki).

### [Tests](https://github.com/logsdon-lab/CenMAP/wiki/6.-Test)
To run tests, refer to the wiki [page](https://github.com/logsdon-lab/CenMAP/wiki/6.-Test).
