# Annotations

### HumAS StV
* `AS-HOR-vs-chm1_cens_v21.stv_row.all.bed`
    * HumAS-HMMER based stv annotations on chm1 reference.
* `AS-HOR-vs-chm13_cens_v18.stv_row.all.bed`
    * HumAS-HMMER based stv annotations on chm13 reference.

> Can be passed to `config["plot_hor_stv"]["ref_stv"]` for visualization.

### RepeatMasker
* `chm13_chm1_cens_v21.trimmed.fa.noheader.out`
    * RepeatMasker output of chm1 and chm13 chromosomes. Includes chrY.
    * No header.

> Can be passed to `config["repeatmasker"]["ref_repeatmasker_output"]` for visualization.

### Alpha-satellite HOR array lengths
* `chm1_hor_length.tsv`
    * chm1 reference alpha-satellite HOR lengths.
* `chm13_hor_length.tsv`
    * chm13 reference alpha-satellite HOR lengths.

> Can be passed to `config["calculate_hor_length"]["ref_hor_lengths"]` for visualization.
