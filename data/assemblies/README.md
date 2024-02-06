# Assemblies
Each subdirectory should contain assemblies for each sample for the following types:
- mito
- rdna
- disconnected
- ebv
- hap1
- hap2
- unassigned
- contaminants

```bash
# tree data/assemblies/ -P "*.fasta.gz"
data/assemblies/
└── HG00171
    ├── HG00171.ps-sseq.exemplar-ebv.fasta.gz
    ├── HG00171.ps-sseq.exemplar-mito.fasta.gz
    ├── HG00171.ps-sseq.exemplar-rdna.fasta.gz
    ├── HG00171.vrk-ps-sseq.asm-disconnected.fasta.gz
    ├── HG00171.vrk-ps-sseq.asm-ebv.fasta.gz
    ├── HG00171.vrk-ps-sseq.asm-hap1.fasta.gz
    ├── HG00171.vrk-ps-sseq.asm-hap2.fasta.gz
    ├── HG00171.vrk-ps-sseq.asm-mito.fasta.gz
    ├── HG00171.vrk-ps-sseq.asm-rdna.fasta.gz
    ├── HG00171.vrk-ps-sseq.asm-unassigned.fasta.gz
    └── HG00171.vrk-ps-sseq.contaminants.fasta.gz

1 directory, 11 files
```

The format of the filenames must follow:
* `{sample}.vrk-ps-sseq.asm-{type}.fasta.gz`
* Contaminants are an exception and are named:
    * `{sample}.vrk-ps-sseq.contaminants.fasta.gz`
