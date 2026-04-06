# autoretro_ab

`autoretro_ab` is a sharable workflow for automated retrotransposon RNA-seq analysis.

By:

Abrar Alghamdi 


It supports:

- downloading SRR runs
- alignment with Bowtie2
- BAM generation with Samtools
- separate Telescope analysis for HERV and L1
- merging Telescope counts
- filtering count matrices
- DESeq2 differential expression
- top hits extraction
- family enrichment
- chromosomal distribution

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/autoretro_ab.git
cd autoretro_ab
pip install -e .
