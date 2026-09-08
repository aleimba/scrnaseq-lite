# scrnaseq-lite: Citations

## Pipeline

> Leimbach A. "scrnaseq-lite: a minimal Nextflow pipeline for droplet-based
> single-cell RNA-seq." https://github.com/aleimba/scrnaseq-lite

## Workflow engine

- [Nextflow](https://doi.org/10.1038/nbt.3820)

> Di Tommaso P, Chatzou M, Floden EW, Prieto Barja P, Palumbo E, Notredame C.
> "Nextflow enables reproducible computational workflows." Nat Biotechnol.
> 2017;35(4):316-319. doi: 10.1038/nbt.3820.

## Pipeline tools

- [SeqKit2](https://doi.org/10.1002/imt2.191)

> Shen W, Sipos B, Zhao L. "SeqKit2: A Swiss army knife for sequence and
> alignment processing." iMeta. 2024;3(3):e191. doi: 10.1002/imt2.191.

> Shen W, Le S, Li Y, Hu F. "SeqKit: A cross-platform and ultrafast toolkit for
> FASTA/Q file manipulation." PLOS ONE. 2016;11(10):e0163962.
> doi: 10.1371/journal.pone.0163962.

- [fastp](https://doi.org/10.1093/bioinformatics/bty560)

> Chen S, Zhou Y, Chen Y, Gu J. "fastp: an ultra-fast all-in-one FASTQ
> preprocessor." Bioinformatics. 2018;34(17):i884-i890.
> doi: 10.1093/bioinformatics/bty560.

- [simpleaf](https://doi.org/10.1093/bioinformatics/btad614)

> He D, Patro R. "simpleaf: a simple, flexible, and scalable framework for
> single-cell data processing using alevin-fry." Bioinformatics.
> 2023;39(10):btad614. doi: 10.1093/bioinformatics/btad614.

- [alevin-fry](https://doi.org/10.1038/s41592-022-01408-3)

> He D, Zakeri M, Sarkar H, Soneson C, Srivastava A, Patro R. "Alevin-fry
> unlocks rapid, accurate and memory-frugal quantification of single-cell
> RNA-seq data." Nat Methods. 2022;19(3):316-322.
> doi: 10.1038/s41592-022-01408-3.

- [piscem](https://github.com/COMBINE-lab/piscem)

> COMBINE-lab. "piscem: mapping sequencing data to De Bruijn graphs, fast."
> https://github.com/COMBINE-lab/piscem
>
> The read mapper and index used by simpleaf. It has no separate publication at
> the time of writing; cite the simpleaf and alevin-fry papers alongside it.

- [QCatch](https://doi.org/10.1093/bioinformatics/btag184)

> Gao Y, He D, Patro R. "QCatch: a framework for quality control assessment and
> analysis of single-cell sequencing data." Bioinformatics. 2026;42(5):btag184.
> doi: 10.1093/bioinformatics/btag184.

- [EmptyDrops](https://doi.org/10.1186/s13059-019-1662-y)

> Lun ATL, Riesenfeld S, Andrews T, Dao TP, Gomes T, participants in the 1st
> Human Cell Atlas Jamboree, Marioni JC. "EmptyDrops: distinguishing cells from
> empty droplets in droplet-based single-cell RNA sequencing data." Genome
> Biol. 2019;20:63. doi: 10.1186/s13059-019-1662-y.
>
> The ambient-RNA model that the cell-calling step is based on.

- [Scanpy](https://doi.org/10.1186/s13059-017-1382-0)

> Wolf FA, Angerer P, Theis FJ. "SCANPY: large-scale single-cell gene
> expression data analysis." Genome Biol. 2018;19:15.
> doi: 10.1186/s13059-017-1382-0.

- [Scrublet](https://doi.org/10.1016/j.cels.2018.11.005)

> Wolock SL, Lopez R, Klein AM. "Scrublet: Computational identification of cell
> doublets in single-cell transcriptomic data." Cell Syst. 2019;8(4):281-291.e9.
> doi: 10.1016/j.cels.2018.11.005.

- [Leiden](https://doi.org/10.1038/s41598-019-41695-z)

> Traag VA, Waltman L, van Eck NJ. "From Louvain to Leiden: guaranteeing
> well-connected communities." Sci Rep. 2019;9:5233.
> doi: 10.1038/s41598-019-41695-z.

- [UMAP](https://arxiv.org/abs/1802.03426)

> McInnes L, Healy J, Melville J. "UMAP: Uniform Manifold Approximation and
> Projection for Dimension Reduction." arXiv:1802.03426 (2018).

- [MultiQC](https://doi.org/10.1093/bioinformatics/btw354)

> Ewels P, Magnusson M, Lundin S, Käller M. "MultiQC: summarize analysis
> results for multiple tools and samples in a single report." Bioinformatics.
> 2016;32(19):3047-3048. doi: 10.1093/bioinformatics/btw354.

- [Quarto](https://doi.org/10.5281/zenodo.5960048)

> Allaire JJ, Teague C, Scheidegger C, Xie Y, Dervieux C, Woodhull G. "Quarto."
> doi: 10.5281/zenodo.5960048. https://quarto.org

## Related pipelines

- [nf-core/scrnaseq](https://nf-co.re/scrnaseq)

> The community-curated single-cell RNA-seq pipeline. It covers a far wider
> range of assays, aligners and options than this pipeline does and is the
> right choice for most production work. Several vendored modules used here
> come from the nf-core module repository, and the barcode whitelists in
> `assets/whitelist/` were taken from nf-core/scrnaseq 4.2.0.

- [nf-core](https://doi.org/10.1038/s41587-020-0439-x)

> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU,
> Di Tommaso P, Nahnsen S. "The nf-core framework for community-curated
> bioinformatics pipelines." Nat Biotechnol. 2020;38(3):276-278.
> doi: 10.1038/s41587-020-0439-x.

## Software packaging and containerisation

- [Bioconda](https://doi.org/10.1038/s41592-018-0046-7)

> Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R,
> Köster J; Bioconda Team. "Bioconda: sustainable and comprehensive software
> distribution for the life sciences." Nat Methods. 2018;15(7):475-476.
> doi: 10.1038/s41592-018-0046-7.

- [BioContainers](https://doi.org/10.1093/bioinformatics/btx192)

> da Veiga Leprevost F, Grüning B, Alves Aflitos S, et al. "BioContainers: an
> open-source and community-driven framework for software standardization."
> Bioinformatics. 2017;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192.

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

> Merkel D. "Docker: lightweight Linux containers for consistent development
> and deployment." Linux Journal. 2014;2014(239):2.

## Data

The `demo` profile uses two public datasets provided by 10x Genomics, both
licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/):

> 1k PBMCs from a Healthy Donor (v3 chemistry), Single Cell Gene Expression
> Dataset by Cell Ranger 3.0.0, 10x Genomics, 19 November 2018.
> https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0

> 10k PBMCs from a Healthy Donor (v3 chemistry), Single Cell Gene Expression
> Dataset by Cell Ranger 3.0.0, 10x Genomics, 19 November 2018.
> https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0

The `test` profile uses data from
[nf-core/test-datasets](https://github.com/nf-core/test-datasets).

Reference genomes and annotation for human data are distributed by
[10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads).
