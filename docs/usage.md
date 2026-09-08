# Usage

How to run `scrnaseq-lite`. This document describes the pipeline as currently
implemented; see [../README.md](../README.md) for an overview.

> **IMPORTANT:** Pre-release. The pipeline is under active development and not yet complete.

## Contents

- [Usage](#usage)
  - [Contents](#contents)
  - [Requirements](#requirements)
  - [Parameter files](#parameter-files)
  - [Samplesheet](#samplesheet)
  - [Reference and index](#reference-and-index)
  - [Running the pipeline](#running-the-pipeline)
  - [Run names and resume](#run-names-and-resume)
  - [Profiles](#profiles)
  - [Parameters](#parameters)
    - [Input and reference](#input-and-reference)
    - [Quantification](#quantification)
    - [Cell calling](#cell-calling)
    - [Cell QC and doublets](#cell-qc-and-doublets)
    - [Normalisation, clustering and markers](#normalisation-clustering-and-markers)
    - [Reproducibility and reporting](#reproducibility-and-reporting)
  - [Troubleshooting](#troubleshooting)

## Requirements

- [Nextflow](https://www.nextflow.io/) 26.04 or newer
- Docker

All pipeline software runs in containers. Nothing needs to be installed on the
host beyond Nextflow and Docker.

Resource needs are modest. The heaviest step is building the index, and the
full human splici index is around 2 GB. The `demo` profile is capped at 2 CPUs
and 8 GB.

## Parameter files

**Every run must pass `-params-file`.** Parameter files are not loaded
automatically, and values left unset are not filled in from schema defaults, so
a run without one will fail or behave unexpectedly.

Three complete parameter sets ship with the pipeline:

| File                    | Purpose                                           |
| ----------------------- | ------------------------------------------------- |
| `params/default.yaml`   | Documented defaults; the starting point for real data |
| `params/test.yaml`      | Small-scale remote data set for testing the pipeline |
| `params/demo.yaml`      | Downsampled PBMC data                         |

Each file is a **complete** parameter, so it is a self-describing record of
one run's parameters. To adapt one, copy it and edit the values; do not delete keys.

Values from `-params-file` outrank anything set in a config file, including an
explicit `null`.

## Samplesheet

A CSV with one row per sample, validated before a run starts (see [nf-schema](https://nextflow-io.github.io/nf-schema/latest/)).

```csv
sample,fastq_1,fastq_2,chemistry,expected_cells
pbmc_1k_v3,data/pbmc_1k_v3_R1.fastq.gz,data/pbmc_1k_v3_R2.fastq.gz,10XV3,1222
pbmc_10k_v3,data/pbmc_10k_v3_R1.fastq.gz,data/pbmc_10k_v3_R2.fastq.gz,10XV3,
```

| Column           | Required | Description                                                   |
| ---------------- | -------- | ------------------------------------------------------------- |
| `sample`         | yes      | Unique identifier; no spaces. Becomes the output prefix        |
| `fastq_1`        | yes      | Gzipped read 1, carrying cell barcode and UMI                  |
| `fastq_2`        | yes      | Gzipped read 2, the cDNA read                                  |
| `chemistry`      | yes      | `10XV2`, `10XV3` or `10XV4`                                    |
| `expected_cells` | no       | Positive integer, or empty                                     |

Notes:

- **Sample names must be unique.** Technical replicates of one library must be
  merged before they reach the pipeline; it currently does not concatenate FASTQs.
- **Read 2 is required.** Droplet-based 10x libraries are always paired, so a
  single-end row cannot be quantified.
- **`expected_cells` is metadata only.** It is not passed to the quantifier.
  The permit-list modes that would consume it are mutually exclusive with the
  unfiltered mode this pipeline uses, and cell calling is a separate step. The
  value is carried for reporting and for comparison against the number of cells
  actually called.
- **`chemistry` drives everything derived from it**: the quantifier chemistry,
  the cell-calling chemistry, the barcode whitelist, and the expected read 1
  length (26 bp for `10XV2`, 28 bp for `10XV3` and `10XV4`). They cannot be set
  independently, because a run whose whitelist disagrees with its chemistry
  produces a plausible-looking but wrong matrix rather than an error.

The observed read 1 length is checked against the declared chemistry and a
mismatch fails the run. That is deliberate: a wrong chemistry is otherwise
silent.

## Reference and index

There is no downloadable pre-built index for this quantifier, so the index is
built once and reused.

1. **Fetch a reference.** For human data:

   ```bash
   bin/download_reference.sh --outdir reference
   ```

   This downloads the 10x Genomics GRCh38 2024-A package, which supplies the
   genome FASTA and the gene annotation GTF. The GTF is used gzipped; it does
   not need decompressing.

2. **Build the index once.** Supplying both a FASTA and a GTF is what selects
   the spliced + intronic (splici) strategy (adapt default.yaml with your settings).

   ```bash
   nextflow run . -profile docker -params-file params/default.yaml \
       --build_index --save_reference \
       --fasta reference/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
       --gtf   reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz \
       --run index_build
   ```

3. **Reuse it in every later run** by pointing `simpleaf_index` at the saved
   directory, either in your parameter file or on the command line:

   ```bash
   --simpleaf_index index/GRCh38-2024-A_splici
   ```

Set `r2_read_length` to the cDNA read length of the data the index will be used
on. It determines the flank length added to intronic sequence, so an index
built for one read length is not ideal for another. Both bundled PBMC datasets
are 91 bp, so one index serves both.

## Running the pipeline

```bash
# Small remote test data, in minutes.
nextflow run . -profile test,docker -params-file params/test.yaml

# The downsampled PBMC demo.
nextflow run . -profile demo,docker -params-file params/demo.yaml --run demo01

# Your own data (adapt default.yaml with your settings)
nextflow run . -profile docker -params-file params/default.yaml \
    --input samplesheet.csv \
    --simpleaf_index index/GRCh38-2024-A_splici \
    --run my_run_01
```

A stub run checks the pipeline workflow without doing any real work:

```bash
nextflow run . -profile test -stub-run -params-file params/test.yaml
```

## Run names and resume

`--run <name>` names the run and fixes its output directory to
`results/<name>`. Without it, the directory is `results/run_<timestamp>` and
**the timestamp is re-evaluated at every launch**.

The practical consequence: **`-resume` only produces a coherent output
directory when `--run <name>` is given.** A resumed run without a name reuses
the cached work but publishes into a brand new directory.

**Best practice**: `--run` is a command-line control only. It should not be placed
in a parameter file, because every run sharing that file would then write into the same
output directory.

Note also that `-output-dir` is a single-dash Nextflow core argument, not a pipeline
parameter, and `--outputDir` does not exist. Use `--run`.

## Profiles

| Profile   | Purpose                                                              |
| --------- | -------------------------------------------------------------------- |
| `docker`  | Run all tools in Docker containers                                    |
| `test`    | Small remote fixture; engine caps only. Results not meaningful        |
| `demo`    | Downsampled PBMC data; engine caps only                               |
| `awsbatch`| AWS Batch execution settings                                          |

Profiles are combined with a comma and no space: `-profile test,docker`.

`test` and `demo` carry engine settings only and set no parameters at all;
their parameters come from the matching file in `params/`.

## Parameters

Every parameter is documented inline in `params/default.yaml`, which is the
authoritative reference. The groups below summarise what each controls.

### Input and reference

| Parameter        | Default | Description                                         |
| ---------------- | ------- | --------------------------------------------------- |
| `input`          | `null`  | Samplesheet CSV                                     |
| `build_index`    | `false` | Build the splici index in this run                  |
| `simpleaf_index` | `null`  | Path to an existing index directory                 |
| `save_reference` | `false` | Publish the built index for reuse                   |
| `fasta`          | `null`  | Genome FASTA; with `gtf`, selects the splici strategy |
| `gtf`            | `null`  | Gene annotation, gzipped or plain                   |
| `r2_read_length` | `91`    | cDNA read length; sets the intron flank size        |

### Quantification

| Parameter          | Default          | Description                             |
| ------------------ | ---------------- | --------------------------------------- |
| `umi_resolution`   | `cr-like`        | UMI resolution strategy                 |
| `permit_list_mode` | `unfiltered-pl`  | Barcode correction against the whitelist, with no cell calling |
| `count_layer`      | `S+A`            | Layers combined into the analysis matrix |

`count_layer` deserves a note. The quantifier stores unspliced, spliced and
ambiguous counts as separate layers and writes their **sum** as the main
matrix. The pipeline reconstructs `count_layer` from the named layers rather
than reading the sum. Use `S+A` for whole cells and `S+A+U` for single nuclei.

### Cell calling

| Parameter             | Default  | Description                                  |
| --------------------- | -------- | -------------------------------------------- |
| `cell_calling`        | `qcatch` | `qcatch` or `threshold`                      |
| `qcatch_n_partitions` | `null`   | Partition count for custom assays only       |

`cell_calling: threshold` applies **no ambient model at all**. It exists so a
fixture too small for one can still exercise the workflow end to end, and it is
used only by the `test` profile. It is not a faster or simpler alternative: a
global cut discards real low-count cells. Results produced with it are not
publishable and should not be described as cell-called.

### Cell QC and doublets

| Parameter              | Default              | Description                        |
| ---------------------- | -------------------- | ---------------------------------- |
| `min_genes`            | `200`                | Debris floor                       |
| `min_cells_per_gene`   | `3`                  | Drop genes seen in fewer cells     |
| `mito_mad`             | `5`                  | MAD cut on mitochondrial percentage |
| `mito_max_pct`         | `20`                 | Hard ceiling on mitochondrial percentage |
| `count_mad`            | `5`                  | MAD cut on log1p total counts      |
| `min_cells_per_sample` | `100`                | Fail readably below this           |
| `mito_gene_prefix`     | `MT-`                | Mitochondrial gene pattern         |
| `ribo_gene_pattern`    | `^RP[SL]`            | Ribosomal gene pattern             |
| `hb_gene_pattern`      | `^HB[ABDEGQZ][0-9]?$`| Haemoglobin gene pattern           |
| `remove_doublets`      | `true`               | Drop predicted doublets before feature selection |
| `expected_doublet_rate`| `0.08`               | Prior doublet rate                 |

Thresholds are MAD-based with absolute floors and ceilings, so no cut-off is
hard-coded in a script. The gene patterns are parameters so a non-human
reference can be used.

`expected_doublet_rate` is a single global prior suiting a capture of roughly
10,000 cells; a 1,000-cell capture sits nearer 0.008. The detector is driven
mainly by its simulated-doublet score distribution, so the prior shifts the
threshold rather than setting it.

### Normalisation, clustering and markers

| Parameter           | Default     | Description                              |
| ------------------- | ----------- | ---------------------------------------- |
| `target_sum`        | `null`      | `null` uses the median library size      |
| `n_hvg`             | `2000`      | Highly variable genes to select          |
| `hvg_flavor`        | `seurat_v3` | Operates on raw counts                   |
| `n_pcs`             | `50`        | Principal components                     |
| `n_neighbors`       | `15`        | Neighbourhood size                       |
| `leiden_resolution` | `1.0`       | Clustering resolution                    |
| `umap_min_dist`     | `0.5`       | UMAP minimum distance                    |
| `marker_method`     | `wilcoxon`  | Marker ranking test                      |
| `n_markers_plot`    | `5`         | Top genes per cluster in the dotplot     |

Marker detection compares each cluster against the rest. It is **not**
differential expression between experimental conditions.

### Reproducibility and reporting

| Parameter               | Default                    | Description                  |
| ----------------------- | -------------------------- | ---------------------------- |
| `seed`                  | `42`                       | Reaches doublet detection, PCA, UMAP and clustering |
| `multiqc_title`         | `scrnaseq-lite QC report`  | QC report title              |
| `skip_multiqc`          | `false`                    | Skip the QC report           |
| `skip_analysis_report`  | `false`                    | Skip the analysis report     |

## Troubleshooting

**The run stops complaining about a read 1 length mismatch.** This is working
as intended. The observed read 1 length does not match the chemistry declared
in the samplesheet. Correct the samplesheet rather than disabling the check;
the wrong chemistry silently produces a wrong matrix.

**A sample yields too few cells and the run fails.** The
`min_cells_per_sample` floor stopped it deliberately, so a sample that lost
nearly everything at QC does not quietly reach clustering. Inspect the QC
report before lowering it.

**Some tools do not run in the container the module declares.** Two container
overrides are in effect, both documented with their reason and their deletion
condition in `conf/containers.config`. One works around an upstream binary that
faults on older CPUs; the other is a straight version upgrade.

**Indexing stalls or fails with I/O errors.** The indexer opens very many
temporary files at once. It is configured to use node-local scratch space for
exactly this reason; on a shared or network filesystem, make sure a real local
scratch directory is available.
