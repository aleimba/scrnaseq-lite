# scrnaseq-lite

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A5%2026.04-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed.svg)](https://www.docker.com/)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE.md)

A Nextflow pipeline that takes droplet-based single-cell RNA-seq
FASTQs to a clustered, annotated AnnData object plus QC and analysis reports.

> **IMPORTANT:** Pre-release. The pipeline is under active development and not yet complete.

## Contents

- [scrnaseq-lite](#scrnaseq-lite)
  - [Contents](#contents)
  - [Purpose](#purpose)
  - [Workflow](#workflow)
  - [Quick start](#quick-start)
  - [Input](#input)
  - [Test and demo data](#test-and-demo-data)
  - [Design decisions](#design-decisions)
  - [Development](#development)
  - [Citations](#citations)
  - [Licence and contributing](#licence-and-contributing)

## Purpose

`scrnaseq-lite` runs **one analysis workflow**, end to end, with no branches to
choose between: reads in, clusters and marker genes out.

It is built for:

- Droplet-based 10x Genomics 3' single-cell RNA-seq
- Chemistries `10XV2`, `10XV3` and `10XV4`
- Paired-end FASTQs, one library per sample, singleplex only
- One species and one reference per run

The design goal is maintainability and reviewability, e.g. there is a single
quantification path rather than a menu of
aligners. For projects that need broader assay support, multiplexed libraries
or a community-maintained codebase, use
[nf-core/scrnaseq](https://nf-co.re/scrnaseq).

## Workflow

```
samplesheet.csv
      |
      +--> read QC ---------------------------------------+
      |      (also asserts R1 length matches chemistry)    |
      v                                                    |
  quantification against a splici index                    |
      |   raw count matrix, no cell calling applied        |
      v                                                    |
  cell calling (ambient-RNA model) --> per-sample report --+
      |   filtered count matrix                            |
      v                                                    |
  cell QC and doublet detection --------------------------+
      |                                                    |
      v                                                    v
  normalise -> HVG -> PCA -> UMAP -> Leiden clustering    QC report
      |
      v
  marker genes  -->  analysis report
```

1. **Read QC.** Per-sample read statistics and a paired-end QC report. Read 1
   is never trimmed or filtered: it carries the cell barcode and UMI at fixed
   offsets. The observed read 1 length is asserted against the chemistry
   declared in the samplesheet, and a mismatch fails the run.
2. **Reference.** A spliced + intronic (splici) index is built once from a
   genome FASTA and a GTF, then reused across runs.
3. **Quantification.** Barcode correction against the full 10x whitelist and
   UMI resolution, producing a raw count matrix with **no cell calling
   applied**.
4. **Cell calling.** Cells are separated from empty droplets with a two-step
   EmptyDrops-style ambient-RNA model, not a knee or a fixed UMI cut-off.
   PBMCs carry little RNA, and a global threshold discards real low-count
   cells.
5. **Cell QC and doublets.** MAD-based outlier filtering with absolute floors
   and ceilings, then seeded doublet detection. Doublets are removed before
   feature selection and PCA.
6. **Normalisation and clustering.** Normalise, log-transform, select highly
   variable genes on raw counts, PCA, neighbours, UMAP and Leiden clustering.
   Raw counts are preserved in a dedicated layer throughout.
7. **Marker genes.** Ranked cluster-vs-rest marker detection. This is currently
   marker detection, **not** differential expression between conditions.
8. **Reports.** A QC report across samples, one cell-calling report per sample,
   and an analysis report carrying the biological result.

Every table shown in a report is also written as a standalone TSV, and every
figure as a PNG. Tool versions, resolved parameters and a run manifest are
written with the results. See [docs/usage.md](docs/usage.md) for parameters and
[CITATIONS.md](CITATIONS.md) for the tools behind each step.

## Quick start

Requires [Nextflow](https://www.nextflow.io/) 26.04 or newer and Docker.

```bash
# Small remote test data. Exercises the workflow in minutes.
nextflow run . -profile test,docker -params-file params/test.yaml

# Your own data.
nextflow run . -profile docker -params-file params/default.yaml \
    --input samplesheet.csv --run my_run_01
```

**IMPORTANT**:

- **`-params-file` is mandatory on every run.** Parameter files are not loaded
  automatically, and unset parameters are not filled in from defaults.
- **`-resume` needs `--run <name>`.** Without a run name the output directory
  is created with a new timestamp at each launch, so a resumed run would write
  into a new directory.

Full instructions, including how to obtain a reference and build the index, are
in [docs/usage.md](docs/usage.md).

## Input

A CSV samplesheet with one row per sample, validated before a run
starts (see [nf-schema](https://nextflow-io.github.io/nf-schema/latest/)). See [assets/samplesheet_example_valid.csv](assets/samplesheet_example_valid.csv).

```csv
sample,fastq_1,fastq_2,chemistry,expected_cells
pbmc_1k_v3,data/pbmc_1k_v3_R1.fastq.gz,data/pbmc_1k_v3_R2.fastq.gz,10XV3,1222
```

| Column           | Required | Description                                                     |
| ---------------- | -------- | --------------------------------------------------------------- |
| `sample`         | yes      | Unique sample identifier; becomes the output prefix              |
| `fastq_1`        | yes      | Gzipped read 1 (cell barcode + UMI)                              |
| `fastq_2`        | yes      | Gzipped read 2 (cDNA)                                            |
| `chemistry`      | yes      | One of `10XV2`, `10XV3`, `10XV4`                                 |
| `expected_cells` | no       | Reporting only; never passed to the quantifier                   |

`chemistry` selects the quantifier chemistry, the cell-calling chemistry,
the barcode whitelist and the expected read 1 length together.

Technical replicates must be merged before they reach the pipeline;
currently it does not concatenate FASTQs.

## Test and demo data

**`test`** uses a small remote fixture from
[nf-core/test-datasets](https://github.com/nf-core/test-datasets) (mouse GRCm38
chromosome 19, 10x 3' v2) and exercises the whole workflow in minutes with no
local data management. It is far too small for an ambient-RNA model, so it
substitutes a crude count threshold for cell calling. **Its results are not
biologically meaningful, but useful to testing a small-scale pipeline run.**

**`demo`** uses two public human PBMC datasets from 10x Genomics, downsampled
to roughly 24 million read pairs each by
`bin/download_and_downsample_testdata.sh`. This is the profile that produces
presentable results.

| Metric              | `pbmc_1k_v3`  | `pbmc_10k_v3` |
| ------------------- | ------------- | ------------- |
| Estimated cells     | 1,222         | 11,769        |
| Mean reads/cell     | 54,502        | 54,286        |
| Median genes/cell   | 1,919         | 1,906         |
| Read 1 / read 2     | 28 bp / 91 bp | 28 bp / 91 bp |

Both come from the same healthy donor at the same chemistry and are nearly
identical in sequencing depth and genes per cell. The real contrast between
them is cell loading, roughly 1,200 against 11,800 recovered cells, which on
the 10x v3 loading curve means markedly different doublet rates. The
downsampling preserves that difference: the 10k sample is subset by cell
barcode rather than by read, so per-cell depth is retained and only the cell
count changes.

Data provided by 10x Genomics:

- [1k PBMCs from a Healthy Donor (v3 chemistry)](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0),
  Single Cell Gene Expression Dataset by Cell Ranger 3.0.0, 10x Genomics,
  19 November 2018.
- [10k PBMCs from a Healthy Donor (v3 chemistry)](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0),
  Single Cell Gene Expression Dataset by Cell Ranger 3.0.0, 10x Genomics,
  19 November 2018.

Both are licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Design decisions

The choices that most shape the codebase, and the reasoning behind each:

- **Publishing is declared once, centrally.** No process writes to the results
  directory itself. Everything is published through the entry workflow, so the
  complete list of outputs can be read in one place. Publish mode is `copy`,
  never `symlink`, because object storage has no symlinks.
- **Every value is defined in exactly one file.** Run-variable science lives in
  `params/*.yaml`, per-process arguments in `conf/modules.config`, resources in
  `conf/base.config`, engine and profile settings in `nextflow.config`. Each
  parameter file is a complete set rather than a sparse overlay, so it doubles
  as a self-describing record of one run.
- **No scripts are embedded in workflow files.** Every Python and shell script
  is a real executable file in `bin/`, so it can be linted, run by hand outside
  Nextflow, and diffed sensibly.
- **Containers are per tool, and overrides are quarantined.** Temporary
  container overrides live alone in `conf/containers.config`, with the reason
  and the deletion condition written next to each, so they stay easy to see and
  to remove.
- **Barcode whitelists are vendored in the repository**
  ([assets/whitelist/](assets/whitelist/)). Otherwise a
  silent network dependency in the middle of a run, which fails on an offline
  node or inside a locked-down private subnet.
- **The count matrix is reconstructed explicitly.** The quantifier stores
  unspliced, spliced and ambiguous counts as separate layers and writes their
  sum as the main matrix. The layer combination used for analysis is a
  parameter and is rebuilt from the named layers rather than read from the sum.

## Development

This pipeline is built spec-first with [Claude Code](https://claude.com/claude-code)
with human-in-the-loop, in an explicit and auditable order: requirements, then
design, then a task list worked one task at a time. Constraints are written down
before any code exists, and findings verified by running the pinned containers
are written back into the specs, so the specs stay the source of truth instead
of drifting behind the implementation.

The specification documents are part of the repository and readable on their
own:

| Document                                                     | Contents                                                        |
| ------------------------------------------------------------ | --------------------------------------------------------------- |
| [01-requirements.md](.claude/specs/01-requirements.md)        | What the pipeline must do, and the scientific constraints        |
| [02-design.md](.claude/specs/02-design.md)                    | Dataflow, module inventory, configuration layering               |
| [03-tasks.md](.claude/specs/03-tasks.md)                      | The implementation task list, with findings recorded per task    |

[CLAUDE.md](CLAUDE.md) is the standing contract the agent works under.

## Citations

Tool and data references are listed in [CITATIONS.md](CITATIONS.md).

## Licence and contributing

Released under the [GPL-3.0 licence](LICENSE.md).
See [CONTRIBUTING.md](CONTRIBUTING.md) and
[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
