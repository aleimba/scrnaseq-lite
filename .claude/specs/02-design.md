# Design — scrnaseq-lite

## 1. High-level dataflow

```
samplesheet.csv
   | nf-schema validation
   v
[ meta, [r1, r2] ]
   |
   +--> SEQKIT_STATS ------> stats.tsv --------------------+
   |         |                                             |
   |         +--> ASSERT_R1_LENGTH  (hard fail on mismatch)
   |                                                       |
   +--> FASTP (report only) --> json + html ---------------+
   |                                                       |
   v                                                       |
SIMPLEAF_QUANT  <-- index (SIMPLEAF_INDEX | --simpleaf_index)
   | *_raw_matrix.h5ad  (unfiltered, full whitelist)       |
   v                                                       |
QCATCH --> *_filtered_matrix.h5ad + report.html + summary.csv
   |                                                       |
   v                                                       |
SCANPY_CELL_QC --> h5ad + PNG     + qc_mqc.tsv ------------+
   |                                                       |
   v                                                       |
SCANPY_DETECT_DOUBLETS --> h5ad + doublet_mqc.tsv ---------+
   |                                                       |
   v                                                       |
SCANPY_NORMALISE_CLUSTER --> h5ad + UMAP PNG               |
   |                                                       |
   v                                                       |
SCANPY_MARKER_GENES --> markers.tsv + dotplot              |
   |                                                       |
   +--> QUARTO_ANALYSIS_REPORT --> analysis_report.html    |
   |                                                       |
   +--> COLLECT_VERSIONS --> versions.tsv                  |
   +--> WRITE_RUN_MANIFEST --> *.yaml                      v
                                              MULTIQC --> multiqc_report.html
```

## 2. Module inventory

| Process                    | Source                            | Notes                                   |
| -------------------------- | --------------------------------- | --------------------------------------- |
| `SEQKIT_STATS`             | nf-core                           | `-a -T`                                 |
| `ASSERT_R1_LENGTH`         | local                             | fails on chemistry mismatch             |
| `FASTP`                    | nf-core                           | report-only; reads not consumed         |
| `SIMPLEAF_INDEX`           | nf-core                           | conditional; `scratch = true`           |
| `SIMPLEAF_QUANT`           | nf-core                           | piscem backend, unfiltered matrix       |
| `QCATCH`                   | nf-core if available, else local  | cell calling + QC HTML                  |
| `SCANPY_CELL_QC`           | local                             | `bin/scanpy_cell_qc.py`                 |
| `SCANPY_DETECT_DOUBLETS`   | local                             | `bin/scanpy_detect_doublets.py`         |
| `SCANPY_NORMALISE_CLUSTER` | local                             | `bin/scanpy_normalise_cluster.py`       |
| `SCANPY_MARKER_GENES`      | local                             | `bin/scanpy_marker_genes.py`            |
| `QUARTO_ANALYSIS_REPORT`   | local                             | `report/analysis_report.qmd`            |
| `COLLECT_VERSIONS`         | local                             | `bin/collect_versions.py`               |
| `WRITE_RUN_MANIFEST`       | local                             | `bin/write_run_manifest.py`             |
| `MULTIQC`                  | nf-core                           | custom content                          |

**Design decision — separate Scanpy processes rather than one monolith.**
Cost: h5ad serialisation between steps. Benefit: `-resume` granularity,
per-step publishing, independent stubbing, and a reviewable DAG. For a
demonstration pipeline that trade-off favours separation.

**Design decision — QCatch for cell calling.** It is the only tool that
reads alevin-fry output natively and models the ambient pool.

**Design decision — no MultiQC module for the quantifier.** MultiQC ships
modules for fastp, SeqKit, Salmon, Kallisto, Bustools and Cell Ranger, but
none for alevin-fry, simpleaf or QCatch. The Salmon module parses
`salmon quant` output, not alevin/piscem mapping logs, so it does not
help. Quantifier QC therefore reaches MultiQC only as custom content
derived from the QCatch summary CSV.

## 3. Channel design

Every channel carries `[ meta, data ]`:

```groovy
meta = [
    id            : 'pbmc_1k_v3',
    single_end    : false,
    chemistry     : '10XV3',
    expected_cells: 1000
]
```

Every operation is annotated with its in/out shapes.

## 4. Configuration layering

| File                     | Contains                                          |
| ------------------------ | ------------------------------------------------- |
| `params/default.yaml`    | run-variable scientific params                    |
| `params/demo.yaml`       | overrides for the bundled PBMC data               |
| `nextflow.config`        | `outputDir`, publish mode, profiles, plugins, manifest |
| `conf/base.config`       | resource labels                                   |
| `conf/modules.config`    | `ext.args`, `ext.prefix`, chemistry mapping, `scratch` |
| `conf/test.config`       | remote nf-core test data                          |
| `conf/demo.config`       | local downsampled PBMC data                       |
| `conf/aws.config`        | Batch region, queue, S3 work dir — no account data |
| `conf/aws.local.config`  | gitignored, environment-specific                  |

A value appears in exactly one file.

## 5. Publishing

`workflow.output.mode = 'copy'` and `outputDir` are set globally in
`nextflow.config`. Each output entry supplies only its relative path.

```groovy
workflow {
    main:
    // ...

    publish:
    seqkit_stats             = ch_seqkit_stats
    fastp_reports            = ch_fastp_json
    simpleaf_quant           = ch_quant_dir_per_sample
    qcatch                   = ch_qcatch_per_sample
    scanpy_cell_qc           = ch_cell_qc_per_sample
    scanpy_detect_doublets   = ch_doublets_per_sample
    scanpy_normalise_cluster = ch_cluster_per_sample
    scanpy_marker_genes      = ch_markers_per_sample
    multiqc                  = ch_multiqc_report
    final_results            = ch_headline_artifacts
}

output {
    seqkit_stats   { path { meta, tsv  -> "seqkit_stats/${meta.id}" } }
    fastp_reports  { path { meta, json -> "fastp/${meta.id}" }
                     index { path 'fastp.csv' ; header true } }
    simpleaf_quant { path { meta, dir  -> "simpleaf_quant/${meta.id}" } }
    qcatch         { path { meta, dir  -> "qcatch/${meta.id}" }
                     index { path 'qcatch.csv' ; header true } }
    // ... one entry per publish: name ...
    final_results  { path { file -> "final_results" } }
}
```

`final_results` is fed by a `mix()` of the headline channels, flattened to
bare file paths so a single simple path closure covers all of them. This
is the most readable of the options considered; the alternative — a second
`>>` publish statement inside each per-sample `path` block — spreads the
definition of "headline artifact" across eight places.

> Verify that mixing heterogeneous channel shapes into a single output
> entry is accepted. If it is not, fall back to per-output `>>` statements
> and record the change in this section.

## 6. `main.nf` structure

```
1. include statements
2. workflow { main: ... publish: ... }
3. output { ... }
```

No `nextflow.enable.dsl = 2`. No processes defined in `main.nf`.

## 7. Scientific defaults and their justification

| Param               | Default          | Why                                                   |
| ------------------- | ---------------- | ----------------------------------------------------- |
| `mapper_backend`    | `piscem`         | smaller index, faster than salmon                     |
| `umi_resolution`    | `cr-like`        | nf-core/scrnaseq default; simple and fast             |
| `count_layer`       | `S+A`            | USA convention for single-cell, not nuclei            |
| cell calling        | QCatch EmptyDrops| depth-independent; a knee discards low-RNA cells      |
| `mito_mad`          | 5                | outlier-based rather than an arbitrary percentage     |
| `mito_max_pct`      | 20               | hard ceiling; PBMC-appropriate                        |
| `min_genes`         | 200              | conventional debris floor                             |
| `n_hvg`             | 2000             | Scanpy/Seurat convention                              |
| `hvg_flavor`        | `seurat_v3`      | operates on raw counts; best for few samples          |
| `n_pcs`             | 50               | comfortably above the elbow for PBMCs                 |
| `n_neighbors`       | 15               | Scanpy default                                        |
| `leiden_resolution` | 1.0              | yields roughly 8-12 clusters on PBMCs                 |
| `marker_method`     | `wilcoxon`       | rank-based and robust                                 |
| `seed`              | 42               |                                                       |

## 8. Container design

`FROM mambaorg/micromamba:<pinned>`. Bioconda pinning gives an auditable
software manifest, which matters more here than saving 200 MB. A single
`micromamba install` layer with every package pinned `=version=build`,
then `micromamba clean --all --yes`. Multi-stage to drop the solver cache.
Non-root user.

## 9. Error handling

- Retry on exit codes 130, 137, 140, 143 with memory escalation via
  `task.attempt`, max 2 attempts.
- Otherwise `errorStrategy 'finish'`, so partial results still publish.
- Scanpy scripts exit non-zero with a readable message when a sample falls
  below `--min_cells_per_sample`.

## 9a. Downloads

Five things are fetched:

| # | Item | Fetched by | Note |
| - | ---- | ---------- | ---- |
| 1 | GRCh38 2024-A reference | `bin/download_reference.sh` | large; supplies genome.fa + genes.gtf.gz |
| 2 | splici index | built from 1, once | ~2 GB; reuse via `--simpleaf_index` |
| 3 | pbmc_1k_v3 FASTQs | `bin/download_and_downsample_testdata.sh` | ~66.6M read pairs |
| 4 | pbmc_10k_v3 FASTQs + filtered barcode matrix | `bin/download_and_downsample_testdata.sh` | ~638.9M read pairs; stream-filter |
| 5 | 10x barcode whitelist | simpleaf, at run time | silent network dependency; see R5.6 |

Item 5 is the one that surprises people. `--unfiltered-pl` without an
explicit path makes simpleaf resolve and cache the whitelist under
`ALEVIN_FRY_HOME`. Pre-seed it in the image or pass it explicitly.

The pbmc FASTQ URLs follow the pattern
`https://cf.10xgenomics.com/samples/cell-exp/3.0.0/<name>/<name>_fastqs.tar`
(a plain tar, not tar.gz, extracted with `--strip-components=1`). That
pattern is CONFIRMED for `pbmc_1k_v3` and INFERRED for `pbmc_10k_v3` and
for the filtered barcode matrix; verify both before use.

## 10. Test data design

Two distinct data paths, described in "Requirements R15.3".

The `test` profile exists so that a real, non-stub run is always cheap
enough to do after every module. It uses remote nf-core test data and a
tiny reference, so no local data management is needed and the index builds
in seconds. Its clustering output is meaningless and must never be shown.

The `demo` profile uses `pbmc_1k_v3` and `pbmc_10k_v3` from the same
donor, downsampled by two different strategies for the reason given in
"Requirements R15.5". Downloading happens once; downsampling is a separate
idempotent step over the downloaded files.

## 11. AWS readiness

- Work directory must be `s3://`. No symlink publishing anywhere.
- No `file()` on relative paths.
- `conf/aws.config` exposes region, queue and work dir as params with no
  values; real values go in the gitignored `conf/aws.local.config`.
- Container in a registry reachable from the VPC.
- `SIMPLEAF_INDEX` needs `scratch = true`; see "Requirements R4.5".
- Verify the current recommended Batch execution path
  (`aws.batch.cliPath` versus Fusion/Wave) before writing
  `conf/aws.config`.
