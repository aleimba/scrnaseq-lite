# Requirements — scrnaseq-lite

## 1. Purpose and scope

### 1.1 Purpose

Take droplet-based scRNA-seq FASTQs to an annotated, clustered AnnData
object plus three reports: a MultiQC QC report, per-sample QCatch reports,
and a Quarto analysis report.

### 1.2 In scope

Read QC -> splici quantification -> cell calling -> cell-level QC ->
doublet detection -> normalisation -> HVG -> PCA -> neighbours -> UMAP ->
Leiden clustering -> marker genes -> reporting.

### 1.3 Currently explicitly out of scope

Spatial transcriptomics; multi-modal assays (CITE-seq, ATAC, V(D)J);
multiplexed libraries (CMO/CellPlex, OCM, FFPE probe barcodes); batch
integration; trajectory inference; RNA velocity; reference-based cell-type
annotation; differential abundance testing; differential gene expression
between conditions.

## 2. Functional requirements

### R1 — Input handling

- R1.1 `--input` is validated against `assets/schema_input.json` using
  `nf-schema` before any compute-heavy task runs.
- R1.2 Samplesheet columns: `sample`, `fastq_1`, `fastq_2`, `chemistry`
  (enum `10XV2` | `10XV3` | `10XV4`), `expected_cells` (integer, optional).
- R1.2a `expected_cells` is METADATA ONLY. It must NOT be passed to
  `simpleaf quant`. simpleaf's permit-list modes `--expect-cells`,
  `--explicit-pl`, `--forced-cells`, `--knee` and `--unfiltered-pl` are
  mutually exclusive, and this pipeline uses `--unfiltered-pl` (R5.2).
  `expected_cells` is carried in `meta` for reporting and for sanity
  comparison against the number of cells QCatch actually calls.
- R1.3 A missing or unreadable FASTQ fails the run, naming the offending
  row.
- R1.4 Duplicate `sample` values fail the run. Technical replicates must
  be merged upstream; document this limitation.
- R1.5 `SEQKIT_STATS` runs on every input FASTQ. The run fails if the
  observed R1 modal length does not match the declared chemistry:
  26 bp for `10XV2`; 28 bp for `10XV3` and `10XV4`.
  *VERIFIED at T2.1 against `simpleaf chemistry lookup` in the pinned
  container, which is also a CORRECTION. The registered geometries are
  `10xv2` -> `1{b[16]u[10]x:}2{r:}` and `10xv3` / `10xv4-3p` ->
  `1{b[16]u[12]x:}2{r:}`. The **v2 UMI is 10 bp, not 12**: 16 + 10 = 26.
  The lengths above were right, the "16 bp barcode + 12 bp UMI" gloss
  repeated in "CLAUDE.md section 4.1" was not, and it is corrected there.
  Nothing in this pipeline writes a geometry string — simpleaf derives it
  from the registered chemistry name — so the T2.3 assertion is unchanged.*

### R2 — Run identity and configuration

- R2.1 All run-variable parameters live in `params/default.yaml` as pure
  YAML with no embedded expressions.
  *Extended at T1.2 and T1.5: `params/default.yaml` is NOT loaded
  automatically, so every run must pass `-params-file <file>`. Each
  profile has its OWN params file, and each is a COMPLETE parameter set
  rather than a sparse overlay: `params/default.yaml`,
  `params/demo.yaml`, `params/test.yaml`. Their key sets must stay
  identical. Correspondingly, `conf/test.config` and `conf/demo.config`
  carry ENGINE settings only and set no `params.*` at all.*
  *Verified by running at T1.5: `-params-file` OUTRANKS params defined
  in any config file, and an explicit `null` in the params file still
  wins. Had a profile config set `params.input`, running it with a
  different params file would have silently nulled the input.*
- R2.2 `nextflow.config` holds only hard-coded settings, profiles,
  plugins, manifest and environment.
- R2.3 Outputs go to `outputDir`, which is `results/<run>` when `--run` is
  given and `results/run_<yyyy-MM-dd_HH-mm-ss>` otherwise.
- R2.4 `-resume` is only supported together with `--run <name>`, because
  the timestamp is re-evaluated at every launch. `README.md` and
  `docs/usage.md` must state this.
- R2.5 Resolved parameters are written to
  `<outputDir>/pipeline_info/params_resolved.yaml`, and a run manifest
  (git commit, revision, command line, profile, container image, Nextflow
  version) to `<outputDir>/pipeline_info/run_manifest.yaml`. YAML, not
  plain text.

### R3 — Read QC

- R3.1 fastp runs once per sample, paired-end, in report-only mode.
- R3.2 fastp must not trim, quality-filter or length-filter R1.
  *Rationale: R1 carries the cell barcode and UMI at fixed offsets.*
- R3.3 fastp read outputs are not consumed downstream; `SIMPLEAF_QUANT`
  reads the original FASTQs. The pipeline performs no read cleaning.
  Document this and its consequence (no poly-G trimming).
  *Settled at T2.2: the module's `discard_trimmed_pass` input makes fastp
  write NO read files at all, which is stronger than not consuming them —
  nothing is produced to consume. The module's meta.yml describes it as
  exactly this use ("to use fastp for the output report only"), so T4.1
  passes `true`. Verified by running, because the mechanism looks like a
  bug: the module emits `def out_fq1 = discard_trimmed_pass ?: "--out1 ..."`,
  so a `true` puts a bare `true` token on the fastp command line. fastp
  1.3.6 tolerates the stray tokens, exits 0 and writes both the JSON and
  the HTML.*
- R3.4 fastp emits JSON and HTML for MultiQC.

### R4 — Reference preparation

- R4.1 A pre-built splici index is accepted via `--simpleaf_index`.
- R4.2 When `--build_index` is set, `simpleaf index` builds the splici
  (spliced + intronic) index in one command:
  `simpleaf index --output <dir> --fasta <fa> --gtf <gtf> --threads <n>
  --rlen <r2_read_length>`. Passing a genome FASTA and a GTF is what
  selects the splici strategy; there is no separate splici flag.
  *CORRECTED at T2.2, read from `simpleaf index --help` in the pinned
  container: **simpleaf 0.30.0 has no `--use-piscem` flag, and no mapper
  flag at all.** piscem is the only backend; the help has a "Piscem Index
  Options" section and no salmon option. This requirement previously
  showed `--use-piscem` in the command line. `--ref-type` defaults to
  `spliced+intronic`, so splici needs no flag either.*
  *`--rlen` IS ours to pass: it sets the flank roers adds to intronic
  sequence and must match the cDNA read length, and the module passes only
  `--threads`, the sequence inputs and `-o`. It comes from
  `params.r2_read_length` via `ext.args`. Verified by building a small
  index with `--rlen 90`.*
  *`--ram-limit-gib` is also passed. It caps SSHash's external minimizer
  sort and defaults to 8 GiB — more than the whole `test` profile
  allowance — so `conf/modules.config` derives it from `task.memory`.*
- R4.2a There is no mapper flag at quant time either, and none is needed:
  piscem is the only backend at this version.
  *AMENDED at T2.2. `params.mapper_backend` was REMOVED from all three
  params files at the same time: it existed only to carry `--use-piscem`,
  and a parameter that reaches no command line is a false record of what
  ran.*
- R4.2b Index path convention: an index built by `simpleaf index` is
  referenced as `<dir>/index/`. (An index built by `piscem build`
  directly would be `<dir>/index/piscem_idx`, which is not our case.)
- R4.3 When `--save_reference` is set, the built index is published so it
  can be reused via `--simpleaf_index` on later runs.
- R4.4 The index build must not re-run on `-resume`.
- R4.5 `SIMPLEAF_INDEX` sets `scratch = true` in `conf/modules.config`.
  The open-file limit must be at least 2048 (`ulimit -n 2048`).
  *Rationale: piscem opens very many intermediate files while indexing.
  Both nf-core and the simpleaf tutorial call this out. A low limit
  produces hangs and opaque errors rather than a clear message.*

  *Settled at T1.5 by reading the module, superseding an earlier and
  incorrect amendment that put the limit in `beforeScript`: the nf-core
  `simpleaf/index` module ALREADY does all of this in its own script
  block, and `simpleaf/quant` does the same minus the `ulimit`:*

  ```bash
  export ALEVIN_FRY_HOME=.
  ulimit -n 2048
  simpleaf set-paths
  ```

  *`ALEVIN_FRY_HOME=.` is the task working directory, so it behaves
  identically on a laptop and on AWS Batch with an S3 work directory.
  Nothing is therefore needed in the T1.4 image, in `beforeScript`, or in
  the `env` scope — only `scratch = true`. T2.1 must confirm this against
  the version actually installed before T2.2 relies on it.*
- R4.6 The reference is the 10x Human GRCh38 2024-A package
  (`refdata-gex-GRCh38-2024-A.tar.gz`), supplying `fasta/genome.fa` and
  `genes/genes.gtf.gz`. simpleaf accepts the gzipped GTF directly.
  `bin/download_reference.sh` fetches it.
- R4.7 The resulting index is roughly 2 GB and piscem maps within about
  3 GB of RAM. Resource labels must not over-provision `SIMPLEAF_QUANT`.

### R5 — Quantification

- R5.1 `simpleaf quant` with the piscem backend.
  *AMENDED at T2.2: this said "(`--mapper_backend`)". That param no longer
  exists — piscem is simpleaf 0.30.0's only backend, selected by nothing.
  See R4.2.*
- R5.2 The permit-list mode is `--unfiltered-pl`, producing a matrix with
  no cell calling applied. No knee, `--expect-cells`, `--forced-cells` or
  fixed-UMI filtering at this step: cell calling is QCatch's job (R6).
  These modes are mutually exclusive in simpleaf, so passing more than one
  is a runtime error.
  *CORRECTED at T1.4-follow-up. This requirement previously read
  "an unfiltered matrix over the full 10x whitelist", which is wrong and
  was repeated in "CLAUDE.md section 4.2" and the design dataflow. Read
  from the alevin-fry `generate-permit-list` docs: under `--unfiltered-pl`
  the whitelist is the barcode-CORRECTION reference, not the row set. Only
  barcodes observed at least `--min-reads` times become matrix rows;
  sub-threshold barcodes are not discarded outright, but their reads are
  corrected onto retained cells rather than forming rows of their own.
  simpleaf passes `--min-reads 10` explicitly. "Unfiltered" is accurate
  only in the sense that no cell calling has happened.*
  *Consequence for T2.2: `--min-reads` is a real knob that decides how
  many barcodes reach the matrix, and it is the most direct lever on the
  R15.3 fixture problem — lowering it admits more barcodes. Decide
  deliberately whether to expose it or leave simpleaf's 10 alone; do not
  leave it unconsidered.*

  *DECIDED at T2.2: left at simpleaf's default of 10, and NOT exposed as a
  param. Confirmed present as `--min-reads <MIN_READS> [default: 10]` in
  `simpleaf quant --help` at 0.30.0. Reason for deciding rather than
  deferring: on real data the default is the documented, comparable
  setting, and the only known reason to move it is to rescue a fixture too
  sparse to call cells on. That case is T2.5's, which has to measure the
  barcode count first; a param added now would invite tuning it on data
  where it should not be tuned. T2.5 may set it through `ext.args` with the
  measurement in hand.*
- R5.3 `--resolution` is a required simpleaf argument. It defaults to
  `cr-like` (`params.umi_resolution`), overridable to `cr-like-em`.
  *Found at T2.1 in the simpleaf 0.30.0 CLI that `conf/containers.config`
  now selects (see R4.5): the option `--small-thresh <N>` resolves cells
  below N with `cr-like` winner-take-all semantics REGARDLESS of
  `--resolution`. It is inert at our `cr-like` default, but it would apply
  silently to small cells under `cr-like-em`. We leave it unset, so
  alevin-fry's own default applies. Anyone switching `umi_resolution` to
  `cr-like-em` should state the effective threshold in the report rather
  than implying every cell used the EM resolution. The full set of
  resolution modes at that version is `cr-like`, `cr-like-em`,
  `parsimony`, `parsimony-em`, `parsimony-gene`, `parsimony-gene-em`.*
- R5.3a `--anndata-out` is passed, producing `af_quant/alevin/quants.h5ad`.
  *PATH CORRECTED at T2.2, by running a real quant on a small synthetic
  fixture: the file lands in `af_quant/alevin/`, NOT directly in
  `af_quant/`. This requirement, "CLAUDE.md section 4.2" and the design
  dataflow all said `af_quant/quants.h5ad`. The full `af_quant/alevin/`
  contents are `quants.h5ad`, `quants_mat.mtx`, `quants_mat_rows.txt` and
  `quants_mat_cols.txt`. T3.1's reader and T4.2's publishing rule must use
  the corrected path.*
  *Also written by simpleaf, one level up: `af_quant/gene_id_to_name.tsv`.
  That is exactly the mapping R6.4a says QCatch otherwise fetches over the
  network, and `QCATCH` receives this whole directory. T2.5 should check
  whether QCatch picks it up on its own before treating R6.4a as open.*
  QCatch accepts mtx input but only the h5ad path carries the doublet and
  mitochondrial metadata; mtx input silently loses it.
  *Settled at T1.5: the nf-core `simpleaf/quant` module HARD-CODES
  `--anndata-out` in its script block, so `conf/modules.config` must not
  pass it again. `QCATCH` likewise hard-codes `--save_filtered_h5ad` and
  `--export_summary_table`, leaving only `--skip_umap_tsne` for our
  `ext.args` (R6.2).*
  *CONFIRMED verbatim at T2.1 by reading the installed modules.*
- R5.4 The downstream count layer is selected by `--count_layer`
  (default `S+A`, the USA convention for single-cell; `S+A+U` for nuclei).
- R5.4a simpleaf writes `X` as the SUM of U + S + A and stores each as a
  separate named layer. `count_layer` must be reconstructed from the named
  layers. Reading `X` and assuming it equals `count_layer` is a silent
  correctness bug.
  *CONFIRMED at T2.2 on a real `quants.h5ad` from the pinned container,
  which is what T5.3 was told to check: `X.sum()` equalled the sum of the
  three layers exactly. The layer names are `spliced`, `unspliced` and
  `ambiguous` — note they are the full words, not U/S/A. `var` carries
  `gene_id` and `gene_symbol`; `obs` carries simpleaf's own per-cell
  metrics (`corrected_reads`, `mapped_reads`, `deduplicated_reads`,
  `mapping_rate`, `dedup_rate`, `mean_by_max`, `num_genes_expressed`,
  `num_genes_over_mean`) plus `barcodes`; `uns` carries `collate_info`,
  `gpl_info`, `quant_info` and `simpleaf_map_info`.*
  *The T1.4 anndata warning REPRODUCED here on genuine simpleaf output:
  `adata.layers.keys()` came back as
  `['ambiguous', 'spliced', 'unspliced', None]`. The spurious `None` is
  real and must be filtered before iterating (T3.1).*
- R5.5 Quantifier output is named `*_raw_matrix.h5ad`.
- R5.6 `--unfiltered-pl` with no explicit path makes simpleaf resolve and
  cache the chemistry's barcode whitelist under `ALEVIN_FRY_HOME` at run
  time. This is a network dependency. Either pre-seed the cache in the
  container image or pass the whitelist path explicitly, and document
  which was chosen. Untreated, this fails on offline nodes and in
  restricted AWS VPCs.
  *Settled at T1.5: the nf-core `simpleaf/quant` module never relies on
  the cache. Its `cellFilteringArgs` helper emits
  `--unfiltered-pl ${cb_list}` from a staged whitelist channel, so the
  whitelist is ALWAYS an explicit file.*

  *SETTLED AT T1.4 — this requirement is now CLOSED. The whitelists are
  VENDORED in this repository at
  `assets/whitelist/10x_V{2,3,4}_barcode_whitelist.txt.gz`, taken from
  nf-core/scrnaseq at tag `4.2.0`. No run reaches the network for a
  whitelist, on any node or in any VPC. Verified by download: 737,280 /
  6,794,880 / 7,372,800 barcodes, byte counts identical to those the
  GitHub contents API reports for that tag, every line a bare 16 bp
  `[ACGT]` string. Provenance and sha256 in `assets/whitelist/README.md`.*

  *V1 is deliberately absent: R1.1's `chemistry` enum is `10XV2`, `10XV3`,
  `10XV4`, so a V1 list could never be selected.*

  *Cost correction: this requirement previously said "roughly 20 MB".
  2.2 + 18.4 + 13.4 MB is **34.0 MB**. The earlier figure was an
  arithmetic slip.*

  *CLOSED AT T2.2: **a GZIPPED permit list works.** No decompression step
  is needed and none may be added.*
  *Proven twice, not inferred. From source, at the version that runs:
  alevin-fry v0.18.1 `src/cellfilter.rs:1914` loads the unfiltered permit
  list through `niffler::from_path`, and `Cargo.toml:55` pins
  `niffler = { version = "3.0.0", features = ["gz"] }`, so the compression
  is auto-detected. (Two other readers exist —
  `cellfilter.rs:850` and `utils.rs:1074 read_filter_list` — that use a
  plain `BufReader`, but they serve the multiplex and `--explicit-pl`
  paths, not `--unfiltered-pl`.) Then by running: a small synthetic
  fixture quantified twice in the pinned container, once with
  `assets/whitelist/10x_V2_barcode_whitelist.txt.gz` and once with the
  same file decompressed, produced 3 barcode rows each — the exact three
  synthesised barcodes — and BYTE-IDENTICAL `quants_mat.mtx` files. A
  binary read of the gzip would have produced a junk whitelist and lost
  those barcodes, so this distinguishes "accepted" from "silently
  wrong".*
- R5.7 simpleaf requires `ALEVIN_FRY_HOME` to be set and
  `simpleaf set-paths` to have been run.
  *Settled at T1.5: the nf-core simpleaf modules do both in their own
  script blocks, with `ALEVIN_FRY_HOME=.` (task-local, so identical on a
  laptop and on AWS Batch). Nothing is needed in the container. See
  R4.5.*

  *CONFIRMED at T2.1 against the INSTALLED modules (simpleaf 0.25.0
  container `community.wave.seqera.io/library/simpleaf:0.25.0--b9f96d8b71a01864`),
  with one correction: `simpleaf/index` does all three lines;
  `simpleaf/quant` does `export ALEVIN_FRY_HOME=.` and
  `simpleaf set-paths` but NO `ulimit`. Only `SIMPLEAF_INDEX` needs
  `scratch = true`.*

  *ALSO FOUND AT T2.1, and it is not a configuration issue: the bioconda
  `alevin-fry 0.15.0 hd612981_0` build inside that container aborts with
  SIGILL on a CPU lacking the instruction set it was compiled for
  (reproduced on Zen 2: avx, avx2, bmi2, no AVX-512). `simpleaf set-paths`
  shells out to `alevin-fry --version`, so it exits 1 and takes the whole
  task with it under `bash -ue`. simpleaf and piscem from the same image
  are fine, and alevin-fry 0.11.2 from the older
  `quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0` is fine.*

  *FIXED at T2.1, in `conf/containers.config`. Root cause: alevin-fry
  `v0.15.0` ships `.cargo/config.toml` with
  `rustflags = ["-C", "target-cpu=native"]` and the bioconda recipe adds
  no RUSTFLAGS of its own, so the published binary is tuned to the build
  worker's CPU. Upstream dropped it in commit `0e13d52` (2026-07-10) for
  precisely this reason, in favour of the portable
  `target-cpu=x86-64-v3 -C target-feature=+avx2`; alevin-fry 0.18.1
  carries the fix. `conf/containers.config` therefore overrides
  `SIMPLEAF_INDEX|SIMPLEAF_QUANT` onto
  `quay.io/biocontainers/simpleaf:0.30.0--hd612981_0` (simpleaf 0.30.0,
  alevin-fry 0.18.1, piscem 0.23.0), whose recipe requires
  `alevin-fry >=0.18.1`. Verified by running on the host that reproduces
  the crash: `simpleaf set-paths` exits 0, and a process named
  `SIMPLEAF_QUANT` completes and reports `alevin-fry 0.18.1` to the
  versions topic. The vendored modules are NOT edited. The override is a
  temporary deviation: delete the file once `nf-core modules update`
  brings a simpleaf module pinning alevin-fry 0.18.1 or newer.*

### R6 — Cell calling

- R6.1 `QCATCH` calls cells from the quantifier output using its two-step
  EmptyDrops-style ambient model.
- R6.2 QCatch runs with `--skip_umap_tsne --export_summary_table
  --save_filtered_h5ad`.
- R6.3 QCatch's doublet flags are not used.
  *Rationale: its Scrublet run exposes no seed; see R7.5.*
- R6.4 QCatch chemistry is derived from the samplesheet `chemistry` column
  through the single mapping shared with the simpleaf chemistry, so the
  two cannot diverge. See "CLAUDE.md section 4.3: Chemistry mapping — one
  source of truth". `--qcatch_n_partitions` is the override for custom
  assays with no mapping.
  *VERIFIED at T2.1 by running `qcatch --help` in the pinned container
  (qcatch 0.2.12): `--chemistry` accepts `10X_3p_v2`, `10X_3p_v3`,
  `10X_3p_v4`, `10X_5p_v3`, `10X_3p_LT` and `10X_HT`, so all three of our
  mappings are real. `params.qcatch_n_partitions` maps to the tool's
  `--n_partitions`, documented as "use only when working with a custom or
  unsupported chemistry"; it overrides the chemistry-based configuration.*
- R6.4a QCatch fetches its gene-id-to-name mapping from a REMOTE registry
  unless `--gene_id2name_file` is supplied, and when that lookup fails it
  omits the mitochondrial plots instead of failing.
  *Found at T2.1 in `qcatch --help`. Same class of hidden run-time network
  dependency the vendored whitelists closed for simpleaf at T1.4, and the
  same nodes are affected: offline hosts and restricted AWS VPCs. T2.5
  decides whether to derive the TSV (`gene_id`, `gene_name`, no header)
  from the GTF the index was built from.*

  *NARROWED at T2.2, and it may already be solved. simpleaf writes
  `gene_id_to_name.tsv` itself, into BOTH the index directory
  (`<idx>/index/` and `<idx>/ref/`) and the quant output
  (`af_quant/gene_id_to_name.tsv`), and `QCATCH` is handed that whole
  `af_quant` directory. The h5ad's `var` also already carries
  `gene_symbol`. T2.5 must check whether QCatch finds the file on its own
  before building anything. If it does not, note the constraint: the
  option CANNOT be supplied through `ext.args`, because the module has no
  input for it and an unstaged path is not visible inside the container.
  It would need the file staged some other way, or the limitation
  documented.*
- R6.5 The QCatch HTML report is published per sample.
- R6.6 QCatch's cell-called matrix is named `*_filtered_matrix.h5ad`.
  *AMENDED at T2.1: the vendored module does not produce that name. It
  renames its outputs to `${prefix}_qcatch_report.html`,
  `${prefix}_filtered_quants.h5ad` and `${prefix}_metrics_summary.csv`,
  and `SIMPLEAF_QUANT` emits a whole `af_quant` directory containing
  `quants.h5ad` rather than an `*_raw_matrix.h5ad` (R5.5). Vendored
  modules are never edited, so both provenance suffixes are applied where
  the files are PUBLISHED, in the T4.2 `output` block. The requirement
  stands; only the place it is satisfied moves.*
- R6.7 `--cell_calling` selects how cells are separated from empty droplets.
  Enum, default `qcatch`:
  - `qcatch` — R6.1's two-step EmptyDrops-style ambient model. The only
    acceptable choice for real data.
  - `threshold` — QCATCH does not run. The unfiltered quantifier matrix
    goes straight to `SCANPY_CELL_QC`, where `min_genes` and
    `min_cells_per_gene` (R7.2) act as a crude cut. NO ambient model.

  *Added at T1.4-follow-up. Rationale: the `test` fixture is ~60,000 read
  pairs against 5,000 declared cells, about 12 reads per cell, and
  nf-core/scrnaseq sets `skip_qcatch = true` on that exact dataset with the
  comment "module does not work on small dataset". Without an escape hatch
  the `test` profile cannot exercise the DAG past quantification (R15.3).*

  *NAMING, deliberate: nf-core calls this `skip_qcatch`, but in THEIR
  pipeline QCatch is an optional QC report and the raw matrix is emitted
  either way. Here QCATCH is on the critical path and is the cell caller,
  so "skip the QC step" would badly understate the consequence. The
  parameter is named for the scientific choice it makes.*

  *CONSTRAINTS on `threshold`, all mandatory:*
  - *It is not scientifically acceptable for real data. A global cut
    discards genuine low-count cells; see "CLAUDE.md section 4.4".*
  - *Only `qcatch` may produce `*_filtered_matrix.h5ad` (R6.6). Output
    from the `threshold` path must never carry that name.*
  - *The method is recorded in `adata.uns['cell_calling_method']`, in
    `run_manifest.yaml` (R2.5), and as a visible banner in the MultiQC and
    Quarto reports whenever it is not `qcatch`.*
  - *An unrecognised value is a hard failure, not a silent default.*

  *RESIDUAL RISK for T2.5 — the mechanism is CONFIRMED, the outcome on
  this fixture is not. Confirmed from the alevin-fry
  `generate-permit-list` docs: under `--unfiltered-pl` only barcodes seen
  at least `--min-reads` times become matrix rows, and simpleaf passes
  `--min-reads 10`. So the matrix handed to the `threshold` path is
  already thinned, and on a ~60,000-read-pair fixture it may hold too few
  barcodes for that path to help either. NOT YET MEASURED: count the rows
  in the fixture's raw matrix at T2.5 before declaring the escape hatch
  works. If it is near-empty, lowering `--min-reads` via `ext.args` is the
  next lever (see R5.2), not a weakening of the pipeline.*

### R7 — Cell QC and doublets

- R7.1 Per cell: total counts, detected genes, % mitochondrial,
  % ribosomal, % haemoglobin.
- R7.2 Filtering thresholds are MAD-based (`--mito_mad`, `--count_mad`)
  with absolute floors and ceilings from params. No hard-coded constants.
- R7.3 Pre- and post-filter violin plots, a barcode-rank plot, and a
  counts-vs-genes scatter coloured by % mitochondrial, as PNG.
- R7.4 A per-sample QC summary TSV formatted as MultiQC custom content.
- R7.5 `scanpy.pp.scrublet` runs per sample, seeded with `--seed`.
- R7.6 Predicted doublets are flagged in `.obs` and removed when
  `--remove_doublets` is true. Removal happens before HVG and PCA.
- R7.7 The run fails with a readable message if a sample retains fewer
  than `--min_cells_per_sample` cells.

### R8 — Normalisation and feature selection

- R8.1 Raw counts preserved in `adata.layers['counts']`.
- R8.2 Shifted-log normalisation: `normalize_total` (`--target_sum`,
  null = median library size) then `log1p`.
- R8.3 HVG selection with `--n_hvg` and `flavor='seurat_v3'` on raw counts.

### R9 — Dimensionality reduction and clustering

- R9.1 PCA (`--n_pcs`), neighbours (`--n_neighbors`).
- R9.2 UMAP (`--umap_min_dist`) and Leiden (`--leiden_resolution`,
  `flavor='igraph'`, `n_iterations=2`).
- R9.3 All stochastic steps use `--seed`.
- R9.4 UMAPs coloured by cluster, sample, total counts and % mitochondrial,
  as PNG.

### R10 — Marker genes

- R10.1 `rank_genes_groups` with `--marker_method` (default `wilcoxon`).
- R10.2 Results as TSV: cluster, gene, logFC, p-value, adjusted p-value.
- R10.3 A dotplot of the top `--n_markers_plot` genes per cluster.
- R10.4 Nothing in this section may be named or described as differential
  expression. It is marker detection.

### R11 — Reporting

- R11.1 MultiQC aggregates seqkit stats, fastp, and custom content from
  the QCatch summary CSV and the Scanpy QC TSV. MultiQC is QC-only.
- R11.2 QCatch HTML reports are published, one per sample.
- R11.3 A Quarto HTML analysis report carries the UMAPs, cluster sizes and
  the marker table.
- R11.4 Every table in a report is also written as a standalone TSV; every
  figure as PNG; the final AnnData as `.h5ad`.

### R12 — Reproducibility

- R12.1 `versions.tsv` (`process`, `tool`, `version`; tab-delimited) in
  the output directory, collected from every process. `process` may repeat
  when one process uses several tools.
- R12.1a Newer nf-core modules have begun migrating version reporting from
  a `versions.yml` output file to the topic-channel convention. Determine
  which convention the ACTUALLY INSTALLED modules use, and make local
  modules match it.

  *Settled at T1.5 by reading the modules: the two conventions ARE mixed
  in nf-core/scrnaseq 4.2.0, so "do not mix them" is not achievable for
  the vendored half of this pipeline.*
  - *`SIMPLEAF_INDEX` and `SIMPLEAF_QUANT` emit `path "versions.yml"`.*
  - *`QCATCH` emits `tuple val("${task.process}"), val('qcatch'),
    eval("qcatch --version ..."), emit: versions_qcatch, topic: versions`.*

  *`COLLECT_VERSIONS` (T3.5) must therefore consume BOTH: the
  `versions.yml` files and the `versions` topic. Our own local modules
  still pick one convention and use it consistently; T2.1 decides which,
  after checking what the installed modules actually do.*

  ***SETTLED AT T2.1, and the T1.5 note above is SUPERSEDED — it described
  the old module snapshot pinned by nf-core/scrnaseq 4.2.0, not what
  `nf-core modules install` fetches. The conventions are NOT mixed. All
  six installed modules — `fastp`, `seqkit/stats`, `multiqc`,
  `simpleaf/index`, `simpleaf/quant`, `qcatch` — report versions as
  `tuple val("${task.process}"), val('<tool>'), eval("<cmd>")` into
  `topic: versions`. No `versions.yml` file exists anywhere under
  `modules/nf-core/`. `COLLECT_VERSIONS` reads the topic ONLY, and every
  local module uses the topic convention.***
  - *ONE EXCEPTION, deliberate and load-bearing: `MULTIQC` emits a plain
    `emit: versions` channel and stays OUT of the topic, because it
    CONSUMES the topic and pushing into it "will let the pipeline hang
    forever" — the module carries that comment itself. So `MULTIQC` runs
    after `COLLECT_VERSIONS`, and its version reaches `versions.tsv`
    through its own channel.*
  - *A failing `eval` does NOT fail the task. Verified at T2.1 with a
    throwaway process whose eval command aborted: the topic received an
    EMPTY string and the run reported success. `COLLECT_VERSIONS` must
    reject or explicitly mark an empty version, never write a blank cell.*
  - *`qcatch --version` prints `qcatch version 0.2.12`, and the module's
    `sed -e 's/qcatch //g'` leaves `version 0.2.12`. Strip a leading
    `version ` token in `COLLECT_VERSIONS`; the module is not edited. The
    other five evals are clean.*
- R12.2 Every tool pinned to an exact version in the container.
  *Partly met at T1.4, and the shortfall is deliberate. `docker/Dockerfile`
  pins the 21 packages it NAMES to exact `version=build` triples, all taken
  from a `micromamba create --dry-run --json` solve. It does NOT pin the
  full 243-package closure, so `numba`, `llvmlite` and `libopenblas` —
  all of which affect numerical output — float on rebuild. A given built
  image is reproducible by digest; a rebuild from source months later is
  not guaranteed bit-identical. Fix, if this is ever wanted: install from
  a generated lockfile instead of a package list. See `docs/container.md`.*
- R12.3 Nextflow `timeline`, `report`, `trace` and `dag` in
  `pipeline_info/`.
- R12.4 Identical inputs, params, seed and container digest give identical
  cluster assignments. Document that PCA may vary across differing BLAS
  thread counts.
  *MEASURED at T1.4, in the built image, seeded, across separate
  processes. Two runs at the same thread count were bit-identical for
  every step. Re-running with `OMP_NUM_THREADS=1` and
  `OPENBLAS_NUM_THREADS=1`:*
  - *PCA CHANGED (`a9d8995ac8386167` -> `cebcc34aa1e99705`) and UMAP
    inherited the change (`23f19efb7c366681` -> `e59b4bd7aaeb86a3`).*
  - *scrublet, `seurat_v3` HVG selection, Leiden assignments and
    `rank_genes_groups` were UNCHANGED.*
  *So the caveat above is correct and belongs in the docs. Write it as
  "PCA and UMAP coordinates may differ across BLAS thread counts". Do NOT
  promise that cluster assignments are stable across thread counts: they
  held on well-separated synthetic data, which is not evidence about real
  samples. Full digests in `docs/container.md`.*

### R13 — Output organisation

- R13.1 Each process publishes to `<outputDir>/<process_name_lowercased>/`
  in copy mode.
- R13.2 `<outputDir>/final_results/` holds copies of the headline
  artifacts: MultiQC HTML, QCatch HTMLs, Quarto analysis report, UMAP
  PNG, marker TSV, final `.h5ad`, `versions.tsv`.
  *Rationale for copies rather than symlinks: S3 has no symlinks and AWS
  Batch is a target.*
- R13.3 Per-sample outputs carry an index file (CSV with header).

### R14 — Portability

- R14.1 Runs with `-profile docker` locally.
- R14.2 Runs with `-profile awsbatch` on AWS Batch.
- R14.3 One PIPELINE-AUTHORED container image, `FROM` a minimal
  micromamba base, target under 2 GB compressed.
  *Amended at T1.5, recording the decision taken at T1.3: this pipeline
  uses PER-MODULE BIOCONTAINERS. Each vendored nf-core module keeps the
  container it pins, and no global `process.container` override is set.
  Our image therefore covers the local modules only — the Scanpy steps,
  Quarto, and the version and manifest scripts.*
  *Rationale: overriding a vendored module's container would leave its
  pinned biocontainer tag as dead text while silently making our image
  responsible for every tool it needs. Cost accepted: several image
  pulls rather than one, and no single digest to pin for R12.4;
  `versions.tsv` (R12.1) is what records what actually ran.*
  *BUILT AND MEASURED at T1.4. Tag `scrnaseq-lite:0.1.0`, no registry
  prefix: built from `docker/Dockerfile` and used locally by
  `-profile docker`. Base `mambaorg/micromamba:2.9.0-debian13`, pinned by
  digest `sha256:ed1bc628...`. Two stages on that one base; the second
  copies `/opt/env` and sets `ENTRYPOINT []`.*
  *Size: 2.66 GB uncompressed, **628 MB gzip-compressed** — the target is
  under 2 GB compressed, so this passes comfortably.*
  *The AWS target is AMAZON ECR. That URI embeds the AWS account ID, so it
  must never appear in a tracked file (N4); it belongs in the gitignored
  `conf/aws.local.config` as a container override.*
  *Two run-time constraints the image MUST keep satisfying, both verified:
  (a) `docker.runOptions = '-u $(id -u):$(id -g)'` gives an arbitrary host
  uid with no `/etc/passwd` entry and no writable home, so `HOME`,
  `MPLCONFIGDIR`, `NUMBA_CACHE_DIR` and the `XDG_*` variables must point
  at world-writable paths or scanpy fails to IMPORT; (b) conda activation
  is skipped in favour of `PATH`, so the eight `QUARTO_*` variables that
  `etc/conda/activate.d/quarto.sh` would have exported must be set
  explicitly or `quarto render` fails. See `docs/container.md`.*
- R14.4 `-profile demo` completes on 2 CPUs / 8 GB RAM in under 20 minutes
  using the bundled downsampled data and a pre-built index.

### R15 — Testing and test data

- R15.1 Every process has a `stub` block. `-stub-run` must complete
  cleanly at every stage of development.
- R15.2 No `nf-test` and no CI in v0.1.0, by decision.
- R15.3 Two data profiles exist and serve different purposes:
  - `test` — small remote data from `nf-core/test-datasets` with the
    reference used by nf-core/scrnaseq 4.2.0 `conf/test.config`. Exercises
    the full DAG in minutes with no local data management. Results are NOT
    biologically meaningful and must not be presented as such.

    *Settled at T1.5. The data is mouse (GRCm38 chr19, gencode vM19) and
    10XV2; verified from the FASTQ stream, R1 is 26 bp and R2 is 90 bp.
    `assets/samplesheet_test.csv` is ours, because nf-core's own
    samplesheet has no `chemistry` column and lists `Sample_Y` twice for
    two lanes, which R1.4's uniqueness rule rejects by design.*

    *KNOWN CEILING: nf-core sets `skip_qcatch = true` on this dataset
    ("module does not work on small dataset") and its 4.2.0 simpleaf test
    emits no QCATCH output at all — roughly 60,000 read pairs against
    5,000 declared cells is about 12 reads per cell.*

    *AMENDED after T1.4. The earlier text concluded "this pipeline cannot
    skip cell calling" and left QCATCH-onward as an open risk. R6.7 now
    supplies the answer: `params/test.yaml` sets
    `cell_calling: "threshold"`, so QCATCH does not run on the fixture and
    the unfiltered matrix goes straight to `SCANPY_CELL_QC`. Cell calling
    is not skipped; it is done crudely, and labelled as such everywhere it
    appears. This is acceptable ONLY because this profile's results are
    already declared meaningless.*

    *So `-profile test -stub-run` is the always-green check, and
    `-profile test,docker` is expected green through `SEQKIT_STATS`,
    `FASTP`, `SIMPLEAF_INDEX`, `SIMPLEAF_QUANT` and the Scanpy chain.
    `params/test.yaml` lowers the QC and clustering values to give it the
    best chance, and marks each as a fixture value.*

    *RESIDUAL RISK, for T2.5: `--unfiltered-pl` applies alevin-fry's
    `--min-reads` cut before a barcode reaches the matrix, so the
    unfiltered matrix may itself be near-empty on a fixture this small.
    Count the barcodes. If there are too few, say so plainly and record
    what `-profile test,docker` does cover; do not weaken the pipeline to
    manufacture a green run.*
  - `demo` — local downsampled human PBMC data (R15.4). This is the run
    that produces presentable results.
- R15.4 The demo dataset is 10x `pbmc_1k_v3` and `pbmc_10k_v3`. Same
  healthy donor, same 10x 3' v3 chemistry, two independent captures, both
  `cellranger count` (singleplex). Multiplexed libraries (OCM, CMO, FFPE)
  are out of scope: demultiplexing requires `cellranger multi` and a
  barcode samplesheet, and running a multiplexed library through simpleaf
  would silently quantify all pooled samples together.

  Verified figures from the 10x metrics summary CSVs (https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_metrics_summary.csv and https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_metrics_summary.csv):

  | Metric              | pbmc_1k_v3   | pbmc_10k_v3   |
  | ------------------- | ------------ | ------------- |
  | Estimated cells     | 1,222        | 11,769        |
  | Mean reads/cell     | 54,502       | 54,286        |
  | Median genes/cell   | 1,919        | 1,906         |
  | Total read pairs    | 66,601,887   | 638,901,019   |
  | Sequencing saturation | 70.8%      | 68.2%         |
  | Median UMI/cell     | 6,628        | 6,521         |
  | Read 1 / Read 2     | 28 bp / 91 bp| 28 bp / 91 bp |

  Both are 28 bp R1 (16 bp barcode + 12 bp UMI) and 91 bp R2, so
  `r2_read_length` is 91 and one splici index serves both.

- R15.4a The two samples are technically near-identical in depth and
  genes per cell. The real contrast is cell loading: ~1,200 versus
  ~11,800 recovered cells, which on the 10x v3 loading curve means
  markedly different doublet rates (roughly 0.8% versus 8%). Barcode
  subsetting preserves that difference. The analysis report should
  highlight the doublet-rate contrast, not depth.

- R15.5 Each dataset is downloaded once and then downsampled, using
  different strategies for a stated reason:
  - `pbmc_1k_v3`: random read-pair subsample, fixed seed, fraction
    approximately 0.36, targeting ~24M read pairs
    (1,222 cells at roughly 19,600 reads per cell).
  - `pbmc_10k_v3`: subsample by CELL BARCODE to ~1,200 of the 11,769
    barcodes, then random-subsample those read pairs to ~24M
    (~1,200 cells at roughly 20,000 reads per cell).

  *Rationale: random read downsampling reduces reads per cell. Applied to
  the 10k dataset at a file size comparable to the 1k dataset it would
  leave roughly 1,200 reads per cell, too shallow for ambient-model cell
  calling. Barcode subsetting preserves per-cell depth and changes only
  cell count. The two subsampled samples end up matched on both cell count
  and depth, so per-sample differences reflect the original experiments,
  not the subsampling.*

- R15.5a `pbmc_10k_v3` is roughly 9.6x the size of `pbmc_1k_v3`
  (638.9M versus 66.6M read pairs). Read the ACTUAL tarball sizes from the
  10x dataset pages before downloading and report them to the user.
  Preferred implementation: stream the source FASTQs from the 10x CDN
  through the barcode filter and the random subsampler in a single pass,
  writing only the subset, so the full file never occupies disk. The
  alternative, downloading in full and deleting afterwards, needs tens of
  GB of transient space and must be offered only if streaming fails.

- R15.6 Both subsampling steps use a fixed seed ('42') and are fully reproducible
  from `bin/download_and_downsample_testdata.sh` alone. R1 and R2 must be
  read in lockstep; assert equal read counts and identical read-name order
  in the outputs before writing.

- R15.7 A `--simple` flag falls back to random subsampling for both
  samples, printing a warning that the 10k sample will be shallow.

- R15.8 The script prints expected download and on-disk sizes and requires
  confirmation before downloading. It is idempotent and re-runnable. Additionally,
  the script does not use curl's silent mode, but shows the progress of the downloads.

- R15.9 Both datasets are licensed CC BY 4.0. `README.md` must attribute
  each dataset by name, credit 10x Genomics, link the dataset page and
  link the CC BY 4.0 licence, following the 10x citation guidelines at
  https://www.10xgenomics.com/support/software/cell-ranger/latest/miscellaneous/cr-citations

## 3. Non-functional requirements

- N1 Pipeline-authored Nextflow under 1500 lines.
- N2 No process requires more than 16 GB RAM under `-profile demo`.
- N3 Every scientific default in `params/default.yaml` carries a one-line
  justification comment.
- N4 No private or environment-specific information in any tracked file.

## 4. Open questions

All questions from earlier drafts are resolved. Recorded here so the
answers are not re-litigated during implementation.

- **Q1 simpleaf flags — RESOLVED.** `simpleaf quant` requires exactly one
  of `--expect-cells | --explicit-pl | --forced-cells | --knee |
  --unfiltered-pl`. This pipeline uses `--unfiltered-pl` (R5.2).
  No mapper flag is passed at all: `--use-piscem` was removed from
  simpleaf, which now has piscem as its only backend (corrected at T2.2;
  see R4.2 and R4.2a). `--resolution` is required (R5.3).
- **Q2 QCatch module — RESOLVED.** An nf-core module named `qcatch`
  exists. Inspect it with `nf-core modules info qcatch`, then
  `nf-core modules install qcatch`. Do not write a local module.
- **Q3 Licence — RESOLVED.** Both datasets are CC BY 4.0; see R15.9.
- **Q4 Read lengths — RESOLVED.** Both are 28 bp R1 / 91 bp R2; see R15.4.
- **Q5 test-profile URLs — RESOLVED by the user.** Still read them from
  nf-core/scrnaseq 4.2.0 `conf/test.config` and the `nf-core/test-datasets`
  scrnaseq branch rather than reconstructing them.

If a NEW ambiguity appears during implementation, add it here with a
`[BLOCKS "<task name>"]` marker and stop.
