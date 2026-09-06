# CLAUDE.md — scrnaseq-lite

## 1. What this repo is

A minimal, scientifically sound Nextflow pipeline for droplet-based
single-cell RNA-seq. Deliberately small: core steps only, optimised for
speed and reviewability rather than feature coverage.

Currently, out of scope by design, not by oversight: spatial transcriptomics;
multi-modal assays (CITE-seq, ATAC, V(D)J); multiplexed libraries
(CMO/CellPlex, OCM, FFPE probe barcodes); batch integration; trajectory
inference; RNA velocity; reference-based cell-type annotation;
differential abundance; differential gene expression between conditions.

## 2. Environment

- Nextflow 26.04.6. Verified in this environment: the workflow `output`
  block is stable; `nextflow.enable.dsl = 2` is NOT required;
  `nextflow lint` exists.
- `nextflow.preview.output` must NOT appear anywhere. It is unsupported
  from Nextflow 25.10 onward.
- `nf-core` tools run from a pyenv virtualenv named `scrnaseq-lite`,
  activated in the working directory.
- All pipeline software lives Docker images. Nothing pipeline-related
  is ever installed on the host.
- `.claude/settings.json` is maintained by hand and is final. Do not read,
  edit, regenerate or comment on it.
- **The repository root stays clean and nf-core-shaped.** Agent-facing
  material lives under `.claude/`: `specs/` (requirements, design, tasks),
  `reference/` (target files to compare your output against, NOT to copy
  into the repo) and `skills/` (slash commands). `CLAUDE.md` is the only
  agent file at the root.
  Never create new top-level directories or scatter planning documents in
  the root; if you need a scratch document, put it in `docs/` and add it
  to `.gitignore`.

## 3. Non-negotiable Nextflow rules

### 3.1 Publishing

- **No `publishDir` in any pipeline-authored process.** Publishing is
  declared in the entry workflow's `publish:` section plus the top-level
  `output` block in `main.nf`.
- Vendored nf-core modules are exempt and are never edited.
- Every name in `publish:` has a matching entry in `output`, and vice
  versa. Nextflow enforces this; do not work around it.
- Closure parameters in an output `path` directive must match the
  published channel's shape exactly.
- Publish mode is set once, globally: `workflow.output.mode = 'copy'`.
  Never `symlink` — S3 has no symlinks and AWS Batch is a target.
- Use `index { path '<name>.csv' ; header true }` on per-sample outputs.

### 3.2 Syntax

- **Never use the implicit `it`.** Always name closure parameters:
  `.map { meta, reads -> ... }`.
- **Be explicit in names.** Process, subworkflow, channel and variable
  names must state what they contain or do: `ch_quant_dir_per_sample`,
  not `ch_q`; `QUANTIFY_WITH_SIMPLEAF`, not `QUANT`.
- Every channel operation carries a shape comment:
  `// in: [ meta, [r1, r2] ]  ->  out: [ meta, path(json) ]`
- Every process has a `stub` block that touches all declared outputs and
  emits its version information.
- **Never embed a python (or R, or awk) script in a `.nf` file — not as a
  heredoc, not as a `#!/usr/bin/env python` script block.** Every script
  is a real, executable file in `bin/`, called from a bash script block.
  A script in `bin/` can be linted, run by hand outside Nextflow, and
  diffed sensibly; a heredoc can do none of those and hides Groovy
  interpolation inside another language's syntax.
- Call those scripts **by bare name** (`my_script.py --arg ...`), as
  nf-core does. Nextflow appends the pipeline's `bin/` to `PATH` for every
  task, so this is a guarantee, not an assumption. Do NOT write
  `${projectDir}/bin/...`: that hard-codes a host path which need not exist
  inside the container on a remote executor.
- Every local process therefore has a **bash** script block, which is also
  a hard requirement: e.g. version topic uses an `eval` output, and
  Nextflow allows `eval` only on bash scripts.

### 3.3 Layout

- Process names are UPPER_SNAKE_CASE.
- **Published subdirectory name = the process name, lower-cased.**
  `SIMPLEAF_QUANT` publishes to `simpleaf_quant/`.
- Local modules: `modules/local/<lowercase_name>/main.nf`.
- Vendored: `modules/nf-core/`, installed with `nf-core modules install`,
  never modified.
- `conf/modules.config` carries `ext.args` and `ext.prefix` only.

### 3.4 Configuration split

- `params/default.yaml` — run-variable parameters. Pure YAML, no Groovy.
  It is NOT loaded automatically: `nextflow.config` does not parse it and
  nf-schema does not apply `default` values from `nextflow_schema.json`
  to unset params. **Every run must pass `-params-file <file>`.**
- Every params YAML is a COMPLETE parameter set, never a sparse overlay.
  `params/demo.yaml` and `params/test.yaml` therefore repeat the defaults
  with their own values filled in. Keep all three key sets identical.
- `-params-file` OUTRANKS params defined in any config file, and an
  explicit `null` in the params file still wins. Therefore
  `conf/test.config` and `conf/demo.config` carry ENGINE
  settings only and set no `params.*`.
- `nextflow.config` — `outputDir`, publish mode, profiles, plugins,
  manifest, env.
- `conf/base.config` — resource labels.
- `conf/modules.config` — per-process `ext.args`, `ext.prefix`, `scratch`,
  and the protocol -> chemistry mapping.
- `conf/containers.config` — container overrides for vendored nf-core
  modules, and nothing else. Keeping these out of `conf/modules.config` is the point: overrides are temporary and must stay easy to see and to remove. Overriding the container is the sanctioned move; editing a vendored module is not.
- `conf/test.config` — engine settings for the `test` profile
  (`resourceLimits`). Small remote nf-core test data. Its parameters lives
  in `params/test.yaml` and the remote data in `assets/samplesheet_test.csv`.
- `conf/demo.config` — engine settings for the `demo` profile
  (`resourceLimits`). Local downsampled human PBMC data. Its parameters lives
  in `params/demo.yaml` and the data is downloaded by a shell script in bin
  and then available in gitignored `/data`.
- A value appears in exactly one file. Never duplicated (except for the
  `params/*.yaml` files, which are complete, see above).
- Validation uses the `nf-schema` plugin only. Do NOT add `nf-validation`;
  it is the superseded predecessor and the two conflict. The schema file
  is `nextflow_schema.json` (JSON).

### 3.5 Run directory and resume

- `params.run` names the run. When null, a timestamp is used. It is
  declared in `nextflow.config` and set ONLY with `--run <name>` on the
  command line. It must never appear in a params YAML: every run sharing
  that file would then write into one output directory.
- `outputDir = params.run ? "results/${params.run}" : "results/run_<ts>"`.
- The timestamp is re-evaluated at every launch. Therefore **`-resume`
  only produces a coherent output directory when `--run <name>` is
  given.** State this in `README.md` and `docs/usage.md`.
- `outputDir` must never feed into any task input, or resume breaks.
- `-output-dir` is a single-dash Nextflow option, not a pipeline param.
  `--outputDir` does not exist. `--run` is the pipeline-level control.

## 4. Scientific constraints — do not "improve" these

### 4.1 Read handling

- **R1 is barcode+UMI**: 26 bp for 10x v2 (16 bp cell barcode + **10 bp**
  UMI); 28 bp for v3 and v4 (16 bp cell barcode + 12 bp UMI). Nothing may
  adapter-trim, quality-trim or length-filter it — that destroys the fixed
  barcode/UMI offsets. The v3/v4 geometry is `1{b[16]u[12]x:}2{r:}`; v2 is
  `1{b[16]u[10]x:}2{r:}`. The pipeline never writes a geometry string —
  simpleaf derives it from the registered chemistry name.
- **fastp is a QC-only branch.** It runs once per sample, paired-end, in
  report-only mode (`--disable_adapter_trimming
  --disable_quality_filtering --disable_length_filtering`), therefore there's
  no fastp read output. `SIMPLEAF_QUANT` reads the original FASTQs.
  Accepted consequence: no poly-G trimming is applied; poly-G
  reads simply fail to map in piscem. Do not add a cleaning pass without
  updating `.claude/specs/01-requirements.md` first.
- `SEQKIT_STATS` asserts the observed R1 modal length matches the
  chemistry declared in the samplesheet. Mismatch is a hard failure. That
  assertion is the only reason seqkit is in this pipeline; fastp covers
  every other read metric.

### 4.2 Quantification

- `simpleaf` with the **piscem** mapping backend and alevin-fry
  quantification.
- Index: `simpleaf index --output <dir> --fasta <fa> --gtf <gtf>
  --threads <n> --rlen <r2_read_length> --ram-limit-gib <n>`. Supplying a
  genome FASTA plus a GTF is what selects the splici (spliced + intronic)
  strategy; there is no separate splici flag, and `--ref-type` already
  defaults to `spliced+intronic`.
  **There is no `--use-piscem` flag.**, simpleaf v0.30.0 has no
  mapper flag at all, because piscem is its only backend. Neither index nor
  quant selects a mapper.
  `--rlen` and `--ram-limit-gib` are ours to pass, from
  `conf/modules.config`; the module passes neither, and `--ram-limit-gib`
  otherwise defaults to 8 GiB regardless of what the task was granted.
- Reference the index as `<dir>/index/`. (`<dir>/index/piscem_idx` would
  be correct only for an index built by `piscem build` directly.)
- The splici index is built once and reused via `--simpleaf_index`. There
  is no downloadable pre-built simpleaf index; 10x distributes Cell Ranger
  references only.
- `SIMPLEAF_INDEX` must set `scratch = true` in `conf/modules.config`.
  simpleaf's mapper opens very many temporary files at once and stalls or
  fails on I/O-limited nodes, notably on AWS.
- The same root cause needs the open-file limit raised to 2048, and
  `ALEVIN_FRY_HOME` set with `simpleaf set-paths` run. **Both nf-core (vendored)
  modules already do all of this in their own script blocks**. `scratch = true` on `SIMPLEAF_INDEX` is the only addition.
- Index size is modest: the full human 2024-A index is around 2 GB and
  piscem maps in under about 3 GB of RAM. Do not over-provision.
- Quant uses `--unfiltered-pl`, emitting a matrix with **no cell calling**
  applied. `--resolution` is required and defaults to `cr-like`.
  Under `--unfiltered-pl` the whitelist is the barcode-CORRECTION reference;
  only barcodes seen at least `--min-reads` times (simpleaf passes 10) become
  rows.
- Quant must pass **`--anndata-out`**, which writes
  **`af_quant/alevin/quants.h5ad`**. QCatch supports both mtx and h5ad
  input, but only the h5ad path carries the doublet and mitochondrial
  metadata we want. Without this flag we silently lose QCatch
  functionality. The nf-core (vendored) `simpleaf/quant` module hard-codes
  the flag, so `ext.args` must not repeat it.
- `--unfiltered-pl` takes an OPTIONAL whitelist path. simpleaf
  resolves and caches the whitelist for the chemistry under
  `ALEVIN_FRY_HOME`, which is a silent NETWORK dependency at run time.
  **Therefore, the whitelists are VENDORED** (now included in the pipeline's
  git repo) at `assets/whitelist/10x_V{2,3,4}_barcode_whitelist.txt.gz` (34.0 MB, from nf-core/scrnaseq 4.2.0) and staged as a channel input, so the path is
  always explicit. Never remove them and never fall back to the bare flag;
  that reintroduces a failure on offline nodes and locked-down AWS VPCs.
  **simpleaf accepts the whitelist GZIPPED.**
- `--expect-cells`, `--explicit-pl`, `--forced-cells`, `--knee` and
  `--unfiltered-pl` are MUTUALLY EXCLUSIVE in simpleaf. Passing more than
  one is a runtime error. Therefore the samplesheet's `expected_cells`
  column is metadata only and must never reach `simpleaf quant`; it exists
  for reporting and for comparison against the cell count QCatch calls.
- Cell calling is a separate explicit step (QCatch). Never add a knee or a
  fixed UMI cut-off to the quant step.
- The h5ad simpleaf writes is **Blosc-compressed**, and no QCatch image can
  read it. `REPACK_H5AD` rewrites it as gzip, so everything downstream
  consumes `*_raw_matrix.h5ad`, never `af_quant/alevin/quants.h5ad`. Do not
  remove that step: the compression is size-dependent, so a small fixture
  reads fine and every real sample does not.
- Matrix files carry a provenance suffix, following nf-core/scrnaseq:
  `*_raw_matrix.h5ad` from the quantifier, `*_filtered_matrix.h5ad` after
  cell calling. `REPACK_H5AD` applies the first one as it writes its
  output; the second is applied when QCatch's file is published.

### 4.3 Chemistry mapping — one source of truth

The samplesheet `chemistry` column drives BOTH the simpleaf chemistry and
the QCatch chemistry, through a single mapping in `conf/modules.config`.
Never let a user set them independently; they would drift apart.

| samplesheet | simpleaf    | QCatch      |
| ----------- | ----------- | ----------- |
| `10XV2`     | `10xv2`     | `10X_3p_v2` |
| `10XV3`     | `10xv3`     | `10X_3p_v3` |
| `10XV4`     | `10xv4-3p`  | `10X_3p_v4` |

VERIFIED with `docker run --rm -v -e ALEVIN_FRY_HOME=/tmp quay.io/biocontainers/simpleaf:0.30.0--hd612981_0 simpleaf chemistry lookup -n ".*"` and
`docker run --rm quay.io/biocontainers/qcatch:0.2.13--pyhdfd78af_0 qcatch --help` in the pinned containers .
The table above is correct as written and needs no further checking.

Note the v4 string: simpleaf registers `10xv4-3p`, NOT `10xv4`; there is
no `10xv4` entry at all. The registry holds exactly 17 chemistries:
`10xv2`, `10xv2-5p`, `10xv3`, `10xv3-5p`, `10xv4-3p`, the four
`10x-flexv*`/`10x-atac`/`10x-arc-atac` entries, and the Visium set. There
is no `smartseq` entry.
QCatch's `--chemistry` accepts `10X_3p_v2`, `10X_3p_v3`, `10X_3p_v4`,
`10X_5p_v3`, `10X_3p_LT` and `10X_HT`.

For custom assays with no QCatch mapping, omit `--chemistry` and set
`params.qcatch_n_partitions` instead, which maps to QCatch's
`--n_partitions` and supplies the partition count for its empty-drops
step.

### 4.4 Cell calling and QC

- **QCatch** performs cell calling with a two-step EmptyDrops-style
  ambient model. A knee threshold is not acceptable here: PBMCs have low
  RNA content (~1 pg/cell) and a global knee discards real low-count cells.
- ONE exception, and only one: `--cell_calling threshold` skips
  QCATCH and lets the Scanpy `min_genes` cut stand in. It exists so a
  fixture too small for an ambient model can still exercise the DAG —
  `params/test.yaml` uses it, `params/demo.yaml` and `params/default.yaml`
  do not. It is NOT a performance option and its output must never be
  described as cell-called. Only `qcatch` may produce
  `*_filtered_matrix.h5ad`. Never set it in a params file that produces
  results anyone will read.
- QCatch runs with `--skip_umap_tsne` (we do our own, seeded, in Scanpy),
  `--export_summary_table` and `--save_filtered_h5ad`.
- **Do NOT pass `--remove_doublets` or `--visualize_doublets` to QCatch.**
  Its Scrublet run exposes no seed. Doublets are detected in an explicit,
  seeded Scanpy step so the parameters are visible in the DAG. Running
  Scrublet twice would be wasteful and confusing.
- QCatch already writes `pct_counts_mt` into `adata.obs`. Recompute it in
  Scanpy anyway so the mitochondrial gene pattern is under our control,
  and note the duplication in a comment.
- Cell QC thresholds are MAD-based with absolute floors and ceilings from
  params. No hard-coded constants in the scripts.

### 4.5 Downstream

- **simpleaf writes `X` as the SUM of unspliced + spliced + ambiguous**,
  and stores U, S and A as separate named layers. `params.count_layer`
  (default `S+A`) must therefore be reconstructed from the named layers.
  Never read `X` directly and assume it matches `count_layer` — that is a
  silent correctness bug, not an error.
- Raw counts preserved in `adata.layers['counts']` before normalisation.
- HVG selection uses `flavor='seurat_v3'` on RAW counts, not log data.
- Doublet removal happens BEFORE HVG and PCA.
- Leiden uses `flavor='igraph'`, `n_iterations=2`.
- `params.seed` reaches scrublet, PCA, UMAP and Leiden. All four.
- `rank_genes_groups` is **marker detection** (cluster vs. rest), NOT
  differential expression. Never label it DE or DGE in code, comments,
  filenames, column headers or reports. A single-cell audience will treat
  that as an error.

## 5. Reporting

Three artifacts with three distinct jobs:

- **MultiQC** — QC only. Sources: seqkit stats, fastp, and custom content
  injected from the QCatch summary CSV and the Scanpy QC TSV. MultiQC has
  no module for alevin-fry, simpleaf or QCatch; verified against its
  supported-tools list. Custom content is the only route. Never put
  analysis results into MultiQC.
- **QCatch HTML** — per-sample cell-calling and mapping QC, shipped as-is,
  one per sample.
- **Quarto analysis report** — the biological result: UMAPs, cluster
  sizes, marker table.

Every table shown in a report is also written as a standalone TSV, and
every figure as PNG, into the results directory.

## 6. Reproducibility

- `versions.tsv`: columns `process`, `tool`, `version`, tab-delimited.
  `process` may repeat when one process uses several tools.
- Resolved parameters as `params_resolved.yaml`; run manifest (git commit,
  revision, command line, profile, container, Nextflow version) as
  `run_manifest.yaml`. YAML, not plain text.
- All software pinned to exact `version=build` in the Dockerfile.
- Nextflow `timeline`, `report`, `trace`, `dag` into `pipeline_info/`.

## 7. Privacy — this repo goes on public GitHub

No personal names, email addresses, institution names, local filesystem
paths, AWS account IDs, S3 bucket names, ARNs or credentials in ANY
tracked file.

One deliberate exception, added on request: `manifest.contributors` in
`nextflow.config` names the pipeline author. That is standard nf-core
practice and is the one place attribution belongs. Name only — no email,
GitHub handle, ORCID or affiliation unless asked for. Environment-specific AWS settings go in `conf/aws.local.config`, which is currently gitignored.

## 8. Style

- Comments explain WHY. Every non-obvious scientific choice gets a reason.
- No dead code, no commented-out blocks, no TODOs in committed code.
- Cross-references between documents quote the full heading text, not a
  bare section number.
- do NOT use emojis in any docs (also not Markdown files like README.md)

## 9. How to work

1. Read `.claude/specs/01-requirements.md`, `.claude/specs/02-design.md`,
   `.claude/specs/03-tasks.md`, in that order.
2. Work ONE task at a time. Mark `[~]` in progress, `[x]` done. Never
   batch tasks.
3. After every module you write, run
   `nextflow run . -profile test -stub-run`. Once `conf/test.config`
   exists, also run `nextflow run . -profile test,docker`, which uses
   small remote nf-core test data and finishes in minutes. They all need
   their respective `-params-file params/<profile>.yaml` files; see section 10. Paste the actual output both times.
4. Show real command output. Never assert something works without
   evidence.
5. If a spec is ambiguous, contradicts the documentation, or needs a
   version detail you have not verified in this session, STOP AND ASK.
   Never resolve by guessing, and never invent module names, container
   tags, CLI flags or URLs.
6. Never fix a failing check by weakening the check.
7. Never install software without asking first.

## 10. Commands

`nextflow lint -exclude` takes a bare directory path; the glob forms `.claude/*` and
`**/.claude/**` are silently NOT honoured.

```bash
nextflow lint . -exclude .claude
nextflow run . -profile test -stub-run -params-file params/test.yaml
nextflow run . -profile test,docker -params-file params/test.yaml
nextflow run . -profile demo,docker -params-file params/demo.yaml --run demo01
nf-core pipeline schema build
nf-core modules list remote | grep -i <tool>
nf-core modules install <tool>/<subtool>
```

## 11. Known-uncertain — verify against live docs, never assume

- KNOWN and FIXED, but easy to mistake for a pipeline bug if the fix is
  dropped: the alevin-fry 0.15.0 that nf-core simpleaf modules pin is a bioconda
  build of a `target-cpu=native` source tree, so it aborts with SIGILL on
  CPUs below the build worker's, taking both simpleaf processes with it.
  `conf/containers.config` overrides them onto a portable bioconda build.
- Whether heterogeneous channel shapes can be mixed into a single
  `output` entry for the curated results directory.
- The actual download size of the `pbmc_10k_v3` FASTQ tarball.
