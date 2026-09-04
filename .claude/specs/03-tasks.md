# Tasks — scrnaseq-lite

Rules:

- One task at a time. Mark `[~]` in progress, `[x]` done. Never batch.
- Tasks in the same wave are independent and may run in parallel git
  worktrees.
- Blocked by an entry under "Requirements section 4: Open questions
  requiring human decision"? STOP and ask, however they should all be
  resolved now.
- After every task: `nextflow run . -profile test -stub-run`, and once
  `conf/test.config` exists also `nextflow run . -profile test,docker`.
  Paste the real output.
- The effort column maps to the Claude Code CLI effort/thinking selector
  (Low / Medium / High / Extra High / Max).

## Wave 0 — Foundations (sequential)

- [ ] **T0.1 Environment probe.** Confirm Nextflow 26.04.6, docker,
      nf-core tools, and the active pyenv `scrnaseq-lite`. Report free
      space on the partition holding the Docker root directory and on the
      one holding this repo. Write `docs/environment.md` (gitignored).
      Do not touch `.claude/`. [Opus, Low]
- [ ] **T0.2 Repo scaffold.** Directory tree, skeleton. Placeholder files only. [Opus, Low]

## Wave 1 — Contracts (parallel)

- [ ] **T1.1 Schemas1.** `assets/schema_input.json`.
      Implements R1.1. Add one valid and one invalid example
      samplesheet. [Opus, High]
- [ ] **T1.2 Params.** `params/default.yaml` and `params/demo.yaml`.
      Implements R2.1 and N3. Pure YAML, one justification comment per
      scientific default, values from "Design section 7: Scientific
      defaults and their justification". [Opus, Medium]
- [ ] **T1.3 Core config.** `nextflow.config` and `conf/base.config`.
      Implements R2.2-R2.4, R12.3, R14.1. Verify that `outputDir` and
      `workflow.output.mode` are settable at top level as written.
      [Opus, High]
- [ ] **T1.4 Container.** `docker/Dockerfile`. Implements R14.3, R12.2,
      R4.5 and R5.7. micromamba base, multi-stage, non-root, every package
      pinned `version=build`. Set `ALEVIN_FRY_HOME` and run
      `simpleaf set-paths`. Raise the open-file limit. Decide and document
      whether the 10x barcode whitelist is pre-seeded into the image
      (R5.6) or fetched at run time. Build it; record the compressed size
      in `docs/container.md` (gitignored). Ask before pulling the base
      image. [Opus, High]
- [ ] **T1.5 Data profiles.** `conf/test.config` and `conf/demo.config`.
      Implements R15.3. `test` uses remote nf-core test data;
      browse `https://github.com/nf-core/test-datasets` on the `scrnaseq`
      branch and nf-core/scrnaseq 4.2.0 `conf/test.config` for the real
      paths and reference files. Do not construct URLs from memory.
      `demo` uses the local downsampled PBMC data from T5.1. Register both
      in the `profiles` block. [Opus, High]

## Wave 2 — Upstream (parallel)

- [ ] **T2.1 Vendored modules.** `nf-core modules install` for fastp,
      seqkit/stats, multiqc, simpleaf/index, simpleaf/quant and qcatch.
      Run `simpleaf chemistry lookup` (or the module's equivalent) inside
      the pinned container to confirm the registered chemistry strings
      before writing the mapping table; not sure if the v4 string is `10xv4-3p`, not
      `10xv4`. Make sure you update `conf/modules.config` and `.claude/specs/01-requirements.md`
      with this information.
      Report the actual installed versions and container tags. Determine
      whether these modules report versions via `versions.yml` or via
      topic channels, and record the answer in `.claude/specs/02-design.md` under
      "Design section 2: Module inventory"; all local modules must then
      match it. Never invent module names. [Opus, Medium]
- [ ] **T2.2 Module args.** `conf/modules.config`: `ext.args`,
      `ext.prefix`, `scratch = true` on `SIMPLEAF_INDEX` (R4.5), and the
      chemistry mapping from "CLAUDE.md section 4.3: Chemistry mapping —
      one source of truth". fastp args must not touch R1; comment the
      reason inline. Implements R3.1-R3.4, R5.1-R5.4, R6.4. [Opus, High]
- [ ] **T2.3 R1 assertion.** `modules/local/assert_r1_length/main.nf`.
      Implements R1.5. [Opus, Medium]
- [ ] **T2.4 Reference subworkflow.**
      `subworkflows/local/prepare_splici_reference.nf`. Implements R4.
      `bin/download_reference.sh` ALREADY EXISTS and fetches the 10x
      GRCh38 2024-A package. Do not rewrite it; read it so the subworkflow
      consumes the same paths (`<ref>/fasta/genome.fa`,
      `<ref>/genes/genes.gtf.gz`) and publishes an index referenced as
      `<dir>/index` (`.gitignore`'d). Confirm whether the installed simpleaf module sets
      `ALEVIN_FRY_HOME` and raises the open-file limit (R4.5, R5.7); if
      not, handle both in the container. [Opus, High]
- [ ] **T2.5 QCatch wiring.** The nf-core module is named `qcatch`.
      Inspect it with `nf-core modules info qcatch` before wiring, so the
      input/output channel shapes are taken from the module rather than
      assumed. Implements R6.1-R6.6. Do NOT pass `--remove_doublets` or
      `--visualize_doublets`. [Opus, High]

## Wave 3 — Scanpy (parallel)

- [ ] **T3.1 Cell QC.** `bin/scanpy_cell_qc.py` and its module.
      Implements R7.1-R7.4. Verify the QCatch h5ad layout empirically
      before writing the reader; do not assume the obs/var structure.
      [Opus, High]
- [ ] **T3.2 Doublets.** `bin/scanpy_detect_doublets.py` and its module.
      Implements R7.5-R7.7. Seeded. [Opus, Medium]
- [ ] **T3.3 Normalise and cluster.** `bin/scanpy_normalise_cluster.py`
      and its module. Implements R8 and R9. Seeded throughout.
      [Opus, High]
- [ ] **T3.4 Marker genes.** `bin/scanpy_marker_genes.py` and its module.
      Implements R10. No filename, header, comment or plot title may say
      DE or DGE. [Opus, Medium]
- [ ] **T3.5 Versions.** `bin/collect_versions.py` and its module ->
      `versions.tsv`. Implements R12.1. [Opus, Medium]
- [ ] **T3.6 Run manifest.** `bin/write_run_manifest.py` and its module ->
      `params_resolved.yaml` and `run_manifest.yaml`. Implements R2.5.
      [Opus, Medium]

## Wave 4 — Assembly (sequential)

- [ ] **T4.1 Subworkflows.** Explicitly named:
      `READ_QC`, `PREPARE_SPLICI_REFERENCE`, `QUANTIFY_AND_CALL_CELLS`,
      `CELL_QC_AND_DOUBLET_FILTERING`, `NORMALISE_CLUSTER_AND_MARKERS`,
      `REPORTING`. Annotate every channel operation with its in/out shape.
      [Opus, High]
- [ ] **T4.2 Entry workflow.** `workflows/scrnaseq_lite.nf` and `main.nf`
      with `publish:` and `output` blocks. Implements R13 and follows
      "Design section 5: Publishing" and "Design section 6: main.nf
      structure". Highest-consequence file in the repo. [Opus, Max]
- [ ] **T4.3 MultiQC.** `assets/multiqc_config.yaml` plus custom-content
      injection of the QCatch summary CSV and the Scanpy QC TSV.
      Implements R11.1. [Opus, High]
- [ ] **T4.4 Analysis report.** `report/analysis_report.qmd` and its
      module. Implements R11.3 and R11.4. [Opus, Medium]

## Wave 5 — Data, run, validate

- [ ] **T5.1 Schemas2.** `nextflow_schema.json` via `nf-core pipeline schema build`
      (pushed back here, as `nf-core pipeline schema build` requires a pipeline with
      `main.nf`).
      Implements R1.2-R1.4. [Opus, High]
- [ ] **T5.2 Test data.** `bin/download_and_downsample_testdata.sh`.
      Implements R15.4-R15.9. All dataset facts are already verified in
      "Requirements R15.4"; do not re-derive them, but DO read the actual
      FASTQ tarball sizes from the two 10x dataset pages and report them
      before downloading anything.
      Targets: both samples end at ~1,200 cells and ~24M read pairs.
      - `pbmc_1k_v3`: random read-pair subsample of 66,601,887 pairs,
        fraction ~0.36, fixed seed.
      - `pbmc_10k_v3`: download the filtered feature-barcode matrix, take
        ~1,200 of the 11,769 barcodes with a fixed seed, keep only read
        pairs whose R1 16 bp prefix is in that set, then random-subsample
        to ~24M pairs. Combine both filters into ONE streaming pass from
        the 10x CDN so the ~9.6x-larger source never lands on disk.
      - Use a short Python helper reading both mates in lockstep; do not
        use `seqkit grep` for the barcode filter.
      - `--simple` falls back to random subsampling for both, warning that
        the 10k sample will be shallow.
      Assert equal read counts and identical read-name order in the R1/R2
      outputs before writing. Record the CC BY 4.0 attribution for both
      datasets. Idempotent and re-runnable. [Opus, High]
- [ ] **T5.3 Green run.** Full `-stub-run`, then `-profile test,docker`,
      then the real `-profile demo,docker` two-sample run. Fix to green.
      [Opus, High]
- [ ] **T5.4 Validation.** Run `/validate-pipeline`. Produce
      `docs/validation-report.md` (gitignored). [Opus, Max]
- [ ] **T5.5 Documentation.** `README.md`, `docs/usage.md`,
      `docs/output.md`. Must state the `--run` / `-resume` requirement
      (R2.4), the difference between the `test` and `demo` profiles
      (R15.3), the two subsampling strategies and why they differ (R15.5),
      and the dataset licence and attribution (R15.9). Add one line to
      `README.md` noting that the pipeline was built spec-driven development
      and that the requirements, design and task documents are in `.claude/specs/`.
      No private or environment-specific information (N4). [Opus, Medium]
- [ ] **T5.6 Changelog.** Finalise `CHANGELOG.md` 0.1.0 with today's date.
      A single "Initial implementation" entry only. [Opus, Low]
