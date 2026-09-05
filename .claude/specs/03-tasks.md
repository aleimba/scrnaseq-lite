# Tasks — scrnaseq-lite

Single source of truth for task status and for every finding that changes
what a later task must do. Marks: `[ ]` not started, `[~]` in progress,
`[x]` done.

Rules:

- One task at a time. Mark `[~]` in progress, `[x]` done. Never batch.
- Tasks in the same wave are independent and may run in parallel git
  worktrees.
- Blocked by an entry under "Requirements section 4: Open questions
  requiring human decision"? STOP and ask, however they should all be
  resolved now.
- After every task: `nextflow run . -profile test -stub-run`, and once
  `conf/test.config` exists also `nextflow run . -profile test,docker`.
  Paste the real output. NOTE: neither is possible before T4.2, which
  creates `main.nf` — until then `nextflow run .` has no entry point.
  `nextflow lint . -exclude .claude` is the check that does work.
- The effort column maps to the Claude Code CLI effort/thinking selector
  (Low / Medium / High / Extra High / Max).

**Recording findings.** A finding that changes a later task goes in the
task entry below, marked `CARRIED FROM <task>` or `SETTLED AT <task>`, and
in the requirement it invalidates in `.claude/specs/01-requirements.md`.
This file is tracked and public, so it must contain NO machine-specific
detail — no local paths, disk figures, or host state (N4). Those belong in
the gitignored `docs/environment.md` and `docs/container.md`.

## Wave 0 — Foundations (sequential)

- [x] **T0.1 Environment probe.** Confirm Nextflow 26.04.6, docker,
      nf-core tools, and the active pyenv `scrnaseq-lite`. Report free
      space on the partition holding the Docker root directory and on the
      one holding this repo. Write `docs/environment.md` (gitignored).
      Do not touch `.claude/`. [Opus, Low]
- [x] **T0.2 Repo scaffold.** Directory tree, skeleton. Placeholder files
      only. [Opus, Low]
      Delivered: dir tree, `.gitkeep` placeholders, `.gitignore`,
      `CHANGELOG.md` skeleton (date filled in at T5.6).
      FINDING: the lint command for this repo is
      `nextflow lint . -exclude .claude`. `-exclude` takes a bare
      directory path; the glob forms `.claude/*` and `**/.claude/**` are
      silently NOT honoured.

## Wave 1 — Contracts (parallel)

- [x] **T1.1 Schemas1.** `assets/schema_input.json`.
      Implements R1.1. Add one valid and one invalid example
      samplesheet. [Opus, High]
      Delivered `assets/schema_input.json` plus
      `assets/samplesheet_example_{valid,invalid}.csv`. Verified with
      nf-schema 2.8.0 via a throwaway `samplesheetToList` probe (removed).
      CARRIED FROM T1.1:
      1. A BLANK optional integer column validates fine against
         `"type": "integer"` and arrives in meta as `[]` — not `null` and
         not absent. Downstream code must treat `meta.expected_cells == []`
         as "not provided" (in Groovy `[]` is falsy).
      2. `uniqueEntries` works at the TOP level of the schema with a
         sibling top-level `errorMessage`. Wrapping it in `allOf`, as the
         nf-schema docs example does, also works but adds a spurious
         "Value does not match against the schemas at indexes [0]" line to
         the error output.
- [x] **T1.2 Params.** `params/default.yaml` and `params/demo.yaml`.
      Implements R2.1 and N3. Pure YAML, one justification comment per
      scientific default, values from "Design section 7: Scientific
      defaults and their justification". [Opus, Medium]
      36 keys each, identical key sets, differing in exactly 3 values
      (`input`, `simpleaf_index`, `multiqc_title`). T1.5 added
      `params/test.yaml` as a third; parity is now three-way. Re-run this
      whenever a param is added in a later wave:
      `python -c "import yaml; ks=[set(yaml.safe_load(open(f))) for f in ('params/default.yaml','params/demo.yaml','params/test.yaml')]; assert ks[0]==ks[1]==ks[2], ks[0]^ks[1]|ks[0]^ks[2]"`
      CARRIED FROM T1.2:
      1. Params from a `-params-file` keep their YAML types (null stays
         null, 1.0 is a Double, 42 an Integer). Params from the CLI arrive
         as **String**: `--leiden_resolution 0.4` binds as `"0.4"`. Safe
         when the value is interpolated into a shell command and cast by
         argparse, but any Groovy arithmetic on a param must cast first.
      2. `--run demo01` overrides correctly and no `run` key leaks from
         either YAML.
      Deviations from the original reference draft, all resolving a
      conflict with CLAUDE.md in favour of CLAUDE.md: dropped `run`
      (`nextflow.config` plus CLI only), dropped `fastp_args` (belongs in
      `conf/modules.config` `ext.args`), added `mito_gene_prefix` /
      `ribo_gene_pattern` / `hb_gene_pattern` so the T3.1 QC script holds
      no hard-coded gene patterns.
- [x] **T1.3 Core config.** `nextflow.config` and `conf/base.config`.
      Implements R2.2-R2.4, R12.3, R14.1. Verify that `outputDir` and
      `workflow.output.mode` are settable at top level as written.
      [Opus, High]
      Also delivered `conf/aws.config`. Verified by running: `nextflow
      config` for the `docker` and `awsbatch` profiles, plus a throwaway
      root `main.nf` stub probe (removed) run bare, with `--run probe01`,
      and with `--run probe01 -resume`. Lint clean on all three files.
      CARRIED FROM T1.3 — Nextflow 26 config-parser behaviour. These are
      not style preferences; each one is a hard parse error:
      1. `def` is rejected outright at config top level: "Variable
         declarations cannot be mixed with config statements". Even a
         config containing nothing but a `def` fails.
      2. `"${outputDir}"` inside a `timeline`/`report`/`trace`/`dag` block
         is a hard parse error (`` `outputDir` is not defined ``), not a
         lint false positive. `outputDir` is write-only: settable, never
         readable from a later statement.
      3. Consequence: `params` is the only scope readable by both
         `outputDir` and the trace blocks, so the run directory is
         computed once into the derived param `run_output_dir`, declared
         in `nextflow.config` only. It must never enter a params file.
         T5.1 should mark it `hidden` in `nextflow_schema.json`.
      4. Inside a `params` block a sibling key must be read as
         `params.run`, not bare `run`.
      5. `file()` is not available to the config parser. The conditional
         `conf/aws.local.config` include uses
         `new File("${projectDir}/...").exists()`. `projectDir` matters:
         a bare relative path resolves against the launch directory, so
         the file would be silently skipped.
      6. A missing `includeConfig` target is a hard error even inside an
         INACTIVE profile. This is why T1.5 added the `test`/`demo`
         profile blocks together with the files they include, and why T2.2
         must add the `conf/modules.config` include together with the file.
      7. `overwrite = true` is required on all four trace blocks. Without
         it `--run <name> -resume` — the supported resume workflow (R2.4)
         — ends with "Report file already exists" warnings and leaves
         stale provenance.
      8. A process with a `stub:` but no `script:` fails to compile
         ("Invalid process definition"). Every module in Waves 2-3 needs
         both.
      9. Lint warnings to obey in T4.2: use lowercase `channel`, not
         `Channel`; prefix unused output-closure parameters with `_`.
      CLOSED at T1.3 (was a CLAUDE.md section 11 open question): both
      `outputDir = ...` and `workflow.output.mode = 'copy'` ARE settable
      at top level. Published files are real copies, not symlinks;
      `pipeline_info/` lands inside the run directory (R12.3);
      `--run probe01 -resume` gave `cached=1` in the same directory (R2.4).
      DECISION, with knock-on effects: per-module biocontainers, so NO
      global `process.container`. R14.3 is therefore only partly met and
      our image covers local modules only. See R14.3 as amended.
      `conf/aws.config` deviates from "Design section 11: AWS readiness"
      by using the native `aws.*`/`workDir` scopes rather than params, so
      deployment settings stay out of the complete-parameter-set files.
      Verified against live Nextflow docs: `aws.batch.cliPath` is
      documented as a legacy requirement, superseded by Wave plus Fusion.
      `nf-schema` pinned at 2.8.0 — the version T1.1 and T1.2 were
      actually verified against, not the template's 2.5.1.
- [x] **T1.4 Container.** `docker/Dockerfile`. Implements R14.3 and
      R12.2. micromamba base, multi-stage, non-root, every package pinned
      `version=build`. Build it; record the compressed size in
      `docs/container.md` (gitignored). Ask before pulling the base
      image. [Opus, High]
      SCOPE (set at T1.3, corrected at T1.5): per-module biocontainers
      were chosen over one global image, so `SIMPLEAF_INDEX`,
      `SIMPLEAF_QUANT` and `QCATCH` never run in this image. It covers the
      LOCAL modules only: the four Scanpy steps, `ASSERT_R1_LENGTH`,
      Quarto, and the version and manifest scripts. Nothing
      simpleaf-related belongs here — no `ALEVIN_FRY_HOME`, no `simpleaf
      set-paths`, no `ulimit`, no whitelist. The nf-core simpleaf modules
      do all of that in their own script blocks; see R4.5, R5.6 and R5.7.
      DELIVERED: `scrnaseq-lite:0.1.0`, base
      `mambaorg/micromamba:2.9.0-debian13` pinned by digest, 243 packages,
      2.66 GB uncompressed and **628 MB gzip-compressed** against R14.3's
      2 GB budget. No registry prefix: built locally, pushed to Amazon ECR
      separately, with the account-bearing URI confined to the gitignored
      `conf/aws.local.config` (N4).
      Verified by running, against the final image with no `-e`
      overrides: non-root smoke test as an arbitrary host uid; the seeded
      scientific paths (S+A layer reconstruction, scrublet, `seurat_v3`
      HVG on RAW counts, PCA, UMAP, Leiden `igraph` `n_iterations=2`,
      `rank_genes_groups`, h5ad round-trip, PNG) BIT-IDENTICAL across two
      separate processes; and a real `quarto render` of a Python `.qmd`
      producing an embedded PNG. That determinism run is what proves
      scanpy 1.12.4 actually works on python 3.14.7 / pandas 3.0.5 /
      numpy 2.5.2 — the solver only proved they were installable together.
      Also delivered here, out of band at user request: the 10x barcode
      whitelists are VENDORED at
      `assets/whitelist/10x_V{2,3,4}_barcode_whitelist.txt.gz`, closing
      the R5.6 run-time network dependency. From nf-core/scrnaseq tag
      `4.2.0`; all three byte counts matched the GitHub contents API
      exactly, every line a bare 16 bp `[ACGT]` barcode, zero malformed.
      737,280 / 6,794,880 / 7,372,800 barcodes, 34.0 MB total. V4 has MORE
      barcodes than V3 yet compresses smaller. V1 omitted: the R1.1 schema
      enum cannot select it. Provenance and sha256 in
      `assets/whitelist/README.md`. Corrected while vendoring: R5.6 said
      the cost was "roughly 20 MB"; 2.2 + 18.4 + 13.4 is 34.0 MB.
      CARRIED FROM T1.4 — container construction, all found by running:
      1. `micromamba install -p <prefix>` REFUSES a prefix that does not
         exist ("Environment must first be created with micromamba
         create"). Use `create`. It is still one layer.
      2. The `mambaorg/micromamba` base ends on `USER mambauser` (uid
         57439). A second build stage inherits that user, so `useradd` and
         any write to `/opt` fail. Stage 2 needs `USER root` first and
         drops to `appuser` at the end. uid 1000 is free in this base.
      3. SKIPPING CONDA ACTIVATION BREAKS QUARTO, in three stages. Putting
         the env's `bin` on `PATH` rather than activating is fine for
         Python and not for quarto: conda-forge ships
         `etc/conda/activate.d/quarto.sh` exporting eight `QUARTO_*`
         variables, and without them quarto hunts for its bundled tools
         under `$QUARTO_BIN_PATH/tools/x86_64/`, which no conda layout
         has. The failures are staggered, so a partial fix looks like
         progress: share path missing -> `quarto --version` dies on
         `cat: .../share/version`; `QUARTO_DENO` missing ->
         `tools/x86_64/deno: No such file or directory`;
         `QUARTO_DART_SASS` missing -> dies in `compileSass` AFTER the
         Python kernel has run, so the log reads like success. Set the
         whole set. One deliberate deviation: `activate.d/quarto.sh` ships
         a BROKEN `QUARTO_DENO_DOM` pointing into the feedstock build
         directory, absent from the image; `deno_dom.sh` has the real path.
      4. `sc.pp.scrublet` requires `scikit-image` and says so only at RUN
         time, not on import: "`threshold` is None and thus scrublet
         requires `scikit-image`". R7.5 mandates that exact call, so the
         package was genuinely missing. Re-solving with it moved no other
         pin (closure 220 -> 243).
      5. fontconfig does not create its own cache directory; it probes a
         fixed list and gives up. Pre-create it in the image.
      6. `docker.runOptions = '-u $(id -u):$(id -g)'` means an arbitrary
         host uid with no `/etc/passwd` entry and no writable home. `HOME`,
         `MPLCONFIGDIR`, `NUMBA_CACHE_DIR` and the `XDG_*` variables must
         point at world-writable paths or scanpy fails to IMPORT.
      7. `ps` reports `ps from procps-ng UNKNOWN` — it runs, which is all
         Nextflow needs for trace metrics, but it carries no version
         string. Relevant to T3.5: do not try to version-report it.
      8. REPRODUCIBILITY GAP, open and deliberate: the Dockerfile pins the
         21 packages it names, not the full 243-package closure. `numba`,
         `llvmlite` and `libopenblas` are transitive and unpinned, and all
         three affect numerical output, so a rebuild months later is not
         guaranteed bit-identical. A given built image IS reproducible by
         digest. Fix if wanted: install from a generated lockfile.
         Recorded in R12.2 so the choice stays visible.
      R12.4 MEASURED, not assumed: with `OMP_NUM_THREADS=1`, PCA changed
      and UMAP followed, while scrublet, HVG, Leiden and markers were
      unchanged. Document "PCA and UMAP may differ across BLAS thread
      counts"; do NOT promise cluster stability, which held only on
      well-separated synthetic data.
- [x] **T1.5 Data profiles.** `conf/test.config` and `conf/demo.config`,
      plus `params/test.yaml`, `assets/samplesheet_test.csv` and
      `assets/samplesheet_demo.csv`. Implements R15.3.
      `conf/*.config` carry ENGINE settings only (`resourceLimits`); all
      `params.*` live in the profile's own complete params file, because
      `-params-file` outranks config-defined params (R2.1). Register both
      profiles in the `profiles` block — note that a missing
      `includeConfig` target is a hard parse error even in an INACTIVE
      profile, so the block and the file must land together.
      Real paths come from the `scrnaseq` branch of
      `https://github.com/nf-core/test-datasets` and from nf-core/scrnaseq
      4.2.0 `conf/test.config`. Do not construct URLs from memory.
      `demo` uses the local downsampled PBMC data from T5.2. [Opus, High]
      Also added the `test`/`demo` profile blocks and
      `manifest.contributors` to `nextflow.config`. Lint clean on all five
      config files.
      CARRIED FROM T1.5, verified by running:
      1. **`-params-file` OUTRANKS params set in a config file**, and an
         explicit `null` in the params file still wins. Probed with a
         throwaway `-c` config setting `params.input = 'CONFIG_WINS'` plus
         `-params-file params/default.yaml`; the result was `null`. This
         is why `conf/test.config` and `conf/demo.config` carry ENGINE
         settings only: a profile config that set `params.input` would be
         silently nulled by any params file that did not.
      2. nf-schema's `exists: true` and `format: file-path` DO accept
         https URLs. `samplesheetToList` on the test samplesheet returned
         full `https://...` URIs intact, so remote staging will work.
      3. Both profiles resolve; `resourceLimits` is 2 CPU / 6 GB / 2 h for
         `test` and 2 CPU / 8 GB / 1 h for `demo` (R14.4). Neither config
         leaks a single `params.*`.
      Test-data facts, read from the live sources and never constructed:
      GRCm38 chr19 FASTA plus gencode vM19 GTF (83 MB together), the two
      `Sample_{X,Y}_S1_L001` FASTQ pairs (~6 MB each), chemistry `10XV2`.
      **R1 is 26 bp and R2 is 90 bp**, read from the FASTQ stream — so
      `r2_read_length` is 90 here, not the PBMC 91 — and the reference is
      MOUSE, so `mito_gene_prefix` is `mt-`.
      `assets/samplesheet_test.csv` is ours because nf-core's has no
      `chemistry` column and lists `Sample_Y` twice for two lanes, which
      our `uniqueEntries` rule rejects by design (R1.4). Only L001 is
      used. `assets/samplesheet_demo.csv` cannot validate until T5.2
      creates the FASTQs it names; that is expected, not a bug.
      KNOWN CEILING on `-profile test,docker`: see R15.3 as amended, and
      T2.5, which settles it.

## Wave 2 — Upstream (parallel)

- [x] **T2.1 Vendored modules.** `nf-core modules install` for fastp,
      seqkit/stats, multiqc, simpleaf/index, simpleaf/quant and qcatch.
      Run `simpleaf chemistry lookup` (or the module's equivalent) inside
      the pinned container to confirm the registered chemistry strings
      before writing the mapping table; not sure if the v4 string is
      `10xv4-3p`, not `10xv4`. Make sure you update `conf/modules.config`
      and `.claude/specs/01-requirements.md` with this information.
      Report the actual installed versions and container tags. Determine
      whether these modules report versions via `versions.yml` or via
      topic channels, and record the answer in `.claude/specs/02-design.md`
      under "Design section 2: Module inventory"; all local modules must
      then match it. Never invent module names. [Opus, Medium]
      CARRIED FROM T1.5 — confirm these against the INSTALLED versions,
      all read from nf-core/scrnaseq 4.2.0:
      (a) `simpleaf/index` and `simpleaf/quant` already do
          `export ALEVIN_FRY_HOME=.`, `ulimit -n 2048` (index only) and
          `simpleaf set-paths` in-script, so T2.2 needs only
          `scratch = true` (R4.5, R5.7).
      (b) `simpleaf/quant` hard-codes `--anndata-out`; `qcatch` hard-codes
          `--save_filtered_h5ad` and `--export_summary_table`. Do not
          repeat them in `ext.args` (R5.3a).
      (c) `simpleaf/quant` takes the whitelist as a staged channel input
          and emits `--unfiltered-pl <file>`. SETTLED AT T1.4: the files
          are vendored at `assets/whitelist/`; R5.6 is closed.
      (d) THE TWO VERSION CONVENTIONS ARE MIXED: simpleaf modules emit
          `versions.yml`, `qcatch` emits to the `versions` topic (R12.1a).
          CARRIED FROM T1.4 — (d) IS PROBABLY OBSOLETE. It was read from
          nf-core/scrnaseq 4.2.0, which pins an OLD module snapshot:
          simpleaf 0.19.5, container
          `quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0`, output
          `path "versions.yml"`. Current nf-core/modules ships simpleaf
          **0.25.0**, container
          `community.wave.seqera.io/library/simpleaf:0.25.0--b9f96d8b71a01864`,
          emitting `topic: versions` for alevin-fry, piscem and simpleaf
          alike — the same convention as `qcatch`. `nf-core modules
          install` fetches the latter, so the conventions are likely NOT
          mixed. Confirm against what actually installs, then correct
          R12.1a and T3.5 rather than assuming either way.
      (e) CARRIED FROM T1.4: nf-core/scrnaseq 4.2.0 `assets/protocols.json`
          independently confirms the whole mapping table in "CLAUDE.md
          section 4.3", `10xv4-3p` included — `10XV2`->`10xv2`/`10X_3p_v2`,
          `10XV3`->`10xv3`/`10X_3p_v3`, `10XV4`->`10xv4-3p`/`10X_3p_v4`.
          That is a second independent source, not a substitute for
          running `simpleaf chemistry lookup` in the pinned container.
      DELIVERED: all six modules installed from `nf-core/modules` `master`
      with nf-core tools 4.1.0. Versions and containers, read from the
      installed `main.nf` and `environment.yml` and confirmed by running
      each container:

      | Module          | tool versions                              | container (docker)                                              |
      | --------------- | ------------------------------------------ | --------------------------------------------------------------- |
      | `fastp`         | fastp 1.3.6                                | `community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde` |
      | `seqkit/stats`  | seqkit 2.13.0                              | `community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751` |
      | `multiqc`       | multiqc 1.35                               | `community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc` |
      | `simpleaf/index`| simpleaf 0.25.0, alevin-fry 0.15.0, piscem 0.20.0 | `community.wave.seqera.io/library/simpleaf:0.25.0--b9f96d8b71a01864` |
      | `simpleaf/quant`| simpleaf 0.25.0, alevin-fry 0.15.0, piscem 0.20.0 | same as `simpleaf/index`                                    |
      | `qcatch`        | qcatch 0.2.12 (python 3.14.6)              | `community.wave.seqera.io/library/pip_qcatch:03b88593a5cca75b`  |

      The pinned `git_sha` of each is in `modules.json`. CLOSED (was a
      CLAUDE.md section 11 open question): every one of the six is a
      Seqera Wave community image, none is a `quay.io/biocontainers` image.
      CARRIED FROM T2.1 — findings, every one produced by running:
      1. **THE VERSION CONVENTIONS ARE NOT MIXED.** All six modules report
         versions with `eval(...)` into `topic: versions`. Not one
         `versions.yml` exists in `modules/nf-core/`. Note (d) above is
         confirmed obsolete and R12.1a is corrected accordingly. Local
         modules use the topic convention, and T3.5 reads the topic only —
         with ONE exception, next.
      2. `MULTIQC` deliberately does NOT publish to the topic. It carries
         the comment "MultiQC should not push its versions to the
         `versions` topic. Its input depends on the versions topic to be
         resolved thus outputting to the topic will let the pipeline hang
         forever", and emits a plain `emit: versions` channel instead.
         T3.5 and T4.2 must therefore merge the topic with MULTIQC's own
         channel, and MULTIQC must run AFTER `COLLECT_VERSIONS` consumes
         the topic. Wiring MULTIQC into the topic deadlocks the run.
      3. `simpleaf chemistry lookup -n ".*"` in the pinned container lists
         exactly 17 registered chemistries. `10xv2`, `10xv3` and
         `10xv4-3p` are registered; **`10xv4` is NOT**. The mapping table
         in "CLAUDE.md section 4.3" is confirmed against the tool itself,
         which is the third independent source after nf-core/scrnaseq
         `assets/protocols.json` (note (e)).
      4. SCIENTIFIC CORRECTION, from the same lookup: the registered
         geometries are `10xv2` -> `1{b[16]u[10]x:}2{r:}` and `10xv3` /
         `10xv4-3p` -> `1{b[16]u[12]x:}2{r:}`. **The v2 UMI is 10 bp, not
         12.** R1 lengths are unchanged (v2 26 bp = 16+10, v3/v4 28 bp =
         16+12), so R1.5 and the T2.3 assertion still hold, but the
         "16 bp barcode + 12 bp UMI" gloss in "CLAUDE.md section 4.1" is
         wrong for v2 and is corrected there and in R1.5.
      5. `qcatch --help` confirms `--chemistry` accepts `10X_3p_v2`,
         `10X_3p_v3`, `10X_3p_v4` (also `10X_5p_v3`, `10X_3p_LT`,
         `10X_HT`), and confirms `--n_partitions` and `--skip_umap_tsne`.
         `params.qcatch_n_partitions` maps to `--n_partitions` (R6.4).
      6. **BLOCKER on this host, not in the pipeline.** The bioconda
         `alevin-fry 0.15.0 hd612981_0` build in the pinned simpleaf
         container dies with SIGILL ("Illegal instruction (core dumped)",
         exit 132) on `alevin-fry --version`. The host CPU is Zen 2 —
         `avx`, `avx2` and `bmi2` but no AVX-512 — so the build uses an
         instruction set this CPU lacks. `simpleaf` and `piscem` from the
         same image run fine, and `alevin-fry 0.11.2` from
         `quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0` runs fine, so
         it is that one build. CONSEQUENCE: `simpleaf set-paths` itself
         exits 1 ("Error running alevin-fry ... signal: 4 (SIGILL)"), and
         both simpleaf modules run it in their script blocks under
         `bash -ue`, so `SIMPLEAF_INDEX` and `SIMPLEAF_QUANT` CANNOT RUN on
         a CPU without that instruction set. Verified end to end with a
         throwaway Nextflow process (removed): `[ERROR] SET_PATHS_PROBE
         exit: 1`. Affects T2.4, T2.5 and above all T5.3, whose real
         `-profile test,docker` and `demo,docker` runs need a host whose
         CPU the build supports. `-stub-run` is unaffected. Do NOT edit the
         vendored modules to work around it. **FIXED at T2.1; see 6a.**
      6a. **FINDING 6 IS FIXED, in `conf/containers.config`.**
         Root cause, read from the sources rather than inferred: alevin-fry
         at tag `v0.15.0` ships a `.cargo/config.toml` containing
         `rustflags = ["-C", "target-cpu=native"]`, and the bioconda recipe
         builds it with a plain `cargo install --locked --path .` setting
         no RUSTFLAGS of its own, so the published binary is tuned to the
         build worker's CPU. Upstream fixed exactly this in alevin-fry
         commit `0e13d52` (2026-07-10), "make .cargo/config.toml portable
         by default (drop target-cpu=native)", whose stated reason is that
         it "can fault with SIGILL on machines with fewer CPU features".
         `v0.18.1` and master build the portable baseline
         `target-cpu=x86-64-v3 -C target-feature=+avx2` instead.
         nf-core/scrnaseq 4.2.0 never hits this because it pins
         `quay.io/biocontainers/simpleaf:0.19.5--ha6fb395_0`, carrying
         alevin-fry 0.11.2 — OLDER than the offending config, not newer
         than the fix. Copying their pin would cost four simpleaf minor
         versions; we go forward instead.
         THE FIX: bioconda's simpleaf 0.30.0 recipe requires
         `alevin-fry >=0.18.1` and `piscem >=0.23.0`, and
         `quay.io/biocontainers/simpleaf:0.30.0--hd612981_0` resolves to
         simpleaf 0.30.0 / alevin-fry 0.18.1 / piscem 0.23.0.
         `conf/containers.config` overrides `container` for
         `SIMPLEAF_INDEX|SIMPLEAF_QUANT` to that image and `nextflow.config`
         includes it. The vendored modules stay untouched. It lives in its
         own file rather than in `conf/modules.config`, which carries
         `ext.args` and `ext.prefix` only, so the deviation stays visible.
         VERIFIED by running, on the host that reproduces the SIGILL:
         `simpleaf set-paths` exits 0 in the new image; a throwaway process
         NAMED `SIMPLEAF_QUANT` (removed), running the modules' own
         `export ALEVIN_FRY_HOME=. ; simpleaf set-paths` sequence under the
         repo's `process.shell`, gave `[SUCCESS] completed=1 failed=0` and
         `TOPIC: SIMPLEAF_QUANT alevin-fry '0.18.1'`;
         `nextflow config -profile docker` resolves the override; lint
         clean at 18 files.
         NO CLI DRIFT between 0.25.0 and 0.30.0 for anything the modules
         use: `--fasta`, `--gtf`, `--threads`, `-o` on index and
         `--chemistry`, `--index`, `--reads1`, `--reads2`, `--t2g-map`,
         `--resolution`, `--output`, `--threads`, `--anndata-out`,
         `--unfiltered-pl`, `--map-dir`, `--min-reads` on quant all still
         exist, and the chemistry registry is identical — same 17 entries,
         `10xv4-3p`, no `10xv4`.
         ONE NEW OPTION in 0.30.0, recorded in R5.3: `--small-thresh <N>`
         resolves cells below N with `cr-like` (winner-take-all) semantics
         REGARDLESS of `--resolution`. Inert at our `cr-like` default; it
         would silently apply under `cr-like-em`. Left unset, alevin-fry's
         own default applies.
         DELETE `conf/containers.config` and its `includeConfig` line once
         `nf-core modules update` brings a simpleaf module pinning
         alevin-fry 0.18.1 or newer. Trap while editing that file: a `*/`
         inside a block comment — as in a path glob like
         `simpleaf/*/environment.yml` — closes the comment and is a hard
         parse error.
         STILL TRUE, and why the override earns its cost: any bioconda
         package built from a `target-cpu=native` source tree can fail this
         way, silently and only on some hosts, and `-stub-run` never
         exercises it.
      7. A FAILING `eval:` version command does NOT fail the task. The
         probe in finding 6 emitted `TOPIC: EVAL_PROBE alevin-fry ''` and
         the run reported `[SUCCESS] completed=1 failed=0`. So a broken
         tool silently becomes an EMPTY version string. T3.5 must treat an
         empty version as an error or an explicit `unknown`, never write a
         blank cell into `versions.tsv`.
      8. `qcatch --version` prints `qcatch version 0.2.12`, and the
         module's `sed -e 's/qcatch //g'` leaves **`version 0.2.12`** —
         verified with `cat -A`. Every other module's eval is clean
         (`1.3.6`, `2.13.0`, `0.20.0`, `0.25.0`). T3.5 must strip a
         leading `version ` token; we cannot fix the vendored module.
      9. `nf-core modules install` also wrote eight `conf/containers_*.config`
         files (docker / singularity-https / singularity-oras / conda-lock
         x amd64 / arm64), each listing only `FASTP` and `MULTIQC`. Nothing
         includes them and they duplicate the container already declared in
         each module's `main.nf`, which breaks "a value appears in exactly
         one file". They were DELETED. They will reappear on the next
         `nf-core modules install`; delete them again.
      10. Module input shapes, taken from the installed files, for T2.4,
          T2.5 and T4.1 — do not assume the older shapes:
          - `FASTP`: `tuple val(meta), path(reads), path(adapter_fasta)`
            plus three `val` flags (`discard_trimmed_pass`,
            `save_trimmed_fail`, `save_merged`).
          - `SEQKIT_STATS`: `tuple val(meta), path(reads)`; `ext.args`
            defaults to `--all`.
          - `MULTIQC`: `tuple val(meta), path(multiqc_files),
            path(multiqc_config), path(multiqc_logo), path(replace_names),
            path(sample_names)` — one tuple, meta included.
          - `SIMPLEAF_INDEX`: four tuples (genome fasta+gtf, transcript
            fasta, probe csv, feature csv); supplying fasta+gtf yields
            `--fasta ... --gtf ...`, which is the splici path (R4.2).
            Outputs `${prefix}/index`, optional `${prefix}/ref` and
            optional `${prefix}/ref/{t2g,t2g_3col}.tsv`.
          - `SIMPLEAF_QUANT`: `tuple val(meta), val(chemistry), path(reads)`;
            `tuple val(meta2), path(index), path(txp2gene)`;
            `tuple val(meta3), val(cell_filter), val(number_cb), path(cb_list)`;
            `val resolution`; `tuple val(meta4), path(map_dir)`. Reads must
            arrive as a flat R1,R2,R1,R2 list — the module does
            `reads.collate(2).transpose()`. It ADDS a `filtered` key to
            meta. `cell_filter = 'unfiltered-pl'` with the vendored
            whitelist as `cb_list` produces `--unfiltered-pl <file>`;
            note (c) confirmed.
          - `QCATCH`: `tuple val(meta), val(chemistry), path(quant_dir)`.
      11. Note (a) confirmed with one correction: `simpleaf/index` does
          `export ALEVIN_FRY_HOME=.`, `ulimit -n 2048` and
          `simpleaf set-paths`; `simpleaf/quant` does the first and the
          third but **no `ulimit`**. Note (b) confirmed verbatim:
          `--anndata-out` is hard-coded in quant, `--save_filtered_h5ad`
          and `--export_summary_table` in qcatch.
      12. `QCATCH` renames its outputs to `${prefix}_qcatch_report.html`,
          `${prefix}_filtered_quants.h5ad` and
          `${prefix}_metrics_summary.csv`. That is `_filtered_quants`, NOT
          the `*_filtered_matrix.h5ad` R6.6 and "CLAUDE.md section 4.2"
          require, and the quantifier's matrix is
          `af_quant/alevin/quants.h5ad` inside a directory (path corrected
          at T2.2), not `*_raw_matrix.h5ad`. The vendored modules
          cannot be edited, so the provenance suffix must be applied where
          the file is PUBLISHED (T4.2 `output` block) or by the local
          Scanpy module that consumes it. Decide once, at T4.2, and use the
          same rule for both matrices.
      13. QCatch has a run-time NETWORK dependency of its own: without
          `--gene_id2name_file` it "will attempt to retrieve the mapping
          from a remote registry", and if that fails the mitochondrial
          plots are dropped from the HTML — quietly, not as an error. The
          same class of problem the vendored whitelists solved at T1.4.
          T2.5 must decide whether to feed it the GTF-derived mapping.
- [x] **T2.2 Module args.** `conf/modules.config`: `ext.args`,
      `ext.prefix`, `scratch = true` on `SIMPLEAF_INDEX` (R4.5 — and ONLY
      that; the module already sets `ALEVIN_FRY_HOME` and the open-file
      limit itself, so no `beforeScript` is needed. This supersedes an
      earlier, incorrect amendment of this task), and the chemistry
      mapping from "CLAUDE.md section 4.3: Chemistry mapping — one source
      of truth". Do not repeat flags the modules hard-code; see T2.1 note
      (b). Also add `includeConfig 'conf/modules.config'` to
      `nextflow.config`, which T1.3 deliberately left out — and note T1.3
      finding 6: the include and the file must land together.
      CARRIED FROM T1.4 — the whitelist SOURCING question is closed; the
      files are vendored at
      `assets/whitelist/10x_V{2,3,4}_barcode_whitelist.txt.gz` (R5.6).
      Two things remain here:
      (i)  extend the chemistry mapping to a THIRD column so one table
           gives the simpleaf chemistry, the QCatch chemistry and the
           whitelist filename. They must not be settable independently.
      (ii) PROVE, by running it in the pinned simpleaf container, whether
           `--unfiltered-pl` accepts a GZIPPED file. It is undocumented;
           see R5.6. If it does not, decompress in the module. Do not
           infer this from the fact that nf-core/scrnaseq passes `.txt.gz`.
      fastp args must not touch R1; comment the reason inline.
      Implements R3.1-R3.4, R5.1-R5.4, R6.4. [Opus, High]
      CARRIED FROM T2.1:
      - The chemistry mapping table is now VERIFIED against
        `simpleaf chemistry lookup` and `qcatch --help` in the pinned
        containers (T2.1 findings 3 and 5). Write it here, three columns as
        planned, and map `params.qcatch_n_partitions` to QCatch's
        `--n_partitions`. T2.1 deliberately did NOT create
        `conf/modules.config`: this task owns the file and its
        `includeConfig` line, and they must land together (T1.3 finding 6).
      - The v2 UMI is 10 bp, the v3/v4 UMI 12 bp (T2.1 finding 4). Nothing
        in `ext.args` sets a geometry — simpleaf takes it from the
        registered chemistry name — but do not restate "16+12" in a comment.
      - `simpleaf/quant` has no `ulimit` of its own; only `index` does.
        `scratch = true` still belongs on `SIMPLEAF_INDEX` only (R4.5).
      - VERIFIED, so it needs no probe here: fastp 1.3.6 accepts
        `--disable_adapter_trimming --disable_quality_filtering
        --disable_length_filtering` alongside the `--detect_adapter_for_pe`
        the module hard-codes into its paired-end branch. Exit 0, JSON and
        HTML both written. The hard-coded flag is inert in report-only mode.
      DELIVERED: `conf/modules.config` plus its `includeConfig` line in
      `nextflow.config`, and `mapper_backend` removed from all three params
      files. Lint clean at 19 files; three-way params parity holds at 36
      keys. `.claude/reference/modules.config.snippet` was NOT copied — it
      is stale in five separate ways, listed in finding 6 below.
      CARRIED FROM T2.2 — findings, each produced by running:
      1. **`--use-piscem` DOES NOT EXIST in simpleaf 0.30.0**, the only
         simpleaf this pipeline runs. `simpleaf index --help` has no mapper
         flag at all: piscem is the sole backend, under a "Piscem Index
         Options" heading. R4.2, R4.2a, R5.1, design section 7 and
         "CLAUDE.md section 4.2" all showed it and are corrected.
         `params.mapper_backend` is REMOVED from `params/default.yaml`,
         `params/demo.yaml` and `params/test.yaml` — a param that reaches
         no command line misrepresents what ran. T5.1: it must not appear
         in `nextflow_schema.json`.
      2. `--rlen` is OURS to pass and previously reached nothing. It sets
         roers' intron flank length and must track the cDNA read length, so
         `ext.args` passes `--rlen ${params.r2_read_length}`. Verified by
         building a real index with `--rlen 90`.
      3. `--ram-limit-gib` (new at this version) caps SSHash's external
         minimizer sort and DEFAULTS TO 8 GiB — larger than the `test`
         profile's entire 6 GB `resourceLimits` allowance. `ext.args`
         derives it from `task.memory` at 75%, floored at 1, so it tracks
         the profile instead of the host. Verified accepted by running
         `simpleaf index --rlen 90 --ram-limit-gib 2`.
      4. **THE GZIPPED WHITELIST WORKS — T2.2 item (ii) is CLOSED**, and no
         decompression step may be added. Two independent proofs, in R5.6:
         alevin-fry v0.18.1 reads the unfiltered permit list through
         `niffler::from_path` with the `gz` feature; and a synthetic
         fixture quantified with the gzipped vendored whitelist and with
         the same file decompressed gave the same 3 barcodes and
         BYTE-IDENTICAL `quants_mat.mtx`. The second proof is the one that
         matters: it rules out "accepted but silently garbage".
      5. **`--anndata-out` writes `af_quant/alevin/quants.h5ad`, not
         `af_quant/quants.h5ad`.** Found by running a real quant. R5.3a,
         "CLAUDE.md section 4.2" and the design dataflow all had the wrong
         path. R5.4a's `X` = spliced + unspliced + ambiguous convention was
         confirmed on that same file, along with the full obs/var/uns
         layout and a reproduction of the T1.4 spurious `None` layer key —
         all recorded in R5.4a for T3.1.
      6. `.claude/reference/modules.config.snippet` MUST NOT be copied. It
         is stale in five ways, each of which would break the run: a
         top-level `def` (hard parse error, T1.3 finding 1); `beforeScript`
         for `ulimit` (the module does it); `params.fastp_args` (dropped at
         T1.2); `--use-piscem` (gone, finding 1); and repeats of
         `--anndata-out`, `--unfiltered-pl`, `--resolution`,
         `--export_summary_table` and `--save_filtered_h5ad`, every one of
         which the modules hard-code or take as an input, so passing them
         again is a runtime error. The reference directory is a comparison
         target, not a source.
      7. `SIMPLEAF_QUANT` gets NO `ext.args`, deliberately: every argument
         it needs is either hard-coded in the module or arrives as a module
         input. `--min-reads` stays at simpleaf's 10, decided rather than
         inherited (R5.2); T2.5 may lower it with a measured barcode count.
      8. No `ext.prefix` is set anywhere. Every module already defaults to
         `${meta.id}`; setting it would duplicate an existing value.
      9. `params.multiqc_title` is interpolated at config-parse time, so an
         unguarded `--title '${params.multiqc_title}'` puts the literal
         string `null` on the report of anyone who forgets `-params-file`
         — params files are not auto-loaded. The `ext.args` is guarded with
         a ternary; verified by resolving the config without a params file
         and seeing `args = ''`.
      10. `nextflow config` in this version does NOT accept `-params-file`
          ("Unknown option"). To check how a params file interacts with a
          config-declared param, use a throwaway root `main.nf` and
          `nextflow run` — the T1.3 and T1.5 technique. Done here, and
          removed: `chemistry_map` SURVIVES a `-params-file` that never
          mentions it, so T1.5's "`-params-file` outranks config params"
          applies per key, not to the whole params map.
      11. `nextflow run` warns "There's no process matching config
          selector: X" for every `withName` that matches nothing, and lint
          does NOT catch it. Right now every selector warns, because no
          process exists yet. T4.2 must confirm the warnings are GONE:
          that message is how a mistyped process name shows up.
- [x] **T2.3 R1 assertion.** `modules/local/assert_r1_length/main.nf`.
      Implements R1.5. [Opus, Medium]
      DELIVERED: `ASSERT_R1_LENGTH` calling `bin/assert_r1_length.py`,
      plus an `r1_len` column added to `params.chemistry_map` in
      `conf/modules.config` (26 / 28 / 28), so the required length shares
      the chemistry's single source of truth and the module reads no params
      itself. Lint clean.
      Input shape: `[ meta, path(seqkit_stats_tsv), val(r1_filename),
      val(expected_r1_length) ]`. Output: `[ meta, path(tsv) ]`, one row —
      sample, chemistry, required, observed, min, max, num_reads,
      uniform_length, result — plus `python` to the versions topic.
      THE CHECK IS "AT LEAST", not "equal to" (R1.5 as amended): it applies
      simpleaf's own per-chemistry leniency, because the registered
      geometries end in `x:` and discard whatever follows the barcode and
      UMI. Too short is a hard failure; longer passes with a note.
      VERIFIED by running all five cases against real
      `seqkit stats --all` TSVs, in the T1.4 image:
      26 bp/`10XV2` pass; **28 bp/`10XV2` pass** — the case the earlier
      strict rule failed wrongly; 26 bp/`10XV3` fail, "R1 is TOO SHORT";
      24 bp/`10XV2` fail; ragged 24-26 bp/`10XV2` pass on the median with
      `uniform_length = no`. Plus the wiring-error path and `-stub-run`.
      CARRIED FROM T2.3:
      1. **NO script may be embedded in a `.nf` file** — no heredoc, no
         python shebang. Scripts are executable files in `bin/`, called
         from a bash script block; `bin/` is on `PATH` for every task. Now
         a standing rule in "CLAUDE.md section 3.2".
      2. **An `eval` output REQUIRES a bash script block.** A
         `#!/usr/bin/env python` script block fails at runtime with
         "Process output of type 'eval' is only allowed with Bash process
         scripts". Since every module reports versions through the topic,
         which is an eval, no local module may use a non-bash shebang.
      3. **`-stub-run` STILL EXECUTES THE CONTAINER** for eval outputs. The
         stub run reported `python 3.14.7`, the image's python, not the
         host's 3.14.4. A stub run is not image-free.
      4. A `bin/` script is only on `PATH` when the REPO is the Nextflow
         project directory. A throwaway probe workflow placed elsewhere
         fails with exit 127; put the probe at the repo root instead, as
         T1.3 and T1.5 did, and delete it afterwards.
      5. `seqkit stats --all` emits `Q1`, `Q2`, `Q3` LENGTH quartiles; `Q2`
         is the median and is what "modal length" is implemented as. There
         is no mode column. Full header:
         `file format type num_seqs sum_len min_len avg_len max_len Q1 Q2
         Q3 sum_gap N50 N50_num Q20(%) Q30(%) AvgQual GC(%) sum_n`.
- [ ] **T2.4 Reference subworkflow.**
      `subworkflows/local/prepare_splici_reference.nf`. Implements R4.
      `bin/download_reference.sh` ALREADY EXISTS and fetches the 10x
      GRCh38 2024-A package. Do not rewrite it; read it so the subworkflow
      consumes the same paths (`<ref>/fasta/genome.fa`,
      `<ref>/genes/genes.gtf.gz`) and publishes an index referenced as
      `<dir>/index` (`.gitignore`'d). `ALEVIN_FRY_HOME` and the open-file
      limit need no handling here: the module does both in-script
      (R4.5, R5.7, and T2.1 note (a)). [Opus, High]
      CARRIED FROM T2.1: `SIMPLEAF_INDEX` takes FOUR input tuples — genome
      fasta+gtf, transcript fasta, probe csv, feature csv — and picks the
      splici path only when fasta and gtf are both present and the other
      three are empty. It emits `${prefix}/index` plus an OPTIONAL
      `${prefix}/ref` and `${prefix}/ref/{t2g,t2g_3col}.tsv`; the t2g file
      is what `SIMPLEAF_QUANT` wants as `txp2gene`, so carry it through.
      A real index build needs the `conf/containers.config` override to be
      in place (T2.1 findings 6 and 6a); without it `simpleaf set-paths`
      aborts before any work happens. With it, the image is simpleaf
      0.30.0 / alevin-fry 0.18.1 / piscem 0.23.0.
      CARRIED FROM T2.2: `conf/modules.config` already supplies `--rlen`,
      `--ram-limit-gib` and `scratch = true`, so this subworkflow adds no
      arguments — it only wires channels. Do NOT pass `--use-piscem`; it
      does not exist at this version (T2.2 finding 1). A small index built
      here from a genome FASTA plus a GTF produced `<dir>/index/` (matching
      R4.2b), `<dir>/ref/t2g_3col.tsv` for `SIMPLEAF_QUANT`'s `txp2gene`,
      and a `gene_id_to_name.tsv` in both — see R6.4a before T2.5.
- [ ] **T2.5 QCatch wiring.** The nf-core module is named `qcatch`.
      Inspect it with `nf-core modules info qcatch` before wiring, so the
      input/output channel shapes are taken from the module rather than
      assumed. Implements R6.1-R6.6. Do NOT pass `--remove_doublets` or
      `--visualize_doublets`. [Opus, High]
      CARRIED FROM T1.5: this task settles the known ceiling on the
      `test` profile — whether QCATCH can run at all on the small nf-core
      fixture, which nf-core itself skips it on. See R15.3. If it cannot,
      say so plainly and record what `-profile test,docker` does cover;
      do not weaken the pipeline to make the test pass.
      WIRE `--cell_calling` HERE (R6.7). The parameter already exists in
      all three params YAMLs; `params/test.yaml` sets `"threshold"`.
      (i)   Branch the subworkflow on `params.cell_calling`:
            `qcatch` -> `SIMPLEAF_QUANT -> QCATCH -> SCANPY_CELL_QC`;
            `threshold` -> `SIMPLEAF_QUANT -> SCANPY_CELL_QC` directly.
            No new module and no new science param: the crude cut is the
            `min_genes`/`min_cells_per_gene` that `SCANPY_CELL_QC` already
            applies under R7.2.
      (ii)  An unrecognised value must fail fast with a readable message,
            not fall through to a default.
      (iii) MEASURE FIRST, then conclude: count the rows in the fixture's
            raw matrix. Under `--unfiltered-pl` only barcodes seen at
            least `--min-reads` times become rows and simpleaf passes 10
            (R5.2), so the matrix may already be near-empty and the
            `threshold` path may not rescue the profile either. If so,
            lowering `--min-reads` via `ext.args` is the next lever.
            Report the actual barcode count either way.
      (iv)  Only `qcatch` may emit `*_filtered_matrix.h5ad` (R6.6).
      CARRIED FROM T2.1:
      - The module's input is `tuple val(meta), val(chemistry),
        path(quant_dir)` and it emits `report`, `filtered_h5ad` and
        `metrics_summary`, renamed to `${prefix}_qcatch_report.html`,
        `${prefix}_filtered_quants.h5ad` and
        `${prefix}_metrics_summary.csv`. The `_filtered_matrix` suffix
        R6.6 wants is therefore NOT what the module produces; see T2.1
        finding 12 and settle the renaming at T4.2.
      - QCatch reaches the NETWORK for a gene-id-to-name mapping unless
        `--gene_id2name_file` is given, and silently drops the
        mitochondrial plots when that fails (T2.1 finding 13). Decide here
        whether to derive the TSV from the GTF, as the vendored whitelists
        did for simpleaf.
      - Measuring (iii) needs a real `SIMPLEAF_QUANT` run, which is
        possible again now that `conf/containers.config` replaces the
        crashing alevin-fry build (T2.1 findings 6 and 6a). Report the
        actual barcode count; do not guess it.
      CARRIED FROM T2.2:
      - The chemistry comes from `params.chemistry_map[meta.chemistry]`,
        which gives `.simpleaf`, `.qcatch` and `.whitelist` together. Look
        it up ONCE per sample and pass the parts to their modules; never
        let a caller set two of the three independently.
      - An unmapped `chemistry` value must fail with a readable message
        naming the sample and the allowed keys, not with a null-pointer
        error deep inside a module.
      - `ext.args` for QCATCH is already written and adds `--n_partitions`
        only when `params.qcatch_n_partitions` is set. For a custom assay
        the chemistry VALUE passed to the module must then be null, which
        makes the module omit `--chemistry` entirely.
      - The raw matrix is at `af_quant/alevin/quants.h5ad` (T2.2 finding 5).
      CARRIED FROM T2.3 — two QCatch findings that change what this task
      builds. Both come from running the pinned containers against real
      simpleaf output; see R6.4a and R6.4b.
      - **R6.4a is CLOSED, and the repack keeps it closed.** The repacked
        file is still an h5ad carrying `var['gene_symbol']`, so the
        no-network branch is the one taken; `--gene_id2name_file` is not
        needed on this path. Confirmed at T2.6 by running stock qcatch on
        the repacked real file with `--network none`.
        The detail, in case the input ever changes: with h5ad input QCatch
        makes NO network call
        (it returns early when `var` has `gene_symbol`, which simpleaf's
        h5ad does). With MTX input it always calls the registry and,
        offline, dies with an unhandled `ConnectionError`, exit 1 — not
        the graceful degradation `qcatch --help` describes. Passing
        `--gene_id2name_file <staged_dir>/gene_id_to_name.tsv` fixes it,
        and simpleaf writes that file into the very directory QCATCH is
        handed. This is possible from `ext.args` because the path is
        INSIDE the staged input, which makes the staged directory's name
        load-bearing: stage it under a fixed name.
      - **R6.4b is RESOLVED at T2.6, so this task only WIRES it.**
        simpleaf's h5ad is Blosc-compressed and no QCatch image can read
        it; `REPACK_H5AD` rewrites it as gzip first. So the chain here is
        `SIMPLEAF_QUANT -> REPACK_H5AD -> QCATCH`, and QCATCH's
        `path(quant_dir)` input receives the repacked FILE, not the
        `af_quant` directory. Verified at T2.6 on the real pbmc1k data:
        stock qcatch 0.2.13 reads it with `--network none` and calls 1,176
        cells from 25,705 barcodes.
        Passing the DIRECTORY instead would put QCatch back on simpleaf's
        Blosc file, which it prefers over the mtx when both are present.
      - The QCATCH container is already overridden to
        `quay.io/biocontainers/qcatch:0.2.13--pyhdfd78af_0` in
        `conf/containers.config` — current release, and faster. It does
        not change any of the above.
      - A too-small matrix fails inside cell calling with
        `AssertionError: Invalid selection of 0-count barcodes!`
        (`cell_calling.py:298`). That is what item (iii) above will hit if
        the fixture is too sparse, and it is a crash, not a graceful skip
        — another reason the `threshold` escape hatch exists.

- [x] **T2.6 h5ad repack.** `bin/repack_h5ad_blosc_to_gzip.py` and
      `modules/local/repack_h5ad/main.nf`, plus `hdf5plugin` in
      `docker/Dockerfile`. Resolves R6.4b and R7.8. Not in the original
      plan: it exists because simpleaf 0.30.0 writes h5ad with the Blosc
      HDF5 filter and nothing downstream can read it. [Opus, High]
      DELIVERED: `REPACK_H5AD` rewrites `af_quant/alevin/quants.h5ad` as
      gzip — the one filter built into every HDF5 build — and names the
      result `<sample>_raw_matrix.h5ad`, which also satisfies R5.5.
      `docker/Dockerfile` gained `hdf5plugin=7.1.0=py314haeb3a3a_0` and
      `ENV HDF5_PLUGIN_PATH=/opt/env/lib/python3.14/site-packages/hdf5plugin/plugins`.
      Lint clean at 21 files.
      ROUTE CHOSEN BY THE USER over the alternative of building our own
      QCatch image: that would be a SECOND pipeline-authored image, which
      R14.3 does not allow, and our image needed `hdf5plugin` anyway for
      R7.8 — two image changes instead of one, for no run-time gain.
      IMAGE TAG STAYS `scrnaseq-lite:0.1.0`, on the user's instruction:
      still initial development, so the local image is deleted and rebuilt
      rather than versioned. **Any run after this change needs that
      rebuild.**
      VERIFIED end to end on REAL data — the user's
      `reference/pbmc1k_quant` matrix, 25,713 x 38,606, filters
      `{'32001': 20, 'gzip': 5}`:
      - the script repacked it in **2.5 s**, 23 MB -> 28.6 MB, verifying
        that every dataset name, shape and dtype survived and that the sums
        of `X` and all three layers were unchanged;
      - the MODULE ran it in the rebuilt image with no ephemeral install,
        real and `-stub-run`, reporting `python 3.14.7` and `h5py 3.16.0`
        to the versions topic;
      - stock `quay.io/biocontainers/qcatch:0.2.13--pyhdfd78af_0` then
        processed the module's own output with `--network none`, exit 0,
        writing `QCatch_report.html`, `summary_table.csv` and
        `filtered_quants.h5ad`, and calling **1,176 cells from 25,705
        barcodes** — the right order of magnitude for pbmc1k.
      CARRIED FROM T2.6:
      1. **`HDF5_PLUGIN_PATH` alone is enough.** Copying plugins into
         HDF5's compiled-in default directory is NOT needed, and neither is
         `import hdf5plugin` in code we control. Verified both ways against
         a real Blosc file.
      2. The repack is an h5py-level copy, deliberately NOT an AnnData
         round-trip: it preserves `uns`, dtypes and attributes exactly,
         needs no AnnData semantics, and sidesteps the spurious `None`
         layer key of T1.4. It copies large datasets in slices, so peak
         memory does not scale with a sample's non-zero count.
      3. It VERIFIES itself and fails loudly: dataset set, shapes and
         dtypes must match, the sums of `X` and every layer must match, and
         no non-gzip filter may survive into the output. This file feeds
         every downstream step, so a silent truncation here would surface
         as biology rather than as an error.
      4. `blosc` and `c-blosc2` were already in the image and are NOT the
         same thing as the filter plugin: those are the compression
         libraries, `hdf5plugin` is what teaches HDF5 to use them.
      5. Gzip level 4 was chosen as a size/time compromise; the file grows
         about 22% against Blosc. That is the price of being readable
         everywhere.

## Wave 3 — Scanpy (parallel)

**CARRIED FROM T2.3 — no embedded scripts.** Every python script is an
executable file in `bin/`, called from a bash script block. No heredocs,
no `#!/usr/bin/env python` script blocks. `bin/assert_r1_length.py` with
`modules/local/assert_r1_length/main.nf` is the worked example, and the
rule is in "CLAUDE.md section 3.2".

**SETTLED AT T2.6, was carried from T2.3 — the Blosc problem no longer
reaches this wave.** `REPACK_H5AD` rewrites simpleaf's Blosc h5ad as gzip
before anything else sees it, so BOTH paths into `SCANPY_CELL_QC` — QCatch's
output under `--cell_calling qcatch`, the repacked raw matrix under
`threshold` — are plain gzip. Read `<sample>_raw_matrix.h5ad` on the
threshold path, never `af_quant/alevin/quants.h5ad` directly. The T1.4 image
now carries `hdf5plugin` as well, so reading the Blosc original would work,
but nothing should rely on that.

**CARRIED FROM T2.3, applies to every module in this wave.** A process
that reports its version through the `versions` topic uses an `eval`
output, and Nextflow allows `eval` outputs ONLY on bash script blocks:
a `#!/usr/bin/env python` script block dies with "Process output of type
'eval' is only allowed with Bash process scripts". Every local module here
therefore keeps a bash script block that CALLS its `bin/*.py` script —
which is the layout these tasks already specify — and never switches to a
python shebang. `modules/local/assert_r1_length/main.nf` is the worked
example. Note also that `-stub-run` still runs the container, because the
eval is evaluated there.

- [ ] **T3.1 Cell QC.** `bin/scanpy_cell_qc.py` and its module.
      Implements R7.1-R7.4. Verify the QCatch h5ad layout empirically
      before writing the reader; do not assume the obs/var structure.
      [Opus, High]
      CARRIED FROM T1.4, applies to all of Wave 3: anndata 0.13.3 returns
      a spurious `None` among `adata.layers.keys()` — observed as
      `['ambiguous', 'counts', 'spliced', 'unspliced', None]` after an
      h5ad round-trip. Never blindly iterate or `sorted()` layer keys;
      filter `None` first. This bites the R5.4a `S+A` reconstruction
      directly.
      Also CARRIED FROM T1.2: the gene patterns are params
      (`mito_gene_prefix`, `ribo_gene_pattern`, `hb_gene_pattern`), so
      this script holds no hard-coded gene patterns.
      CARRIED FROM R6.7: this step receives EITHER QCatch's cell-called
      matrix or, under `--cell_calling threshold`, the raw quantifier
      matrix. Do not assume QCatch's obs/var columns are present — verify
      the layout of both inputs empirically. Write
      `adata.uns['cell_calling_method']` so provenance travels with the
      h5ad after it leaves the results directory, and make the
      barcode-rank plot (R7.3) work on whichever matrix arrives; it is
      more informative on the raw one.
- [ ] **T3.2 Doublets.** `bin/scanpy_detect_doublets.py` and its module.
      Implements R7.5-R7.7. Seeded. [Opus, Medium]
      CARRIED FROM T1.4: `sc.pp.scrublet` needs `scikit-image`, which is
      now in the image. It raises only at RUN time, so an import check
      will not catch a regression here.
- [ ] **T3.3 Normalise and cluster.** `bin/scanpy_normalise_cluster.py`
      and its module. Implements R8 and R9. Seeded throughout.
      [Opus, High]
      CARRIED FROM T1.4: `flavor='seurat_v3'` needs `scikit-misc` and
      `flavor='igraph'` needs `leidenalg`; both are in the image and both
      were exercised in the T1.4 determinism run.
- [ ] **T3.4 Marker genes.** `bin/scanpy_marker_genes.py` and its module.
      Implements R10. No filename, header, comment or plot title may say
      DE or DGE. [Opus, Medium]
- [ ] **T3.5 Versions.** `bin/collect_versions.py` and its module ->
      `versions.tsv`. Implements R12.1. [Opus, Medium]
      SETTLED AT T2.1, superseding the T1.5 note that said the two
      conventions are mixed: ALL SIX installed modules emit `eval(...)`
      into `topic: versions` and no `versions.yml` exists anywhere in
      `modules/nf-core/`. `COLLECT_VERSIONS` reads the topic only, and all
      local modules use the topic convention too. Three consequences, all
      from T2.1:
      - `MULTIQC` is the ONE exception and must NOT be wired into the
        topic: it consumes the topic, and feeding it back in hangs the run
        forever (the module says so in a comment). Merge its plain
        `emit: versions` channel separately, and order `COLLECT_VERSIONS`
        before `MULTIQC`.
      - A failing eval yields an EMPTY string and the task still
        SUCCEEDS (finding 7). Treat an empty version as an error or an
        explicit `unknown`; never write a blank cell.
      - `qcatch`'s eval yields `version 0.2.12`, not `0.2.12` (finding 8).
        Strip a leading `version ` token. The vendored module is not edited.
      CARRIED FROM T1.4: do not try to version-report `ps`; the
      `procps-ng` build in the image reports `UNKNOWN`.
- [ ] **T3.6 Run manifest.** `bin/write_run_manifest.py` and its module ->
      `params_resolved.yaml` and `run_manifest.yaml`. Implements R2.5.
      [Opus, Medium]
      CARRIED FROM R6.7: `run_manifest.yaml` must record `cell_calling`.
      A run whose cells were called by a crude threshold has to be
      distinguishable from one that used the ambient model, from the
      results directory alone.

## Wave 4 — Assembly (sequential)

- [ ] **T4.1 Subworkflows.** Explicitly named:
      `READ_QC`, `PREPARE_SPLICI_REFERENCE`, `QUANTIFY_AND_CALL_CELLS`,
      `CELL_QC_AND_DOUBLET_FILTERING`, `NORMALISE_CLUSTER_AND_MARKERS`,
      `REPORTING`. Annotate every channel operation with its in/out shape.
      [Opus, High]
      CARRIED FROM T2.2:
      - `READ_QC` passes `discard_trimmed_pass = true` to `FASTP`, so fastp
        writes no read files at all. Verified at T2.2 that this works
        despite the module putting a bare `true` token on the command line;
        it is the module's documented report-only mode and it satisfies
        R3.3 more strongly than merely ignoring the outputs.
      - The chemistry, the QCatch chemistry and the whitelist FILENAME all
        come from `params.chemistry_map`; the whitelist path is built
        against `assets/whitelist/` and staged as `cb_list`. Pass the
        gzipped file as-is — R5.6 is closed, and decompressing would just
        add a process.
      - `SIMPLEAF_QUANT` needs `cell_filter = 'unfiltered-pl'`, the
        `t2g_3col.tsv` from `SIMPLEAF_INDEX` as `txp2gene`, and
        `params.umi_resolution` as the `val resolution` input. None of
        those are `ext.args`.
      CARRIED FROM T2.3 — wiring `ASSERT_R1_LENGTH` in `READ_QC`. It takes
      `[ meta, path(seqkit_stats_tsv), val(r1_filename),
      val(expected_r1_length) ]`:
      - `r1_filename` must be the STAGED name of the R1 file as
        `SEQKIT_STATS` saw it — `reads[0].name`, not a path. The module
        matches the stats row on that exact string and fails loudly if it
        finds no match, so a mis-wire is caught rather than silently
        skipped.
      - `expected_r1_length` is `params.chemistry_map[meta.chemistry].r1_len`.
      - `SEQKIT_STATS` runs once per sample over BOTH mates, so its TSV has
        an R1 row and an R2 row; the module needs the R2 row to exist and
        ignores it.
- [ ] **T4.2 Entry workflow.** `workflows/scrnaseq_lite.nf` and `main.nf`
      with `publish:` and `output` blocks. Implements R13 and follows
      "Design section 5: Publishing" and "Design section 6: main.nf
      structure". Highest-consequence file in the repo. [Opus, Max]
      This task is what first makes `nextflow run .` possible; the
      after-every-task run check in the Rules above only becomes real
      here. Obey T1.3 finding 9 (lowercase `channel`, `_`-prefixed unused
      output-closure parameters) and finding 8 (every process needs both
      `script:` and `stub:`).
      CARRIED FROM R6.7: under `--cell_calling threshold` QCATCH does not
      run, so the `qcatch` entry in `publish:` receives an EMPTY channel
      while still needing its matching `output` entry, which carries
      `index { ... header true }`. VERIFY by running that an empty channel
      publishes cleanly rather than erroring or leaving a header-only CSV
      behind. `-profile test` exercises this path, so it is not a corner
      case.
      CARRIED FROM T2.1: the provenance suffixes R6.6 and "CLAUDE.md
      section 4.2" mandate are NOT what the vendored modules produce —
      QCatch writes `*_filtered_quants.h5ad` and simpleaf writes
      `af_quant/alevin/quants.h5ad` inside a directory — note the `alevin/`
      level, corrected at T2.2. Apply `*_raw_matrix.h5ad`
      and `*_filtered_matrix.h5ad` in the `output` block's `path`
      directive, one rule for both, rather than editing the modules.
      Also: `MULTIQC` must be ordered AFTER `COLLECT_VERSIONS` and must
      never be wired into the `versions` topic, which deadlocks the run
      (T2.1 finding 2).
- [ ] **T4.3 MultiQC.** `assets/multiqc_config.yaml` plus custom-content
      injection of the QCatch summary CSV and the Scanpy QC TSV.
      Implements R11.1. [Opus, High]
      CARRIED FROM R6.7: under `--cell_calling threshold` the QCatch
      summary CSV does not exist, and it is MultiQC's ONLY quantifier-QC
      source (design section 2). The report must then cover reads and the
      Scanpy QC TSV only — and must SAY so, with a visible banner naming
      the cell-calling method, rather than letting the section silently
      vanish as though quantifier QC had passed.
      CARRIED FROM T2.2: the report TITLE is already set, as
      `--title '${params.multiqc_title}'` in `conf/modules.config`
      `ext.args`. Do not also put a title in `assets/multiqc_config.yaml` —
      that would be the same value in two files.
- [ ] **T4.4 Analysis report.** `report/analysis_report.qmd` and its
      module. Implements R11.3 and R11.4. [Opus, Medium]
      CARRIED FROM T1.4 — the `.qmd` must NOT call
      `matplotlib.use("Agg")`. Under quarto's jupyter engine that silently
      suppresses every figure: the render SUCCEEDS, prints no warning, and
      the HTML simply contains zero `<img>` tags. Verified both ways in
      the T1.4 image. Rendering is otherwise confirmed working —
      `quarto render --to html --embed-resources` produced an embedded PNG
      and the executed cell's stdout.
      CARRIED FROM R6.7: the report carries a visible banner whenever
      `cell_calling != 'qcatch'`, stating that cells were called by a
      crude threshold with no ambient model and that the result is not
      publishable. Read it from `adata.uns['cell_calling_method']`, which
      travels with the h5ad, rather than from the param.

## Wave 5 — Data, run, validate

- [ ] **T5.1 Schemas2.** `nextflow_schema.json` via
      `nf-core pipeline schema build` (pushed back here, as it requires a
      pipeline with `main.nf`). Implements R1.2-R1.4. [Opus, High]
      CARRIED FROM T1.3: mark the derived param `run_output_dir` as
      `hidden`.
      CARRIED FROM T2.2: `chemistry_map` is the second derived param and
      must also be `hidden` — it is a lookup table, not a user knob, and
      `nf-core pipeline schema build` will offer to add it. `mapper_backend`
      must NOT appear at all; it was removed from every params file.
      CARRIED FROM R6.7: `cell_calling` is an enum,
      `["qcatch", "threshold"]`, default `qcatch`. Its schema description
      must state that `threshold` performs no ambient modelling and is not
      acceptable for real data — the schema is what `nf-core pipeline
      schema build` surfaces to users, so the warning belongs there too.
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
      This is what makes `assets/samplesheet_demo.csv` validate; until
      then its failure is expected (T1.5).
- [ ] **T5.3 Green run.** Full `-stub-run`, then `-profile test,docker`,
      then the real `-profile demo,docker` two-sample run. Fix to green.
      [Opus, High]
      CARRIED FROM T2.3: `-stub-run` is NOT image-free. Eval outputs are
      evaluated in the container even when the script body is stubbed, so
      the first stub run pulls every image and needs `scrnaseq-lite:0.1.0`
      to exist. Budget for that instead of treating the stub pass as a
      no-dependency smoke test.
      CARRIED FROM T2.1 — the SIGILL blocker is FIXED, but check the fix is
      still in place before blaming anything else. `conf/containers.config`
      overrides the two simpleaf processes onto
      `quay.io/biocontainers/simpleaf:0.30.0--hd612981_0` (alevin-fry
      0.18.1, the portable build) because the module's pinned alevin-fry
      0.15.0 is a `target-cpu=native` bioconda build that SIGILLs on any
      CPU below the build worker's; see T2.1 findings 6 and 6a. If a run
      dies with "Error running alevin-fry ... signal: 4 (SIGILL)", that
      override went missing. Two things this task must still confirm on
      REAL data. T2.2 has since settled both on a synthetic fixture:
      R5.4a's `X` = spliced + unspliced + ambiguous convention HOLDS at
      simpleaf 0.30.0, and `--anndata-out` writes
      `af_quant/alevin/quants.h5ad` — one level deeper than every spec
      previously said. What remains here is confirming the same on real
      data at scale. Do not call the pipeline green on stubs alone.
- [ ] **T5.4 Validation.** Run `/validate-pipeline`. Produce
      `docs/validation-report.md` (gitignored). [Opus, Max]
- [ ] **T5.5 Documentation.** `README.md`, `docs/usage.md`,
      `docs/output.md`. Must state the `--run` / `-resume` requirement
      (R2.4), the difference between the `test` and `demo` profiles
      (R15.3), the two subsampling strategies and why they differ (R15.5),
      and the dataset licence and attribution (R15.9). Add one line to
      `README.md` noting that the pipeline was built spec-driven
      development and that the requirements, design and task documents are
      in `.claude/specs/`. No private or environment-specific information
      (N4). [Opus, Medium]
      CARRIED FROM T1.4: word the reproducibility claim as "PCA and UMAP
      coordinates may differ across BLAS thread counts" (R12.4). Do not
      promise cluster stability across thread counts.
- [ ] **T5.6 Changelog.** Finalise `CHANGELOG.md` 0.1.0 with today's date.
      A single "Initial implementation" entry only. [Opus, Low]
- [ ] **T5.7 CITATIONS.md.** Include a `CITATIONS.md` with the relevant
      citations for this pipeline, from
      `https://raw.githubusercontent.com/nf-core/scrnaseq/refs/tags/4.2.0/CITATIONS.md`.
      Keep only the tools this pipeline actually runs. [Opus, Low]
