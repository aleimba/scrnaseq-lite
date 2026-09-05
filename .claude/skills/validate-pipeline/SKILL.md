Perform a full validation pass on this pipeline. Produce
`docs/validation-report.md` with one section per check and a
PASS / FAIL / UNCERTAIN verdict for each. Be conservative: mark UNCERTAIN
rather than guessing.

## A. Syntax

1. Run `nextflow lint .` and report every finding.
2. Run `nextflow run . -profile test -stub-run`. It must complete cleanly.
3. Run `nextflow run . -profile test,docker`. It must complete cleanly.
4. Grep for bare `it` inside closures in pipeline-authored `.nf` files.
   Report file:line for each hit. Report `modules/nf-core/` separately.
5. Grep for `publishDir` outside `modules/nf-core/`.
6. Grep for `nextflow.preview.output` and `nextflow.enable.dsl`.
7. Grep for `nf-validation`.

## B. Dataflow

8. Produce the DAG. Report orphan channels and any process whose output is
   neither consumed nor published.
9. Confirm every name in `publish:` has a matching `output` entry, and
   vice versa.
10. Confirm every output `path` closure's parameters match the published
    channel's actual shape.
11. Confirm every `include` resolves to an existing file and process.
12. Confirm each channel operation's shape comment matches reality.
13. Confirm every process reports tool versions into the `versions` topic
    — settled at T2.1: all six installed nf-core modules use the topic and
    none writes a `versions.yml`, so every local module must match. Confirm
    `COLLECT_VERSIONS` receives all of them, that it runs BEFORE `MULTIQC`,
    and that `MULTIQC` is NOT wired into the topic it consumes (doing so
    hangs the run). Confirm no row of `versions.tsv` has an empty version
    (a failing `eval` yields an empty string without failing the task) and
    that qcatch's row reads `0.2.12`, not `version 0.2.12`.
14. Flag any `combine` that could produce a cartesian explosion.
15. Confirm `outputDir` does not feed into any task input.

## C. Scientific correctness

16. Quote the resolved fastp command. Confirm it cannot alter R1.
17. Confirm fastp read outputs are not consumed downstream.
18. Confirm `simpleaf quant` uses `--unfiltered-pl` and NO other
    permit-list mode. `--expect-cells`, `--explicit-pl`, `--forced-cells`
    and `--knee` are mutually exclusive with it and must not appear.
18a. Confirm the samplesheet `expected_cells` value never reaches
    `simpleaf quant`. It is metadata only.
18b. Confirm `--resolution` is passed and matches `params.umi_resolution`.
18b2. Confirm `--anndata-out` is passed, so QCatch receives h5ad rather
    than mtx.
18b3. Confirm the v4 chemistry maps to `10xv4-3p`, not `10xv4` — verified
    at T2.1 against `simpleaf chemistry lookup`, which registers no
    `10xv4`. Also confirm the QCatch column resolves to `10X_3p_v{2,3,4}`.
18b4. Confirm `ALEVIN_FRY_HOME` is set and the open-file limit is raised
    for the simpleaf processes. Both come from the vendored modules'
    script blocks (`index` sets `ulimit -n 2048`, `quant` does not need
    it); the pipeline adds only `scratch = true` on `SIMPLEAF_INDEX`.
18b5. State whether the barcode whitelist is pre-seeded in the container
    or fetched at run time, and confirm the choice is documented.
18c. Confirm `--use-piscem` is passed at INDEX time and not duplicated at
    quant time.
18d. Confirm the index is referenced as `<dir>/index/`.
19. Confirm `SIMPLEAF_INDEX` sets `scratch = true`.
20. Confirm the samplesheet `chemistry` value drives BOTH the simpleaf
    chemistry and the QCatch chemistry through a single mapping, and that
    the two cannot be set independently.
21. Confirm QCatch runs with `--skip_umap_tsne --export_summary_table
    --save_filtered_h5ad` and WITHOUT `--remove_doublets` or
    `--visualize_doublets`.
22. Confirm raw counts are in `layers['counts']` before normalisation and
    that `seurat_v3` HVG receives raw counts.
23. Confirm `params.seed` reaches scrublet, PCA, UMAP and Leiden.
24. Confirm doublet removal precedes HVG and PCA.
25. Confirm the matrix orientation is cells x genes, and that
    `params.count_layer` is reconstructed from the NAMED U/S/A layers
    rather than read from `X`. simpleaf writes `X` as U+S+A, so reading
    `X` for an `S+A` request is a silent correctness bug.
26. Grep all code, comments, filenames, plot titles and report templates
    for "DE", "DGE" and "differential expression". Any hit near
    marker-gene code is a FAIL.
27. Confirm the R1 length assertion enforces 26 bp for `10XV2` and 28 bp
    for `10XV3` and `10XV4`.

## D. Outputs

28. Confirm each process publishes to its own lower-cased name.
29. Confirm `final_results/` receives copies of the set in R13.2.
30. Confirm `mode: 'symlink'` appears nowhere.
31. Confirm matrices carry the `_raw_matrix` / `_filtered_matrix` suffix.

## E. Repository layout

32a. Confirm the repository root contains no planning or agent documents
    other than `CLAUDE.md`, and no top-level `specs/` or `reference/`
    directory.
32b. Confirm nothing under `.claude/` has been modified.
32c. Run `nf-core pipelines lint` and report every finding. If it objects
    to the layout, report the specific check name rather than changing
    the layout.

## F. Portability and privacy

32. Confirm no absolute local paths anywhere.
33. Confirm every process has a container directive.
34. Confirm `conf/aws.config` contains no account IDs, bucket names or
    ARNs, and that `conf/aws.local.config` is gitignored.
35. Grep every tracked file for personal names, email addresses,
    institution names, home directory paths, AWS account IDs, bucket names
    and ARNs. Any hit is a FAIL.

## G. Documentation

36. Confirm `README.md` and `docs/usage.md` state that `-resume` requires
    `--run <name>`.
37. Confirm the docs state that `-profile test` results are not
    biologically meaningful.
38. Confirm both datasets are attributed by name, credited to 10x
    Genomics, and linked to their dataset pages and to CC BY 4.0.

## H. Report

State clearly which checks you could not perform, and why.
