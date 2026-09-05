# 10x cell-barcode whitelists

Vendored, gzipped 10x Genomics cell-barcode whitelists, one per chemistry
supported by `assets/schema_input.json`.

## Why these files are in the repository

`simpleaf quant` is run in `--unfiltered-pl` mode, which corrects observed
cell barcodes against the full whitelist for the declared chemistry. Given
the bare flag, simpleaf resolves that whitelist over the network at run
time and caches it under `ALEVIN_FRY_HOME`. That is a silent run-time
network dependency: it fails on an offline compute node and inside a
locked-down AWS VPC, and it fails late, in the middle of a run.

Vendoring the whitelists removes the dependency outright. The file is
staged as an ordinary channel input and passed explicitly as
`--unfiltered-pl <file>`, so no run ever reaches for the network.

The cost is 34 MB in the repository. That is the whole price of a pipeline
that runs unattended in a private subnet.

## Contents

| File | Chemistry (samplesheet) | Barcodes | Bytes |
| ---- | ----------------------- | -------- | ----- |
| `10x_V2_barcode_whitelist.txt.gz` | `10XV2` | 737,280 | 2,238,617 |
| `10x_V3_barcode_whitelist.txt.gz` | `10XV3` | 6,794,880 | 18,350,152 |
| `10x_V4_barcode_whitelist.txt.gz` | `10XV4` | 7,372,800 | 13,390,930 |

Every line is a bare 16 bp `[ACGT]` barcode. There is no header.

10x V1 is deliberately absent: `assets/schema_input.json` restricts the
samplesheet `chemistry` column to `10XV2`, `10XV3` and `10XV4`, so a V1
list could never be selected.

## Provenance

Retrieved 2026-09-04 from nf-core/scrnaseq at the pinned tag `4.2.0`:

```
https://raw.githubusercontent.com/nf-core/scrnaseq/4.2.0/assets/whitelist/<file>
```

Downloaded sizes were checked against the byte counts the GitHub contents
API reports for that tag; all three matched exactly.

```
4101687b6cbb947b8ace340c38eecf872a1a59f230eab23becacd038a46c6fb5  10x_V2_barcode_whitelist.txt.gz
6f5c08bd6c0c63e7cd4f62bc5ec47024baaf6b230524f313799b9f5431df7b37  10x_V3_barcode_whitelist.txt.gz
649c5cc14a2b449db16867e9772cc7a115ac971563a6a4a9097156ff8c84940a  10x_V4_barcode_whitelist.txt.gz
```

The lists themselves originate with 10x Genomics, which distributes them
inside the Cell Ranger software package (`737K-august-2016`,
`3M-february-2018` and `3M-3pgex-may-2023` respectively). They are
redistributed here as nf-core/scrnaseq redistributes them.

## Which file goes with which run

The mapping from the samplesheet `chemistry` column to the simpleaf
chemistry string, the QCatch chemistry string and the whitelist file is a
single table in `conf/modules.config`. Never select a whitelist by hand:
the three values must not be settable independently or they will drift
apart. See "CLAUDE.md section 4.3: Chemistry mapping — one source of
truth".
