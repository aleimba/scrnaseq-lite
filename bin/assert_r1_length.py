#!/usr/bin/env python
"""Assert that observed R1 length is compatible with the declared chemistry.

Implements requirement R1.5. R1 is not a biological read: it is a fixed
layout of cell barcode and UMI at known offsets, and the samplesheet
`chemistry` column declares that layout. An R1 that is too SHORT for the
declared chemistry means the chemistry is wrong, the mates are swapped, or
R1 has been trimmed -- and each of those yields a plausible but WRONG
matrix rather than an error further downstream. So it is a hard failure.

Length is read from `seqkit stats --tabular --all` output. seqkit has no
mode column, so the MEDIAN length (its `Q2` column) stands in for the modal
length: for 10x R1 the two are the same number, because the read is a fixed
length by construction, and the median is robust to a minority of ragged
reads in a way `max_len` is not.

The comparison uses the SAME leniency as simpleaf, per chemistry, rather
than a stricter rule of its own. simpleaf's registered geometries end in
`x:` -- `1{b[16]u[10]x:}` for 10xv2 and `1{b[16]u[12]x:}` for 10xv3 and
10xv4-3p -- which means "barcode, UMI, then discard whatever remains". A
longer R1 is therefore perfectly usable, and a v2 library sequenced with a
28-cycle R1 is a real and valid case. Only a read too short to contain the
barcode and UMI is an error, so the test is `observed >= expected`.
"""

import argparse
import csv
import sys


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stats", required=True, help="seqkit stats --all TSV")
    parser.add_argument("--r1-filename", required=True, help="staged R1 filename, as seqkit saw it")
    parser.add_argument("--expected-length", required=True, type=int, help="minimum R1 length for the chemistry")
    parser.add_argument("--sample", required=True, help="sample id, for messages")
    parser.add_argument("--chemistry", required=True, help="samplesheet chemistry, for messages")
    parser.add_argument("--output", required=True, help="TSV to write")
    return parser.parse_args(argv)


def read_r1_row(stats_path, r1_filename, sample):
    with open(stats_path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for row in rows:
        if row["file"] == r1_filename:
            return row

    # Not a data problem: the workflow passed a name seqkit never saw.
    seen = ", ".join(sorted(row["file"] for row in rows)) or "nothing"
    sys.exit(
        f"ASSERT_R1_LENGTH [{sample}]: no seqkit stats row for R1 file "
        f"'{r1_filename}'. The stats file lists: {seen}. This is a pipeline "
        f"wiring error, not a data problem."
    )


def main(argv=None):
    args = parse_args(argv)

    row = read_r1_row(args.stats, args.r1_filename, args.sample)
    observed = int(float(row["Q2"]))
    min_len = int(float(row["min_len"]))
    max_len = int(float(row["max_len"]))
    n_seqs = int(float(row["num_seqs"]))

    ragged = min_len != max_len
    passed = observed >= args.expected_length
    longer = observed > args.expected_length

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "chemistry",
                "required_r1_length",
                "observed_r1_length",
                "min_length",
                "max_length",
                "num_reads",
                "uniform_length",
                "result",
            ]
        )
        writer.writerow(
            [
                args.sample,
                args.chemistry,
                args.expected_length,
                observed,
                min_len,
                max_len,
                n_seqs,
                "no" if ragged else "yes",
                "pass" if passed else "fail",
            ]
        )

    if not passed:
        sys.exit(
            f"ASSERT_R1_LENGTH [{args.sample}]: R1 is TOO SHORT for the "
            f"declared chemistry.\n"
            f"  samplesheet chemistry : {args.chemistry}\n"
            f"  required R1 length    : at least {args.expected_length} bp "
            f"(cell barcode + UMI)\n"
            f"  observed (median)     : {observed} bp   "
            f"[min {min_len}, max {max_len}, over {n_seqs} reads]\n"
            f"  R1 file               : {args.r1_filename}\n"
            f"An R1 shorter than the barcode and UMI cannot be decoded. Check "
            f"the chemistry column for this sample, confirm R1 and R2 are not "
            f"swapped, and confirm R1 has not been trimmed. R1 must never be "
            f"trimmed or length-filtered: that shifts the barcode and UMI "
            f"offsets."
        )

    if longer:
        # Valid and common: a v2 library sequenced with a 28-cycle R1.
        # simpleaf's geometry discards the surplus bases.
        print(
            f"ASSERT_R1_LENGTH [{args.sample}]: R1 median length {observed} bp "
            f"exceeds the {args.expected_length} bp required by {args.chemistry}. "
            f"The surplus bases are discarded by simpleaf's geometry. Proceeding.",
            file=sys.stderr,
        )

    if ragged:
        print(
            f"ASSERT_R1_LENGTH [{args.sample}]: R1 reads are not all the same "
            f"length [min {min_len}, max {max_len}]. The median is used. "
            f"Proceeding.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
