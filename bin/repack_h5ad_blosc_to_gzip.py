#!/usr/bin/env python
"""Rewrite an h5ad's datasets with gzip compression instead of Blosc.

WHY THIS EXISTS. `simpleaf quant --anndata-out` writes `quants.h5ad` using
the Blosc HDF5 filter (id 32001) once the matrix is above a small size
threshold. Reading that file needs the Blosc filter plugin, which the QCatch
container does not ship -- QCatch fails while parsing its arguments, before
any analysis, with "Can't synchronously read data (can't open directory
.../hdf5/plugin)". The threshold is why this is easy to miss: a toy fixture
is written uncompressed and reads fine everywhere, while every real sample
is compressed and reads nowhere.

Gzip is the one filter HDF5 always has built in, so a gzip-compressed copy is
readable by any HDF5 build. QCatch itself writes its outputs with gzip.

This is a raw h5py copy, NOT an AnnData round-trip. It preserves `uns`,
dtypes, attributes and structure exactly, needs no AnnData semantics, and
avoids the spurious `None` layer key that anndata 0.13.3 returns.
"""

import argparse
import sys

import h5py
import numpy as np

# Copy large datasets in slices rather than reading them whole, so peak
# memory tracks this constant instead of a sample's non-zero count.
SLICE_ROWS = 1_000_000


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="h5ad written by simpleaf quant")
    parser.add_argument("--output", required=True, help="gzip-compressed h5ad to write")
    parser.add_argument(
        "--compression-level", type=int, default=4,
        help="gzip level; 4 is a deliberate size/time compromise (default: 4)",
    )
    return parser.parse_args(argv)


def copy_dataset(name, source, dest_group, key, level):
    """Recreate one dataset under dest_group, gzip-compressed where useful."""
    # Scalars, object dtypes (variable-length strings) and empty datasets are
    # copied verbatim: they cannot be chunked, and chunking is what
    # compression requires.
    if source.shape == () or source.dtype == object or source.size == 0:
        dest = dest_group.create_dataset(key, data=source[()])
        dest.attrs.update(source.attrs)
        return

    dest = dest_group.create_dataset(
        key,
        shape=source.shape,
        dtype=source.dtype,
        compression="gzip",
        compression_opts=level,
        chunks=source.chunks or True,
    )

    if source.shape[0] <= SLICE_ROWS:
        dest[...] = source[...]
    else:
        for start in range(0, source.shape[0], SLICE_ROWS):
            stop = min(start + SLICE_ROWS, source.shape[0])
            dest[start:stop] = source[start:stop]

    dest.attrs.update(source.attrs)


def copy_group(source_group, dest_group, level):
    for key, item in source_group.items():
        if isinstance(item, h5py.Group):
            child = dest_group.create_group(key)
            child.attrs.update(item.attrs)
            copy_group(item, child, level)
        else:
            copy_dataset(item.name, item, dest_group, key, level)


def collect_datasets(handle):
    """name -> (shape, dtype) for every dataset in the file."""
    found = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            found[name] = (obj.shape, obj.dtype)

    handle.visititems(visit)
    return found


def checksum(handle, name):
    """Sum of a numeric dataset, or None when it is not summable."""
    if name not in handle:
        return None
    dataset = handle[name]
    if not isinstance(dataset, h5py.Dataset) or dataset.dtype == object:
        return None
    if not np.issubdtype(dataset.dtype, np.number):
        return None
    total = 0.0
    for start in range(0, max(dataset.shape[0], 1), SLICE_ROWS):
        total += float(np.asarray(dataset[start:start + SLICE_ROWS]).sum())
    return total


def verify(source_path, dest_path):
    """Fail loudly unless the copy is faithful.

    This file feeds every downstream step, so a silently truncated matrix
    would surface as biology rather than as an error.
    """
    with h5py.File(source_path, "r") as source, h5py.File(dest_path, "r") as dest:
        src_sets = collect_datasets(source)
        dst_sets = collect_datasets(dest)

        if set(src_sets) != set(dst_sets):
            missing = sorted(set(src_sets) - set(dst_sets))
            extra = sorted(set(dst_sets) - set(src_sets))
            sys.exit(
                f"REPACK_H5AD: dataset set changed during repack.\n"
                f"  missing from output: {missing or 'none'}\n"
                f"  unexpected in output: {extra or 'none'}"
            )

        for name, (shape, dtype) in src_sets.items():
            if dst_sets[name] != (shape, dtype):
                sys.exit(
                    f"REPACK_H5AD: '{name}' changed shape or dtype: "
                    f"{shape}/{dtype} -> {dst_sets[name][0]}/{dst_sets[name][1]}"
                )

        # Sum the count data itself, not just its shape: matching shapes with
        # different values would be the worst possible outcome here.
        summable = ["X/data"] + [
            f"layers/{layer}/data" for layer in source.get("layers", {})
        ]
        for name in summable:
            src_total = checksum(source, name)
            dst_total = checksum(dest, name)
            if src_total != dst_total:
                sys.exit(
                    f"REPACK_H5AD: '{name}' sum changed during repack: "
                    f"{src_total} -> {dst_total}"
                )
            if src_total is not None:
                print(f"  verified {name}: sum {src_total:g}", file=sys.stderr)


def main(argv=None):
    args = parse_args(argv)

    with h5py.File(args.input, "r") as source, h5py.File(args.output, "w") as dest:
        dest.attrs.update(source.attrs)
        copy_group(source, dest, args.compression_level)

    verify(args.input, args.output)

    with h5py.File(args.output, "r") as dest:
        filters = set()

        def note(name, obj):
            if isinstance(obj, h5py.Dataset):
                filters.update(obj._filters or {})

        dest.visititems(note)

    leftover = sorted(f for f in filters if f != "gzip")
    if leftover:
        # Belt and braces: the point of this step is that nothing needing an
        # external plugin survives into the output.
        sys.exit(
            f"REPACK_H5AD: output still uses non-gzip HDF5 filters {leftover}. "
            f"Downstream containers without those plugins cannot read it."
        )

    print(f"REPACK_H5AD: wrote {args.output} (gzip level {args.compression_level})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
