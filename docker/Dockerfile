# scrnaseq-lite - image for the pipeline's LOCAL modules only.
#
# Build:
#   docker build -t scrnaseq-lite:0.1.0 docker/
#
# SCOPE. This pipeline uses per-module biocontainers: every vendored
# nf-core module keeps the container it pins and there is no global
# process.container override. This image therefore covers only what this
# repository itself authors -- the four Scanpy steps, the R1-length
# assertion, the Quarto analysis report, and the version and manifest
# scripts.
#
# Deliberately ABSENT, because each ships in its own biocontainer:
# simpleaf, alevin-fry, piscem, qcatch, fastp, seqkit, multiqc. In
# particular nothing simpleaf-related belongs here -- no ALEVIN_FRY_HOME,
# no `simpleaf set-paths`, no raised open-file limit, no barcode
# whitelist. The nf-core simpleaf modules do all of that in their own
# script blocks.

ARG MICROMAMBA_IMAGE=mambaorg/micromamba@sha256:ed1bc6280b1eadf303610c00befaf9ba52dbacbc988ae4d7d9765801465f129e

# --- Stage 1: solve and install into a self-contained prefix ----------------
#
# The base is pinned by digest, not by the mutable `2.9.0-debian13` tag, so
# a rebuild cannot silently move underneath the package pins below.
FROM ${MICROMAMBA_IMAGE} AS build

# /opt is root-owned in the base image, whose default user is mambauser.
USER root

# One install layer. Every package this pipeline names is pinned to an
# exact version=build triple, taken from a `micromamba create --dry-run`
# solve rather than written from memory. The full 220-package closure that
# this resolves to is recorded in docs/container.md.
#
# `create`, not `install`: micromamba refuses to install into a prefix
# that does not exist yet, with "Environment must first be created with
# micromamba create". It is still a single layer.
#
# Why each non-obvious package is here:
#   scikit-misc   - required by sc.pp.highly_variable_genes(flavor='seurat_v3')
#   scikit-image  - required by sc.pp.scrublet to pick its own threshold.
#                   Without it scrublet raises "`threshold` is None and thus
#                   scrublet requires `scikit-image`" at run time, not at
#                   import, so nothing catches it until a real sample runs.
#   leidenalg     - required by sc.tl.leiden(flavor='igraph')
#   hdf5plugin    - supplies the HDF5 filter plugins, Blosc above all.
#                   simpleaf writes quants.h5ad with the Blosc filter (id
#                   32001) once the matrix is over a small size threshold,
#                   and without this package h5py cannot read that file at
#                   all: "Can't synchronously read data (can't open
#                   directory ...)". It is needed by REPACK_H5AD, which
#                   rewrites the matrix as gzip so QCatch's own container
#                   can read it too. The threshold is why the gap went
#                   unnoticed: a toy fixture is written uncompressed and
#                   reads fine anywhere, while any real sample does not.
#   quarto        - renders report/analysis_report.qmd
#   ipykernel     - quarto's Python engine needs a kernel to execute the .qmd
#   procps-ng     - supplies `ps`, without which Nextflow cannot collect
#                   per-task CPU and memory metrics for the trace report
RUN micromamba create -y -p /opt/env \
        -c conda-forge -c bioconda \
        python=3.14.7=hcd007b5_106_cp314 \
        scanpy=1.12.4=pyhd8ed1ab_1 \
        anndata=0.13.3.post0=pyhd8ed1ab_0 \
        numpy=2.5.2=py314hd3a7d6b_1 \
        scipy=1.18.0=py314hf07bd8e_0 \
        pandas=3.0.5=py314hb4ffadd_1 \
        h5py=3.16.0=nompi_py314h77fd174_104 \
        hdf5plugin=7.1.0=py314haeb3a3a_0 \
        scikit-learn=1.9.0=np2py314hf09ca88_0 \
        scikit-misc=0.5.2=py314hacf7208_1 \
        scikit-image=0.26.0=np2py314hda1ea4c_0 \
        python-igraph=1.0.0=py314ha160325_0 \
        leidenalg=0.12.0=py314ha160325_0 \
        umap-learn=0.5.12=py314hdafbbf9_0 \
        matplotlib-base=3.11.1=py314h9a9c090_2 \
        seaborn=0.13.2=hd8ed1ab_3 \
        pyyaml=6.0.3=py314h67df5f8_1 \
        quarto=1.9.38=hd55b112_0 \
        jupyter_core=5.9.1=pyhc90fa1f_0 \
        nbclient=0.11.0=pyhcf101f3_1 \
        ipykernel=7.3.0=pyha191276_0 \
        procps-ng=4.0.6=h18c060e_0 \
    && micromamba clean --all --yes

# --- Stage 2: runtime ------------------------------------------------------
#
# Same base, same digest. Copying only /opt/env leaves the solver's package
# cache behind in stage 1, which is the whole point of the split. Using one
# base rather than a slimmer OS keeps a single image to audit and avoids
# guessing which system libraries the conda-forge builds and quarto's
# bundled deno runtime expect to find.
FROM ${MICROMAMBA_IMAGE} AS runtime

# The base image ends on USER mambauser, which cannot write /etc/passwd or
# /opt. Take root back for the setup below and drop privileges at the end.
USER root

COPY --from=build /opt/env /opt/env

# Put the prefix on PATH directly and clear the inherited ENTRYPOINT.
#
# The base image's entrypoint is an activation wrapper. Nextflow appends
# its own command to a container's entrypoint, so leaving an activating
# entrypoint in place makes every task depend on that wrapper behaving.
# Setting PATH removes the need for activation at all.
ENV PATH=/opt/env/bin:$PATH \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

ENTRYPOINT []

# Where HDF5 looks for filter plugins. Setting this is what makes h5py able
# to read simpleaf's Blosc-compressed quants.h5ad; the alternative,
# `import hdf5plugin` before h5py, would only work in code we control, and
# the variable also covers anything that opens the file indirectly.
# Verified: with this set, plain h5py reads a real simpleaf h5ad; without
# it, the same read fails with "can't open directory (.../hdf5/plugin)".
ENV HDF5_PLUGIN_PATH=/opt/env/lib/python3.14/site-packages/hdf5plugin/plugins

# nextflow.config runs this image with `-u $(id -u):$(id -g)`, so at run
# time the process has the HOST uid, which matches no entry in /etc/passwd
# and owns no home directory. Every library that wants to write a cache --
# matplotlib, numba, quarto's deno runtime -- must therefore be pointed at
# a world-writable location, or it warns on import and can fail outright.
ENV HOME=/tmp \
    XDG_CACHE_HOME=/tmp/.cache \
    XDG_CONFIG_HOME=/tmp/.config \
    XDG_DATA_HOME=/tmp/.local/share \
    MPLCONFIGDIR=/tmp/.config/matplotlib \
    NUMBA_CACHE_DIR=/tmp/.cache/numba

# fontconfig is the exception: it will not create its own cache directory,
# it only probes a fixed list and gives up. matplotlib pulls it in, so the
# directory has to exist up front or every plotting step prints
# "Fontconfig error: No writable cache directories".
RUN mkdir -p /tmp/.cache/numba /tmp/.cache/fontconfig /tmp/.config/matplotlib /tmp/.local/share \
    && chmod -R 1777 /tmp/.cache /tmp/.config /tmp/.local

# Replacement for conda activation, which we deliberately skip.
#
# Putting /opt/env/bin on PATH is enough for Python, but NOT for quarto.
# conda-forge ships /opt/env/etc/conda/activate.d/*.sh, which an activated
# environment would source and we never do. quarto's launcher then falls
# back to locating its bundled tools relative to its own directory: it
# sees the basename `bin`, concludes it is running from a source checkout,
# and looks under $QUARTO_BIN_PATH/tools/x86_64/. Nothing is there in a
# conda layout.
#
# The failures are staggered, which makes this easy to half-fix:
#   share path missing -> `quarto --version` dies on `cat: .../share/version`
#   QUARTO_DENO missing -> `.../tools/x86_64/deno: No such file or directory`
#   QUARTO_DART_SASS missing -> render dies in compileSass, AFTER the
#                               Python kernel has run, so the log looks
#                               like a successful execution
# Setting the whole set is the only way to be sure. Values mirror
# activate.d/quarto.sh, deno.sh and deno_dom.sh.
#
# QUARTO_DENO_DOM deviates on purpose: activate.d/quarto.sh ships a
# BROKEN value pointing at the feedstock build directory
# (/home/conda/feedstock_root/...), which does not exist in the image.
# deno_dom.sh has the real path, so use that.
ENV QUARTO_SHARE_PATH=/opt/env/share/quarto \
    QUARTO_CONDA_PREFIX=/opt/env \
    QUARTO_DENO=/opt/env/bin/deno \
    QUARTO_DENO_DOM=/opt/env/lib/deno_dom.so \
    QUARTO_PANDOC=/opt/env/bin/pandoc \
    QUARTO_ESBUILD=/opt/env/bin/esbuild \
    QUARTO_TYPST=/opt/env/bin/typst \
    QUARTO_DART_SASS=/opt/env/bin/sass \
    DENO_INSTALL_ROOT=/opt/env \
    DENO_DOM_PLUGIN=/opt/env/lib/deno_dom.so \
    GSETTINGS_SCHEMA_DIR=/opt/env/share/glib-2.0/schemas \
    TYPST_PACKAGE_PATH=/opt/env/share/typst/packages \
    TYPST_FONT_PATHS=/opt/env/fonts \
    XML_CATALOG_FILES="file:///opt/env/etc/xml/catalog file:///etc/xml/catalog"

# PYTHONNOUSERSITE is NOT set here. nextflow.config already sets it in its
# `env` scope, and a value belongs in exactly one file.

# Non-root by default. The `docker` profile overrides this with the host
# uid, but the `awsbatch` profile does not, so the default still matters.
# uid 1000 is free in this base; its own user is mambauser at 57439.
#
# /work is 1777 rather than owned by appuser because the `docker` profile
# runs as an arbitrary host uid that owns nothing in the image.
RUN useradd --create-home --shell /bin/bash --uid 1000 appuser \
    && mkdir -p /work \
    && chmod 1777 /work

USER appuser
WORKDIR /work

LABEL org.opencontainers.image.title="scrnaseq-lite" \
      org.opencontainers.image.description="Local-module runtime for the scrnaseq-lite pipeline: Scanpy, Quarto and the version/manifest scripts" \
      org.opencontainers.image.version="0.1.0" \
      org.opencontainers.image.source="https://github.com/aleimba/scrnaseq-lite" \
      org.opencontainers.image.licenses="MIT"
