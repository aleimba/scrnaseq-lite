#!/usr/bin/env bash
#
# download_reference.sh
#
# Downloads the 10x Genomics Human GRCh38 2024-A reference package, which
# supplies the genome FASTA and gene annotation GTF that the pipeline needs
# to build its simpleaf/piscem splici index.
#
# This script ONLY downloads. Index building is always done by the pipeline,
# so that the pinned container is the single source of truth for the simpleaf
# version, ALEVIN_FRY_HOME setup and open-file limits. Duplicating the index
# command here would create a second code path that could silently drift.
#
# Reference and layout follow the COMBINE-lab simpleaf/piscem tutorial:
#   https://combine-lab.github.io/alevin-fry-tutorials/2023/simpleaf-piscem/
#
# Usage:
#   bin/download_reference.sh [--outdir DIR] [--yes] [--force]
#
#   --outdir DIR   Where to place the reference (default: reference)
#   --yes, -y      Skip the confirmation prompt
#   --force        Re-download even if the reference is already present
#   --help, -h     Show this message
#
set -Eeuo pipefail

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

# 10x Genomics Human reference, GRCh38, version 2024-A.
# Unpacks to fasta/genome.fa and genes/genes.gtf.gz among other files.
# simpleaf accepts the gzipped GTF directly, so it is not decompressed here.
readonly REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
readonly REF_NAME="refdata-gex-GRCh38-2024-A"

OUTDIR="reference"
ASSUME_YES=0
FORCE=0

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

usage() {
    sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
}

require() {
    command -v "$1" >/dev/null 2>&1 || die "'$1' is required but not on PATH."
}

confirm() {
    local prompt="$1" reply
    if [[ "${ASSUME_YES}" -eq 1 ]]; then
        log "--yes given, proceeding: ${prompt}"
        return 0
    fi
    read -r -p "${prompt} [y/N] " reply
    [[ "${reply}" =~ ^[Yy]$ ]] || die "Aborted by user."
}

# Report the remote size without downloading, so the user can decide before
# committing disk space.
remote_size() {
    local url="$1" bytes
    bytes="$(curl -sIL "${url}" \
        | awk 'BEGIN{IGNORECASE=1} /^content-length:/ {v=$2} END{gsub(/\r/,"",v); print v}')"
    if [[ -n "${bytes}" && "${bytes}" =~ ^[0-9]+$ ]]; then
        awk -v b="${bytes}" 'BEGIN{printf "%.1f GB", b/1073741824}'
    else
        echo "unknown"
    fi
}

free_space_gb() {
    df -Pk "$1" | awk 'NR==2 {printf "%.1f", $4/1048576}'
}

# --------------------------------------------------------------------------
# Argument parsing
# --------------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)  OUTDIR="$2"; shift 2 ;;
        --yes|-y)  ASSUME_YES=1; shift ;;
        --force)   FORCE=1; shift ;;
        --help|-h) usage ;;
        *)         die "Unknown argument: $1 (try --help)" ;;
    esac
done

require curl
require tar
require awk

REF_DIR="${OUTDIR}/${REF_NAME}"
FASTA="${REF_DIR}/fasta/genome.fa"
GTF="${REF_DIR}/genes/genes.gtf.gz"

# --------------------------------------------------------------------------
# Download
# --------------------------------------------------------------------------

mkdir -p "${OUTDIR}"

log "Reference     : ${REF_URL}"
log "Target dir    : ${REF_DIR}"
log "Download size : $(remote_size "${REF_URL}")"
log "Free space on $(readlink -f "${OUTDIR}") : $(free_space_gb "${OUTDIR}") GB"
log ""
log "The archive is extracted as it streams, so peak usage is roughly the"
log "unpacked size rather than download plus unpacked. Allow headroom."

if [[ -f "${FASTA}" && -f "${GTF}" && "${FORCE}" -eq 0 ]]; then
    log "Reference already present, skipping download. Use --force to redo."
else
    confirm "Download the GRCh38 2024-A reference now?"
    mkdir -p "${REF_DIR}"
    log "Downloading and extracting..."
    # --strip-components=1 drops the archive's own top-level directory so the
    # contents land directly in ${REF_DIR}.
    curl -fL "${REF_URL}" | tar xzf - --strip-components=1 -C "${REF_DIR}"
    log "Done."
fi

[[ -f "${FASTA}" ]] || die "Expected FASTA not found: ${FASTA}"
[[ -f "${GTF}"   ]] || die "Expected GTF not found: ${GTF}"

log "FASTA: ${FASTA} ($(du -h "${FASTA}" | cut -f1))"
log "GTF  : ${GTF} ($(du -h "${GTF}" | cut -f1))"

# --------------------------------------------------------------------------
# Next step
# --------------------------------------------------------------------------

cat >&2 <<EOF

Reference ready. Build the splici index once, through the pipeline:

    nextflow run . -profile demo,docker \\
        --build_index --save_reference \\
        --fasta ${FASTA} \\
        --gtf   ${GTF} \\
        --r2_read_length 91 \\
        --run   build_index

The index is published under that run's results directory. Reuse it on
every later run:

    --simpleaf_index <run_dir>/simpleaf_index/index

Note the trailing 'index' directory rather than a file prefix: an index
built by 'simpleaf index' is referenced by its directory. A prefix such as
'index/piscem_idx' would apply only to an index built by calling
'piscem build' directly.

EOF
