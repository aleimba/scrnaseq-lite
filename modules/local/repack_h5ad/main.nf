process REPACK_H5AD {
    tag "${meta.id}"
    label 'process_single'

    // Local module, so it runs in the pipeline-authored image. That image is
    // the only one here carrying the HDF5 Blosc filter (via hdf5plugin and
    // HDF5_PLUGIN_PATH), which is precisely what makes this step possible.
    container 'scrnaseq-lite:0.1.0'

    input:
    // in: [ meta, path(quants_h5ad) ]  - af_quant/alevin/quants.h5ad from SIMPLEAF_QUANT
    tuple val(meta), path(quants_h5ad)

    output:
    // out: [ meta, path(h5ad) ] - same data, gzip-compressed, readable anywhere.
    // The name carries the provenance suffix R5.5 requires, so no later step
    // has to rename it.
    tuple val(meta), path("*_raw_matrix.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('h5py'), eval('python -c "import h5py; print(h5py.__version__)"'), topic: versions, emit: versions_h5py

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // simpleaf writes this file with the Blosc HDF5 filter (id 32001) on any
    // real sample. QCatch's container has no Blosc plugin and fails while
    // parsing its arguments, before any analysis; gzip is built into every
    // HDF5 build, so a gzip copy is readable everywhere. See R6.4b.
    """
    repack_h5ad_blosc_to_gzip.py \\
        --input '${quants_h5ad}' \\
        --output '${prefix}_raw_matrix.h5ad'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_raw_matrix.h5ad
    """
}
