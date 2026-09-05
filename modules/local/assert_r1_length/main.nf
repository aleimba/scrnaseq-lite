process ASSERT_R1_LENGTH {
    tag "${meta.id}"
    label 'process_single'

    // Local module, so it runs in the pipeline-authored image. Vendored
    // modules keep their own containers; see "Design section 2".
    container 'scrnaseq-lite:0.1.0'

    input:
    // in: [ meta, path(seqkit_stats_tsv), val(r1_filename), val(expected_r1_length) ]
    //
    // The R1 file NAME is passed rather than inferred. SEQKIT_STATS reports
    // one row per input file, keyed by the staged filename, and picking the
    // R1 row by row order or by a `_R1_` pattern would be a guess: nothing
    // guarantees either. The workflow already knows which file it passed.
    //
    // The expected length comes from `params.chemistry_map[...].r1_len`, so
    // this module reads no params of its own and the chemistry cannot drift
    // from the length it implies.
    tuple val(meta), path(seqkit_stats_tsv), val(r1_filename), val(expected_r1_length)

    output:
    // out: [ meta, path(tsv) ] - one row, the observed lengths and the verdict
    tuple val(meta), path("*.r1_length_check.tsv"), emit: check
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //"'), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A bash script block calling bin/assert_r1_length.py. Bash is required,
    // not stylistic: Nextflow rejects an `eval` output on any other
    // interpreter ("Process output of type 'eval' is only allowed with Bash
    // process scripts"), and the version above is reported through the
    // versions topic, which is an eval.
    """
    assert_r1_length.py \\
        --stats '${seqkit_stats_tsv}' \\
        --r1-filename '${r1_filename}' \\
        --expected-length ${expected_r1_length} \\
        --sample '${meta.id}' \\
        --chemistry '${meta.chemistry}' \\
        --output '${prefix}.r1_length_check.tsv'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.r1_length_check.tsv
    """
}
