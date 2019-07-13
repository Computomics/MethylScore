class ParameterChecks {
  static void checkParams(params) {
    assert params.SAMPLE_SHEET, "samplesheet.tsv has to be specified!"
    assert params.GENOME, "reference genome in fasta format has to be specified!"

    assert params.HUMAN instanceof Boolean, "HUMAN must be set to either false (off) or true (on)"
    assert params.IGV instanceof Boolean, "IGV must be set to either false (off) or true (on)"
    assert params.STATISTICS instanceof Boolean, "STATISTICS must be set to either false (off) or true (on)"
    assert params.DO_DEDUP instanceof Boolean, "DO_DEDUP must be set to either false (off) or true (on)"
    assert params.DMRS_PER_CONTEXT instanceof Boolean, "DMRS_PER_CONTEXT must be set to either false (off) or true (on)"

    assert params.DMR_CONTEXTS.tokenize(',').size() in 1..3 && params.DMR_CONTEXTS =~ "C[G,H]{1,2}", "DMR_CONTEXTS must be CG, CHG or CHH!"
    assert params.MR_FREQ_CHANGE in 0..100, "MR_FREQ_CHANGE must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH_DIFF in 0..100, "CLUSTER_MIN_METH_DIFF must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH_DIFF_CG in 0..100 || DMRS_PER_CONTEXT == false, "CLUSTER_MIN_METH_DIFF_CG must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH_DIFF_CHG in 0..100 || DMRS_PER_CONTEXT == false, "CLUSTER_MIN_METH_DIFF_CHG must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH_DIFF_CHH in 0..100 || DMRS_PER_CONTEXT == false, "CLUSTER_MIN_METH_DIFF_CHH must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH_CG in 0..100 || DMRS_PER_CONTEXT == false, "CLUSTER_MIN_METH_CG must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH_CHG in 0..100 || DMRS_PER_CONTEXT == false, "CLUSTER_MIN_METH_CHG must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH_CHH in 0..100 || DMRS_PER_CONTEXT == false, "CLUSTER_MIN_METH_CHH must be between 0 and 100!"
    assert params.CLUSTER_MIN_METH in 0..100, "CLUSTER_MIN_METH must be between 0 and 100!"
    assert params.MR_FREQ_DISTANCE instanceof Integer && params.MR_FREQ_DISTANCE >= 0, "MR_FREQ_DISTANCE must be a non-negative integer!"
    assert params.SLIDING_WINDOW_SIZE instanceof Integer && params.SLIDING_WINDOW_SIZE >= 0, "SLIDING_WINDOW_SIZE must be a non-negative integer!"
    assert params.SLIDING_WINDOW_STEP instanceof Integer && params.SLIDING_WINDOW_STEP >= 0, "SLIDING_WINDOW_STEP must be a non-negative integer!"
    assert params.DMR_MIN_C instanceof Integer && params.DMR_MIN_C >= 0, "DMR_MIN_C must be a non-negative integer!"
    assert params.DMR_MIN_COV instanceof Integer && params.DMR_MIN_COV >= 0, "DMR_MIN_COV must be a non-negative integer!"
    assert params.MR_BATCH_SIZE instanceof Integer && params.MR_BATCH_SIZE >= 0, "MR_BATCH_SIZE must be a non-negative integer!"
    assert params.HDMR_FOLD_CHANGE >= 0, "HDMR_FOLD_CHANGE must be a non-negative number!"
    assert params.FDR_CUTOFF > 0 && params.FDR_CUTOFF < 1, "FDR_CUTOFF must be between 0 and 1!"

    assert params.MIN_QUAL instanceof Integer && params.MIN_QUAL in 1..40, "MIN_QUAL must be between 1 and 40!"
    assert params.IGNORE_FIRST_BP instanceof Integer && params.IGNORE_FIRST_BP >= 0, "IGNORE_FIRST_BP must be a non-negative integer!"
    assert params.IGNORE_LAST_BP instanceof Integer && params.IGNORE_LAST_BP >= 0, "IGNORE_LAST_BP must be a non-negative integer!"

    assert params.MIN_COVERAGE instanceof Integer && params.MIN_COVERAGE > 0, "MIN_COVERAGE must be a non-negative integer!"
    assert params.DESERT_SIZE instanceof Integer && params.DESERT_SIZE > 0, "DESERT_SIZE must be a non-negative integer!"
    assert params.MERGE_DIST instanceof Integer && params.MERGE_DIST > 0, "MERGE_DIST must be a non-negative integer!"
    assert params.MR_MIN_C instanceof Integer, "MR_MIN_C must be an integer! Negative value turns on permutation test."
    assert params.TRIM_METHRATE in 0..100, "TRIM_METHRATE must be between 0 and 100!" 
  }
}
