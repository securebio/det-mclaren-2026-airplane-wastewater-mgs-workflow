/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORT_TSV as SORT_READS } from "../../../modules/local/sortTsv"
include { MARK_SIMILARITY_DUPLICATES } from "../../../modules/local/markSimilarityDuplicates"

/***********
| WORKFLOW |
***********/

workflow MARK_VIRAL_DUPLICATES {
    take:
        groups // Labeled viral hit TSVs partitioned by group
        max_error_frac // Mismatches allowed per mate in similarity marking, as a fraction of the overlap
    main:
        // 1. Sort by read ID, so similarity marking sees reads in a deterministic order.
        reads_sorted_ch = SORT_READS(groups, "seq_id").sorted
        // 2. Mark similarity duplicates among all reads. Alignment-coordinate marking is
        // skipped: it merges distinct molecules at the same coordinates, and its grouping
        // is quadratic in the depth of identical-coordinate stacks. Without a
        // prim_align_dup_exemplar column in its input, the tool treats every read as a
        // candidate and adds the column with each read as its own exemplar.
        sim_dup_ch = MARK_SIMILARITY_DUPLICATES(reads_sorted_ch, max_error_frac).output
    emit:
        // Hits annotated by similarity duplicate marking
        hits = sim_dup_ch
        // Extra outputs for testing
        test_in = groups
        test_sorted = reads_sorted_ch
}
