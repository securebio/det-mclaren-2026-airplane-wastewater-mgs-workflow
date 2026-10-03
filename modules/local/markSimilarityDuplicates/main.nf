// Tool source: rust-tools/mark_duplicates_similarity/
process MARK_SIMILARITY_DUPLICATES {
    label "mark_similarity_duplicates_resources"
    label "rust_tools"
    tag "id=${sample}"
    input:
        tuple val(sample), path(tsv)
        val(max_error_frac) // Mismatches allowed per mate, as a fraction of the overlap
    output:
        tuple val(sample), path("similarity-duplicate-marked_${tsv}"), emit: output
        tuple val(sample), path("input_${tsv}"), emit: input
    script:
    """
    mark_duplicates_similarity -i "${tsv}" -o "similarity-duplicate-marked_${tsv}" \\
        --max-error-frac ${max_error_frac}
    ln -s "${tsv}" "input_${tsv}"
    """
}
