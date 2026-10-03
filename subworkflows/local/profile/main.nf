/**********************************************
| SUBWORKFLOW: HIGH-LEVEL TAXONOMIC PROFILING |
**********************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { SORTMERNA } from "../../../modules/local/sortmerna"
include { MINIMAP2 } from "../../../modules/local/minimap2"
include { TAXONOMY as TAXONOMY_RIBO } from "../../../subworkflows/local/taxonomy"
include { TAXONOMY as TAXONOMY_NORIBO } from "../../../subworkflows/local/taxonomy"
include { ADD_FIXED_COLUMN as ADD_KRAKEN_RIBO } from "../../../modules/local/addFixedColumn"
include { ADD_FIXED_COLUMN as ADD_BRACKEN_RIBO } from "../../../modules/local/addFixedColumn"
include { ADD_FIXED_COLUMN as ADD_KRAKEN_NORIBO } from "../../../modules/local/addFixedColumn"
include { ADD_FIXED_COLUMN as ADD_BRACKEN_NORIBO } from "../../../modules/local/addFixedColumn"
include { CONCATENATE_TSVS_LABELED as CONCATENATE_KRAKEN_PER_SAMPLE } from "../../../modules/local/concatenateTsvs"
include { CONCATENATE_TSVS_LABELED as CONCATENATE_BRACKEN_PER_SAMPLE } from "../../../modules/local/concatenateTsvs"

/****************
| MAIN WORKFLOW |
****************/

workflow PROFILE {
    take:
        reads_ch
        single_end
        params_map // Uses: ribo_suffix, bracken_threshold, platform, db_download_timeout, ref_dir, sortmerna_ref, sortmerna_index_dir
    main:
        kraken_db_ch = "${params_map.ref_dir}/results/kraken_db"
        // Separate ribosomal reads
        if (params_map.platform == "ont") {
            ribo_ref = "${params_map.ref_dir}/results/mm2-ribo-index"
            ribo_minimap2_params = [
                suffix: params_map.ribo_suffix,
                remove_sq: false,
                alignment_params: "",
                db_download_timeout: params_map.db_download_timeout
            ]
            ribo_ch = MINIMAP2(reads_ch, ribo_ref, ribo_minimap2_params)
            ribo_in = ribo_ch.reads_mapped
            noribo_in = ribo_ch.reads_unmapped
        } else {
            // SortMeRNA against a prebuilt index of its standard sensitive database
            // (bin/build_sortmerna_index.sh); a pair is ribosomal if either mate aligns
            if (!params_map.sortmerna_ref || !params_map.sortmerna_index_dir) {
                error "Set params.sortmerna_ref and params.sortmerna_index_dir to a prebuilt SortMeRNA reference and index (see bin/build_sortmerna_index.sh)."
            }
            ribo_smr_params = single_end.map { v -> [
                suffix: params_map.ribo_suffix,
                interleaved: !v
            ] }
            ribo_ch = SORTMERNA(reads_ch, params_map.sortmerna_ref, params_map.sortmerna_index_dir, ribo_smr_params)
            ribo_in = ribo_ch.match
            noribo_in = ribo_ch.nomatch
        }
        // Run taxonomic profiling separately on ribo and non-ribo reads
        taxonomy_params = params_map + [classification_level: "D"]
        tax_ribo_ch = TAXONOMY_RIBO(ribo_in, kraken_db_ch, single_end, taxonomy_params)
        tax_noribo_ch = TAXONOMY_NORIBO(noribo_in, kraken_db_ch, single_end, taxonomy_params)
        // Add ribosomal status to output TSVs
        kr_ribo = ADD_KRAKEN_RIBO(tax_ribo_ch.kraken_reports, "ribosomal", "TRUE", "ribo")
        kr_noribo = ADD_KRAKEN_NORIBO(tax_noribo_ch.kraken_reports, "ribosomal", "FALSE", "noribo")
        br_ribo = ADD_BRACKEN_RIBO(tax_ribo_ch.bracken, "ribosomal", "TRUE", "ribo")
        br_noribo = ADD_BRACKEN_NORIBO(tax_noribo_ch.bracken, "ribosomal", "FALSE", "noribo")
        // Concatenate ribo + noribo for each sample
        kr_combined = kr_ribo.output.join(kr_noribo.output)
            .map { sample, ribo_file, noribo_file -> [sample, [ribo_file, noribo_file]] }
        kr_per_sample = CONCATENATE_KRAKEN_PER_SAMPLE(kr_combined, "kraken")
        br_combined = br_ribo.output.join(br_noribo.output)
            .map { sample, ribo_file, noribo_file -> [sample, [ribo_file, noribo_file]] }
        br_per_sample = CONCATENATE_BRACKEN_PER_SAMPLE(br_combined, "bracken")
    emit:
        bracken = br_per_sample.output
        kraken = kr_per_sample.output
        // Extra outputs for testing
        test_ribo = ribo_in
        test_noribo = noribo_in
}
