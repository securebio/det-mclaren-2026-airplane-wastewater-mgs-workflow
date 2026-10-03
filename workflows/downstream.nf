/***********************************************************
| WORKFLOW: DOWNSTREAM ANALYSIS OF PRIMARY WORKFLOW OUTPUT |
***********************************************************/

/***************************
| MODULES AND SUBWORKFLOWS |
***************************/

include { LOAD_DOWNSTREAM_DATA } from "../subworkflows/local/loadDownstreamData"
include { DISCOVER_RUN_OUTPUT } from "../subworkflows/local/discoverRunOutput"
include { CONCAT_RUN_OUTPUTS_BY_GROUP } from "../subworkflows/local/concatRunOutputsByGroup"
include { MARK_VIRAL_DUPLICATES } from "../subworkflows/local/markViralDuplicates"
include { COUNT_READS_PER_CLADE } from "../modules/local/countReadsPerClade"
include { COPY_FILE_BARE as COPY_PYPROJECT } from "../modules/local/copyFile"
include { COPY_FILE_BARE as COPY_INPUT } from "../modules/local/copyFile"
include { WRITE_SENTINEL_DOWNSTREAM } from "../modules/local/writeSentinelDownstream"

/*****************
| MAIN WORKFLOWS |
*****************/

workflow DOWNSTREAM {
    main:
        // Prepare channels from input CSV file
        load_ch = LOAD_DOWNSTREAM_DATA(params.input_file, params.input_base_dir ?: projectDir)
        start_time_str = load_ch.start_time_str
        // Discover all per-sample output files and match to groups
        pipeline_pyproject_path = file("${projectDir}/pyproject.toml")
        discover_ch = DISCOVER_RUN_OUTPUT(load_ch.run_dirs, load_ch.groups, pipeline_pyproject_path, params.platform).output
        // Concatenate per-sample outputs into per-group TSVs
        concat_ch = CONCAT_RUN_OUTPUTS_BY_GROUP(discover_ch)
        // Prepare inputs for clade counting
        viral_db_path = "${params.ref_dir}/results/total-virus-db-annotated.tsv.gz"
        viral_db = channel.value(viral_db_path)
        // Short-read: mark duplicates by sequence similarity among all reads, then count
        // reads per clade. ONT skips both, so its DOWNSTREAM output has no viral hits.
        if (params.platform == "ont") {
            clade_counts_ch = channel.empty()
        }
        else {
            mark_dup_ch = MARK_VIRAL_DUPLICATES(concat_ch.hits, params.sim_dup_max_error_frac)
            clade_counts_ch = COUNT_READS_PER_CLADE(mark_dup_ch.hits, viral_db).output
        }
        // Viral assignments are not validated with BLAST, so no validation_hits are produced.
        // The duplicate-marked hits stay in the MARK_SIMILARITY_DUPLICATES and
        // COUNT_READS_PER_CLADE work directories.
        // Prepare publishing channels
        params_str = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
        params_ch = channel.of(params_str).collectFile(name: "params-downstream.json")
        pyproject_ch = COPY_PYPROJECT(channel.fromPath(pipeline_pyproject_path), "pyproject.toml")
        input_file_ch = COPY_INPUT(channel.fromPath(params.input_file), "input_file.csv")

        // Pre-define publish-channel aggregates so we can both emit them and feed them into the sentinel barrier
        input_downstream_ch = params_ch.mix(input_file_ch)
        logging_downstream_ch = pyproject_ch
        results_downstream_ch = clade_counts_ch.mix(
                                    concat_ch.other,
                                    concat_ch.fastp_json)

        // Validate published outputs and write per-group sentinels
        groups_only_ch = load_ch.groups
            .map { _label, _sample, group -> group }
            .unique()
        sentinel_params = params + [output_dir: "${params.base_dir}/output", pyproject_path: "${projectDir}/pyproject.toml"]
        sentinel_ch = WRITE_SENTINEL_DOWNSTREAM(
            groups_only_ch,
            input_downstream_ch.mix(logging_downstream_ch, results_downstream_ch).collect(),
            start_time_str,
            sentinel_params
        )

    emit:
        input_downstream = input_downstream_ch
        logging_downstream = logging_downstream_ch
        intermediates_downstream = channel.empty()
        results_downstream = results_downstream_ch
        experimental_downstream = channel.empty()
        sentinel_downstream = sentinel_ch.sentinel
}
