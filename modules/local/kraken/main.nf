// Perform taxonomic assignment with Kraken2 on streamed data
process KRAKEN {
    label "Kraken2"
    label "kraken_resources"
    tag "id=${sample}"
    input:
        tuple val(sample), path(reads)
        val db_path
        val db_download_timeout
        val label // Distinguishes per-read output files published from one sample, e.g. "ribo"
    output:
        tuple val(sample), path("${sample}_kraken_output_${label}.tsv.gz"), emit: output
        tuple val(sample), path("${sample}.report.gz"), emit: report
        tuple val(sample), path("input_${reads}"), emit: input
    script:
        def extractCmd = reads.toString().endsWith(".gz") ? "zcat" : "cat"
        def out = "${sample}_kraken_output_${label}.tsv"
        def report = "${sample}.report"
        def par = "--use-names --report-minimizer-data --threads ${task.cpus} --report ${report} --memory-mapping"
        """
        # Download Kraken2 database if not already present
        db_local_path=\$(download_db.py "${db_path}" ${db_download_timeout})
        # Run Kraken
        ${extractCmd} ${reads} | kraken2 --db \${db_local_path} ${par} /dev/fd/0 > ${out}
        # Make empty output files if needed
        touch ${out}
        touch ${report}
        # Gzip output and report to save space
        gzip ${out}
        gzip ${report}
        # Link input to output for testing
        ln -s ${reads} input_${reads}
        """
}
