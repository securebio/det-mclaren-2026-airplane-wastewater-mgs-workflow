// Separate ribosomal reads with SortMeRNA against a prebuilt index
// (bin/build_sortmerna_index.sh). For interleaved input, a pair is ribosomal if either
// mate aligns (--paired_in), and both outputs stay interleaved.
process SORTMERNA {
    label "small"
    label "sortmerna"
    tag "id=${sample}"
    input:
        tuple val(sample), path(reads) // Interleaved or single-end
        path(ref_fasta) // Reference FASTA the index was built from
        path(index_dir) // Prebuilt index directory
        val(params_map) // suffix, interleaved
    output:
        tuple val(sample), path("${sample}_${params_map.suffix}_sortmerna_nomatch.fastq.gz"), emit: nomatch
        tuple val(sample), path("${sample}_${params_map.suffix}_sortmerna_match.fastq.gz"), emit: match
        tuple val(sample), path("${sample}_${params_map.suffix}_sortmerna.log"), emit: log
        tuple val(sample), path("input_${reads}"), emit: input
    script:
        def match = "${sample}_${params_map.suffix}_sortmerna_match.fastq.gz"
        def nomatch = "${sample}_${params_map.suffix}_sortmerna_nomatch.fastq.gz"
        def log = "${sample}_${params_map.suffix}_sortmerna.log"
        def paired = params_map.interleaved ? "--paired_in --out2" : ""
        """
        set -euo pipefail
        # SortMeRNA names index files by a hash of the reference FASTA's file name, and the
        # index is built from a FASTA named ref.fasta (bin/build_sortmerna_index.sh), so
        # link the reference under that name
        mkdir smr_ref
        ln -s "\${PWD}/${ref_fasta}" smr_ref/ref.fasta
        # SortMeRNA doesn't check that the index matches the FASTA, so compare their
        # sequence counts
        idx_numseq=\$(awk '\$1 == "numseq:" {print \$2; exit}' ${index_dir}/*.idx_stats.yaml)
        ref_numseq=\$(grep -c '^>' smr_ref/ref.fasta)
        if [[ "\${idx_numseq}" != "\${ref_numseq}" ]]; then
            echo "SortMeRNA index has \${idx_numseq:-no} sequences but ${ref_fasta} has \${ref_numseq}" >&2
            exit 1
        fi
        # Keep SortMeRNA's key-value store and read buffers on local disk
        smr_dir=\$(mktemp -d "\${TMPDIR:-/tmp}/sortmerna.XXXXXX")
        # SortMeRNA crashes on empty input. Read the first byte to test for it, tolerating
        # only the broken pipe that closing the stream early causes.
        gzip -dcf ${reads} | head -c 1 > "\${smr_dir}/first_byte" || [[ \${PIPESTATUS[0]} -eq 141 ]]
        if [[ ! -s "\${smr_dir}/first_byte" ]]; then
            printf "" | gzip > ${match}
            printf "" | gzip > ${nomatch}
            echo "Empty input: SortMeRNA not run" > ${log}
        else
            sortmerna --ref smr_ref/ref.fasta --idx-dir ${index_dir} --index 0 \\
                --reads ${reads} ${paired} --fastx --zip-out 1 --num_alignments 1 \\
                --aligned "\${smr_dir}/match" --other "\${smr_dir}/nomatch" \\
                --workdir "\${smr_dir}/work" --threads ${task.cpus}
            # Check the outputs are intact before reformatting them; failures inside the
            # process substitutions below would not stop the task
            gzip -t "\${smr_dir}"/*.fq.gz
            if [[ "${params_map.interleaved}" == "true" ]]; then
                # Re-interleave the mates, which --out2 writes to separate files in step
                for out in match nomatch; do
                    paste <(gzip -dc "\${smr_dir}/\${out}_fwd.fq.gz" | paste - - - -) \\
                          <(gzip -dc "\${smr_dir}/\${out}_rev.fq.gz" | paste - - - -) \\
                        | tr '\\t' '\\n' | pigz -p ${task.cpus} > "${sample}_${params_map.suffix}_sortmerna_\${out}.fastq.gz"
                done
            else
                mv "\${smr_dir}/match.fq.gz" ${match}
                mv "\${smr_dir}/nomatch.fq.gz" ${nomatch}
            fi
            mv "\${smr_dir}/match.log" ${log}
        fi
        rm -rf "\${smr_dir}"
        ln -s ${reads} input_${reads}
        """
}
