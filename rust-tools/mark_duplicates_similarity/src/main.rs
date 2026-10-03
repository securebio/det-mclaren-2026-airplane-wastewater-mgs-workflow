use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use nao_dedup::{DedupContext, DedupParams, MinimizerParams, ReadPair};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

// ------------------------------------------------------------------------------------------------
// ARGUMENT PARSING
// ------------------------------------------------------------------------------------------------

/// Mark similarity-based duplicates among candidate reads.
///
/// If the input has a prim_align_dup_exemplar column, only alignment-unique reads (their own
/// alignment exemplar) are candidates. If it has none (alignment marking was skipped), every
/// read is a candidate, read IDs must be unique, and the column is added to the output with
/// each read as its own exemplar.
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input gzipped TSV file path
    #[arg(short, long)]
    input: String,
    /// Output gzipped TSV file path
    #[arg(short, long)]
    output: String,
    /// Mismatches allowed per mate, as a fraction of the mates' overlap (rounded down);
    /// an offset between mates spends this budget too, and also shrinks the overlap
    #[arg(long, default_value_t = DedupParams::default().max_error_frac)]
    max_error_frac: f64,
}

// ------------------------------------------------------------------------------------------------
// HELPER FUNCTIONS
// ------------------------------------------------------------------------------------------------

fn find_column(header_fields: &[&str], name: &str) -> Result<usize> {
    header_fields
        .iter()
        .position(|&f| f == name)
        .with_context(|| format!("Missing required column: {}", name))
}

// ------------------------------------------------------------------------------------------------
// MAIN
// ------------------------------------------------------------------------------------------------

fn main() -> Result<()> {
    let args = Args::parse();
    let input_path = &args.input;
    let output_path = &args.output;

    if !(0.0..1.0).contains(&args.max_error_frac) {
        bail!(
            "--max-error-frac must be at least 0 and less than 1, got {}",
            args.max_error_frac
        );
    }

    let start_time = Instant::now();

    // Create deduplication context with default parameters, apart from the error budget
    let dedup_params = DedupParams {
        max_error_frac: args.max_error_frac,
        ..DedupParams::default()
    };
    let minimizer_params = MinimizerParams::default();
    let mut ctx = DedupContext::new(dedup_params, minimizer_params);

    // Pass 1: Process candidate reads
    let mut n_reads = 0;
    let mut alignment_unique_count = 0;
    let mut align_dup_counts: HashMap<String, u32> = HashMap::new();

    eprintln!("Running similarity-based deduplication on candidate reads...");

    let file = File::open(input_path)
        .with_context(|| format!("Cannot open input file: {}", input_path))?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let mut lines = reader.lines();

    // Read header
    let header = lines
        .next()
        .context("Empty input file")?
        .context("Failed to read header")?;

    let header_fields: Vec<&str> = header.split('\t').collect();

    // Find column indices
    let seq_id_idx = find_column(&header_fields, "seq_id")?;
    let query_seq_idx = find_column(&header_fields, "query_seq")?;
    let query_seq_rev_idx = find_column(&header_fields, "query_seq_rev")?;
    let query_qual_idx = find_column(&header_fields, "query_qual")?;
    let query_qual_rev_idx = find_column(&header_fields, "query_qual_rev")?;
    // Absent when alignment marking was skipped: every read is then alignment-unique
    let prim_align_idx = header_fields
        .iter()
        .position(|&f| f == "prim_align_dup_exemplar");

    // Minimum number of fields required to access all required columns
    let min_fields = [
        seq_id_idx,
        query_seq_idx,
        query_seq_rev_idx,
        query_qual_idx,
        query_qual_rev_idx,
        prim_align_idx.unwrap_or(0),
    ]
    .into_iter()
    .max()
    .unwrap()
        + 1;

    // Process reads
    for line_result in lines {
        let line = line_result.context("Failed to read line")?;
        n_reads += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < min_fields {
            bail!(
                "Malformed line {}: expected at least {} fields, got {}",
                n_reads,
                min_fields,
                fields.len()
            );
        }

        let seq_id = fields[seq_id_idx];
        let prim_align_exemplar = prim_align_idx.map_or(seq_id, |i| fields[i]);

        // Count all reads (including alignment dups) per alignment-unique exemplar
        *align_dup_counts
            .entry(prim_align_exemplar.to_string())
            .or_insert(0) += 1;

        // Only process alignment-unique reads
        if seq_id != prim_align_exemplar {
            continue;
        }

        alignment_unique_count += 1;

        let read_pair = ReadPair {
            read_id: seq_id.to_string(),
            fwd_seq: fields[query_seq_idx].to_string(),
            rev_seq: fields[query_seq_rev_idx].to_string(),
            fwd_qual: fields[query_qual_idx].to_string(),
            rev_qual: fields[query_qual_rev_idx].to_string(),
        };

        ctx.process_read(read_pair);
    }

    // Without alignment marking, each read keys its own entry, so a count above 1 means a
    // repeated read ID, which would merge distinct reads in the similarity context
    if prim_align_idx.is_none() {
        if let Some((read_id, _)) = align_dup_counts.iter().find(|(_, &count)| count > 1) {
            bail!("Duplicate read ID in input: {}", read_id);
        }
    }

    let (_total_processed, unique_clusters) = ctx.stats();
    eprintln!(
        "Processed {} candidate reads (out of {} total reads)",
        alignment_unique_count, n_reads
    );
    eprintln!("Found {} unique sequence clusters", unique_clusters);

    // Finalize Pass 1
    ctx.finalize();

    // Compute similarity duplicate group sizes.
    // For each alignment-unique exemplar, map it to its similarity exemplar,
    // then aggregate the alignment dup counts by similarity exemplar.
    let mut group_sizes: HashMap<String, u32> = HashMap::new();
    for (au_id, count) in &align_dup_counts {
        let sim_exemplar = ctx.get_cluster_id(au_id);
        *group_sizes.entry(sim_exemplar).or_insert(0) += count;
    }
    drop(align_dup_counts);

    // Pass 2: Write output with sim_dup_exemplar column
    eprintln!("Pass 2: Writing output with sim_dup_exemplar column...");

    let file_in = File::open(input_path)
        .with_context(|| format!("Cannot open input file: {}", input_path))?;
    let decoder = GzDecoder::new(file_in);
    let reader = BufReader::new(decoder);
    let mut lines = reader.lines();

    let file_out = File::create(output_path)
        .with_context(|| format!("Cannot create output file: {}", output_path))?;
    let encoder = GzEncoder::new(file_out, Compression::default());
    let mut writer = BufWriter::new(encoder);

    // Skip header line and write stored header with new columns
    lines.next();
    let added_prim_align_col = if prim_align_idx.is_none() {
        "\tprim_align_dup_exemplar"
    } else {
        ""
    };
    writeln!(
        writer,
        "{}{}\tsim_dup_exemplar\tsim_dup_group_size",
        header.trim_end(),
        added_prim_align_col
    )
    .context("Failed to write header")?;

    let mut n_prim_align_dups = 0;
    let mut n_sim_dups = 0;

    // Process data rows
    for (line_num, line_result) in lines.enumerate() {
        let line = line_result.context("Failed to read line")?;
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < min_fields {
            bail!(
                "Malformed line {}: expected at least {} fields, got {}",
                line_num + 1,
                min_fields,
                fields.len()
            );
        }

        let seq_id = fields[seq_id_idx];
        let prim_align_exemplar = prim_align_idx.map_or(seq_id, |i| fields[i]);
        // Without alignment marking, each read is its own alignment exemplar
        let added_prim_align_val = if prim_align_idx.is_none() {
            format!("\t{}", seq_id)
        } else {
            String::new()
        };

        if seq_id != prim_align_exemplar {
            // Alignment duplicate - fast path
            writeln!(writer, "{}\tNA\tNA", line.trim_end()).context("Failed to write line")?;
            n_prim_align_dups += 1;
        } else {
            // Alignment-unique - query for similarity exemplar
            let sim_exemplar = ctx.get_cluster_id(seq_id);

            // Only the surviving similarity exemplar gets the nonzero group size;
            // reads that are sim dups (pointing to a different exemplar) get NA.
            let group_size = if sim_exemplar == seq_id {
                let size = group_sizes.get(&sim_exemplar).with_context(|| {
                    format!("sim_exemplar {} missing from group_sizes", sim_exemplar)
                })?;
                format!("{}", size)
            } else {
                "NA".to_string()
            };

            writeln!(
                writer,
                "{}{}\t{}\t{}",
                line.trim_end(),
                added_prim_align_val,
                sim_exemplar,
                group_size
            )
            .context("Failed to write line")?;

            if sim_exemplar != seq_id {
                n_sim_dups += 1;
            }
        }
    }

    writer.flush().context("Failed to flush output")?;

    let elapsed = start_time.elapsed().as_secs();
    eprintln!("Done!");

    let candidates = n_reads - n_prim_align_dups;
    let final_exemplars = candidates - n_sim_dups;
    eprintln!(
        "Processed {n_reads} reads in {elapsed}s. \
         {candidates} were similarity candidates (not alignment duplicates), \
         reduced to {final_exemplars} exemplars via similarity marking."
    );

    Ok(())
}
