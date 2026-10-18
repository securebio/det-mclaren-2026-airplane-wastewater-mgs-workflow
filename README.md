# Nucleic Acid Observatory Viral Metagenomics Pipeline

> [!IMPORTANT]
> This repository is a fork of [nao-mgs-workflow](https://github.com/securebio/nao-mgs-workflow) 3.4.0.1, modified for the analysis in [det-mclaren-2026-airplane-wastewater](https://github.com/securebio/det-mclaren-2026-airplane-wastewater). Its versions carry the local label `+mclaren2026.N` (see `pyproject.toml`), and the top entry of [CHANGELOG.md](CHANGELOG.md) lists how it differs from upstream. The SortMeRNA reference and index used by RUN aren't included: build them with `bin/build_sortmerna_index.sh`.

This Nextflow pipeline is designed to process metagenomic sequencing data, characterize overall taxonomic composition, and identify and quantify reads mapping to viruses infecting certain host taxa of interest. It was developed as part of the [Nucleic Acid Observatory](https://naobservatory.org/) project.

The pipeline currently consists of three workflows:

- [`INDEX`](./docs/index.md): Creates indices and reference files used by the `RUN` and `DOWNSTREAM` workflows[^1].
- [`RUN`](./docs/run.md): Performs the main analysis, including QC, viral identification and taxonomic profiling.
- [`DOWNSTREAM`](./docs/downstream.md): Performs downstream analysis of the results from the `RUN` workflow[^2].

[^1]: The `INDEX` workflow is intended to be run first, after which many instantiations of the `RUN` workflow can use the same index output files.
[^2]: The `DOWNSTREAM` workflow is designed to handle tasks that require cross-read comparisons, including potentially across multiple runs: marking duplicate reads and counting reads per viral clade.

## Documentation

- **Installation and usage:**
    - [Installation instructions](docs/installation.md)
    - [AWS Batch setup](docs/batch.md)
    - [Usage instructions](docs/usage.md)
    - [Troubleshooting](docs/troubleshooting.md)
- **Workflow details:**
    - [INDEX workflow](docs/index.md)
    - [RUN workflow](docs/run.md)
    - [DOWNSTREAM workflow](docs/downstream.md)
- **Configuration and output:**
    - [Configuration files](docs/config.md)
    - [Pipeline outputs](docs/output.md)
- **Other:**
    - [Versioning](docs/versioning.md)
    - [Developer docs](docs/developer.md)
    - [Viral infection status annotation](docs/annotation.md)
    - [Reproducibility](docs/reproducibility.md)
