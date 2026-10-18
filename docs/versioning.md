# Pipeline versioning

From version 2.6.0.0 we're adopting a new 4-number versioning scheme, described below. The primary purpose of this system is to communicate clearly to users of the pipeline what changes they must make when interpreting or using the outputs of the pipeline in downstream applications

1. **Major:** The first number in the version will be incremented each time we seriously rework the pipeline, requiring potentially major changes to downstream code.
2. **Schema:** The second number will be incremented each time we restructure or rename existing pipeline outputs, requiring downstream code to be changed to correctly access output files. This includes renaming or removing files from the `expected-outputs-*` lists in `pyproject.toml`, and any changes to `schemas/*.schema.json` files beyond `title` and `description` fields (which are documentation-only and permitted in point releases). Adding new output files does **not** trigger a schema bump (see Point below).
3. **Results:** The third number will be incremented each time the results are no longer directly comparable to previous versions.
4. **Point:** The fourth number will be incremented any other time the pipeline code changes in a manner that doesn't meet the criteria above, such as changes that impact performance but not results; changes to documentation; options that are off by default; and new output files that don't interfere with existing outputs (including additions to the `expected-outputs-*` lists in `pyproject.toml`).

**Experimental outputs:** The `experimental/` and `experimental_downstream/` output directories are exempt from the above guarantees. Files in these directories may change or be removed in any release, including point releases. See [output.md](output.md) for details.

Users relying on pipeline outputs should take the following actions in response to pipeline changes:

1. **Point change or higher:** Review changes for new outputs and options that could be relevant to user's application.
2. **Results change or higher:** Review changes for effects on output interpretation; show caution in comparing outputs across the version boundary.
3. **Schema change or higher:** Update downstream code to reflect new output schema.
4. **Major change:** Review new code and outputs thoroughly; be prepared for major changes to downstream code.

Note that, when a version has not yet been merged to `master`, it should have the suffix `-dev`.

A downstream fork that modifies a release marks its own versions with a [local version label](https://peps.python.org/pep-0440/#local-version-identifiers): `+` and dot-separated letters and digits after the release it modifies, as in `3.4.0.1+mclaren2026.1`. Use `+`, not a hyphen: `pyproject.toml`'s version must be valid PEP 440 for the project to install. The pipeline's version-compatibility check ignores both kinds of suffix and compares only the four numbers. Increase the label's last number (`+mclaren2026.1` to `+mclaren2026.2`) for any change that could alter a run's outputs, and before any run whose results will be reported; other changes needn't bump it.

## Dev versions and stacking

A `-dev` version represents all unreleased work since the last release. For example, if the last release was `3.2.0.2`, then the current dev version is `3.2.0.3-dev` — representing the next point release in progress. Multiple PRs that each make point-level changes all share the same `-dev` version; they do not each increment the point component. A PR should only change the version from an existing `-dev` version if its changes justify a **higher-level** bump (e.g. point → results would change `3.2.0.3-dev` to `3.2.1.0-dev`).

When updating the pipeline version, both the version in `pyproject.toml` and the topmost heading in `CHANGELOG.md` must be updated to match. See [developer.md](developer.md) for CHANGELOG formatting guidelines.
