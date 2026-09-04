# Analysis Package Versions

`package_versions.yml` records observed runtimes and key direct package
versions for the 44 public R/Python scripts. Figure and source READMEs identify
the analysis family used by each entry point.

Run R entry points with `Rscript` and Python entry points with `python`, using
the commands in the owning Figure README. The file is a version reference, not
a solver lock file or installation specification.

The `spmap.executed_primary_1280` record preserves the Python, PyTorch, and
CUDA versions reported by the final-model manifests. The `spmap.packages`
record lists the observed `OSMOSIS_V2` runtime used by the repository commands.
