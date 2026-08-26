# Peritumoral Niches Stratify Outcomes after Adjuvant Therapy in Colorectal Liver Metastases

This repository contains 38 R/Python scripts for the CRLM spatial
microenvironment study. Figure entry points are organized by main figure, and
shared functions are organized by modality.

## Repository layout

```text
config/       Path-reference worksheet plus shared analysis parameters
docs/         Data, execution, figure/table, and provenance manifests
envs/         Key runtime and package-version reference
figures/      Figure 1-7 entry points and figure-specific documentation
src/          Shared IMC, SpMap, scRNA-seq, and ST functions
```

Each `figures/FigureN/README.md` records the covered panels, execution order,
required inputs, outputs, analysis units, and package-version family.

## Processed data

The data package uses Zenodo concept DOI
[`10.5281/zenodo.17638057`](https://doi.org/10.5281/zenodo.17638057). It contains
the processed FDZS-1 IMC, FDZS-3 scRNA-seq, and FDZS-4 spatial-transcriptomics
objects; Supplementary Tables 1-4; and derived patient-level ISR, ROI-level
cell-composition, and SpMap artifacts. The package `file_manifest.tsv` provides
the complete inventory.

The public fields, object dimensions, and data interfaces are specified in
[docs/data_contract.md](docs/data_contract.md). The
[public-artifact index](docs/public_artifacts.md) maps each derived data or model
artifact to its owning analysis and export code.

## Analysis organization

Run figure entry points from the repository root in the order recorded in each
Figure README. Before direct execution, edit the paths in the script's top-level
configuration and enable its run switch where present. Entry points take no
command-line flags or positional arguments. `config/paths.example.yml` is a
path-reference worksheet for matching the same inputs across scripts; each
script reads its own editable configuration.
`config/analysis_parameters.yml` records shared fixed values, and
`envs/package_versions.yml` records key versions for each analysis family.

The repository contains analytical code for IMC, SpMap, scRNA-seq, and spatial
transcriptomics stages. Figure 7 scripts use processed MatrixMarket bundles,
registered SpMap tile annotations, and offline pathway resources specified in
the path configuration. Figure 6 uses prepared discovery and test patient
tables. Each script documents its required columns and output tables.

The machine-readable [figure/table manifest](docs/figure_table_manifest.tsv)
lists panel and table coverage. The [provenance manifest](docs/provenance_manifest.tsv)
maps analysis code and public data objects to their inputs and outputs.

## License and citation

Code is available under the [MIT License](LICENSE). Cite the software metadata
in [CITATION.cff](CITATION.cff) and the version-specific Zenodo DOI associated
with the processed-data package.
