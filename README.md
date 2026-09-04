# Peritumoral Niches Stratify Outcomes after Adjuvant Therapy in Colorectal Liver Metastases

This repository contains 44 R/Python scripts for the CRLM spatial
microenvironment study. Figure entry points are organized by main figure, and
shared functions are organized by modality.

## Release 2.0.0

Released 2026-09-04. This release documents the Figure 1B supplementary-table
producer, the Figure 6E-F continuous-ISR multivariable Cox entry point, and
the final five-checkpoint `primary_1280` SpMap model archive.

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

Run figure entry points from the repository root after editing the paths in the
script's top-level configuration and enabling its run switch where present.
Entry points take no command-line flags or positional arguments.
`config/paths.example.yml` provides shared path references; each script reads
its own editable configuration.
`config/analysis_parameters.yml` records shared fixed values, and
`envs/package_versions.yml` records key versions for each analysis family.

The repository contains analytical code for IMC, SpMap, scRNA-seq, and spatial
transcriptomics stages. Figure 7 includes the full MatrixMarket workflow and a
direct plot entry point for the released FDZS-3 H5AD; later spatial stages use
registered SpMap tile annotations and offline pathway resources. Figure 1B
reads Supplementary Tables 1 and 5. Figure 6E-F reads those tables together
with the FDZS-1 patient-level H&E-derived ISR table. Each script documents its
required fields and output tables.

The machine-readable [figure/table manifest](docs/figure_table_manifest.tsv)
lists panel and table coverage. The [provenance manifest](docs/provenance_manifest.tsv)
maps analysis code and public data objects to their inputs and outputs.

## License and citation

Code is available under the [MIT License](LICENSE). Cite the software metadata
in [CITATION.cff](CITATION.cff) and cite the processed data through the Zenodo
concept DOI above.
