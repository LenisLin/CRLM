# Peritumoral Niches Stratify Outcomes after Adjuvant Therapy in Colorectal Liver Metastases

This repository contains 35 R/Python scripts for the CRLM spatial
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

The processed-data package uses Zenodo concept DOI
[`10.5281/zenodo.17638057`](https://doi.org/10.5281/zenodo.17638057) and is
organized as nine files:

1. `README.md`: package overview and reading instructions.
2. `LICENSE`: Creative Commons Attribution 4.0 notice for the data package.
3. `metadata.json`: Zenodo dataset metadata.
4. `file_manifest.tsv`: canonical nine-file inventory.
5. `data_dictionary.md`: object schemas, units, identifiers, and worksheet contents.
6. `FDZS-1_IMC/FDZS1_IMC_processed.rds`: 35 markers x 2,018,260 IM/PT/TC cells from 35 patients and 311 ROIs.
7. `FDZS-3_scRNA/FDZS3_scRNA_processed.h5ad`: 35,900 scRNA-seq cells with processed and raw count representations.
8. `FDZS-4_ST/FDZS4_ST_processed.h5ad`: 17,902 Bin100 spots x 56,695 features with 30 normalized RCTD abundance columns.
9. `supplementary_tables/Supplementary_Tables_public.xlsx`: Supplementary Tables 1-4 with patient, ROI, acquisition-channel, and reagent records.

The public fields, object dimensions, and data interfaces are specified in
[docs/data_contract.md](docs/data_contract.md).

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
