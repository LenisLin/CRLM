# Public analysis artifacts

The Zenodo package associated with concept DOI `10.5281/zenodo.17638057`
contains the processed native objects and the derived artifacts below. Paths are
the published Zenodo file names.

| Public artifact | Analysis unit | Owning public code |
|---|---|---|
| `FDZS1_patient_level_ISR.tsv` | FDZS-1 patient | `figures/Figure6/03_prepare_public_isr_tables.py`; `figures/Figure6/04_continuous_isr_multivariable_cox.R`; Supplementary Figure 11D consumer |
| `FDZS2_patient_level_ISR_RFS.tsv` | FDZS-2 patient | `figures/Figure6/03_prepare_public_isr_tables.py`; `figures/Figure6/02_cohort_and_rfs.R` |
| `SpMap_patient_cohort_assignment.tsv` | Patient | `figures/Figure6/03_prepare_public_isr_tables.py` |
| `ISR_definition_and_cutoff.tsv` | Score definition | `figures/Figure6/03_prepare_public_isr_tables.py`; `figures/Figure6/02_cohort_and_rfs.R` |
| `FDZS1_ROI_cell_composition.tsv.gz` | Patient-ROI-tissue | `figures/Figure4/02_export_roi_cell_composition.R` |
| `SpMap_reference_tile_labels.tsv.gz` | Reference tile | `figures/Figure5/03_prepare_public_spmap_artifacts.py` |
| `SpMap_tile_partitions.tsv.gz` | Reference tile/parent group | `figures/Figure5/03_prepare_public_spmap_artifacts.py`; `src/spmap/splits.py` |
| `SpMap_model_performance_5seeds.tsv` | Seed/class/metric | `figures/Figure5/03_prepare_public_spmap_artifacts.py`; Supplementary Figure 11B-C consumer |
| `SpMap_confusion_matrix.tsv` | True/predicted class | `figures/Figure5/03_prepare_public_spmap_artifacts.py`; Supplementary Figure 11B consumer |
| `SpMap_OOF_predictions_C10_5seeds_5folds.tar.gz` | Seed-fold tile prediction | `figures/Figure5/03_prepare_public_spmap_artifacts.py` |
| `SpMap_model_weights_primary_1280_5folds.tar.gz` | Exactly five final `primary_1280` fold checkpoints and manifests | `figures/Figure5/03_prepare_public_spmap_artifacts.py`; `src/spmap/mlp.py` |
| `SpMap_CONCH_features.tar.gz` | Reference-tile feature | `figures/Figure5/03_prepare_public_spmap_artifacts.py`; `src/spmap/conch_features.py` |

The exact fields, identifiers, formulas, and archive contents are defined in
`docs/data_contract.md` and the Zenodo `data_dictionary.md`.
