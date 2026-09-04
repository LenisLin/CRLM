# SpMap

SpMap is a four-class H&E tile classifier for reference-tile construction,
training-tile preprocessing, structure-aware splitting, CUDA MLP training,
WSI inference, and H&E-derived ISR.

## Scientific contract

SpMap uses non-overlapping 256 x 256 H&E tiles and four fixed classes:

| Class ID | Label |
|---:|---|
| 0 | TC |
| 1 | PIR |
| 2 | PSM |
| 3 | OTHER_PT |

Reference PIR/PSM structures are defined by registered IMC cells with explicit
structure IDs. Supply selected TC and OTHER_PT source ROIs in the tile manifest.
`build_structure_bboxes` accepts registered-cell coordinates, image size, and
structure, image, and class columns; `build_complete_tile_manifest` accepts the
selected structure and ROI regions.

The final public model uses paired 768- and 512-dimensional CONCH
representations (`768 + 512 = 1280` input) in a `1280 -> 450 -> 4` MLP with
ReLU, dropout 0.25, and no BatchNorm. Checkpoints are selected by validation
accuracy. The public final-model archive contains exactly five `primary_1280`
fold checkpoints, one per fold, with portable manifests. Training uses
inverse-frequency weighted cross-entropy, AdamW, learning rate 0.001, weight
decay 0.0001, batch size 256, 100 epochs, and seed 24. Splitting uses seed 42,
an explicit structure/source-tile `parent_id`, an 80% development and 20%
common internal holdout split, and five-fold `StratifiedGroupKFold` inside
development.

For patient-level analysis, tile predictions are aggregated within each WSI to
produce one H&E-derived ISR. The score is the predicted PIR tile count divided
by the predicted PSM tile count and is defined only when the predicted PSM tile
count is nonzero.

## Modules and execution order

1. `reference_tiles.py` converts caller-selected registered structures or ROIs
   into clipped H&E regions and complete non-overlapping tile manifests.
2. `preprocessing.py` applies training-side luminosity standardization,
   Vahadane normalization, and the three named offline augmentation policies.
3. `conch_features.py` extracts paired 768/512 CONCH v1 representations on
   CUDA. Supply the CONCH checkpoint identifier, `model_loader`, RGB
   `preprocess` transform, and `forward_representations` adapter.
4. `splits.py` requires explicit, single-class `parent_id` groups and returns
   development, common-holdout, and five train/validation partitions.
5. `mlp.py` preloads all required paired features from the canonical
   uncompressed tar-shard layout into one contiguous RAM matrix before any
   fold or epoch loop. It provides the canonical MLP, training, probability
   output, metrics, and checkpoint-manifest helpers.
6. `wsi.py` reads a selected OpenSlide level on a complete 256-pixel
   grid, filters tiles with strictly more than 80% white observed pixels,
   applies the supplied Macenko normalizer, and zero-pads retained border
   regions on the right and bottom. Supply the WSI level and
   `macenko_normalizer` callable.
7. `inference.py` batches retained WSI tiles through the paired CONCH interface
   and a selected MLP checkpoint. Supply the checkpoint path and `model_id`; it
   returns coordinates, model identity, four named probabilities, and argmax
   class.
8. `isr.py` returns one row per WSI with all four class counts,
   `total_256_tiles`, and ISR.

The tar archive root is supplied by configuration. Its expected layout is:

```text
<feature_archive_root>/
|-- tc/{768,512}/part*.tar
|-- pir_psm/{768,512}/part*.tar
`-- other_pt/{768,512}/part*.tar
```

Archive members are NPZ files keyed by the split manifest's explicit tile key.
Construct exactly one `InMemoryFeatureStore` from all fold and holdout manifests
before entering the fold loop, then pass that store to every `train_fold` call.
Training and inference use CUDA. Key package versions are listed under `spmap`
in `envs/package_versions.yml`.

## Required runtime inputs

For WSI inference, record the selected WSI level or magnification, Macenko
normalizer, CONCH checkpoint identifier, CONCH loader and representation
adapters, and MLP `model_id`. The iterator evaluates white fraction on observed
pixels, normalizes the observed region, then zero-pads the right and bottom
edge to 256 x 256 pixels.

Core split and training modules use Python, NumPy, pandas, scikit-learn, and
PyTorch. Vahadane preprocessing uses `staintools` and OpenCV; stochastic
augmentation uses Albumentations and scikit-image; WSI reading uses OpenSlide.
