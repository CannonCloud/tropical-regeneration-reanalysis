# matters-arising-williams-2024

Reanalysis code and manuscript for a Matters Arising responding to Williams et al. (2024), "Global potential for natural regeneration in deforested tropical regions," *Nature* 636, 131–137. doi: [10.1038/s41586-024-08106-4](https://doi.org/10.1038/s41586-024-08106-4)

Reanalyzing data from Brazil (20.3% of the global estimate), this work shows that the published ~215 Mha estimate is driven by three independent methodological sensitivities: training prevalence assumptions, spatial sampling design, and predictor encoding. Correcting the prevalence imbalance alone shrinks predicted area from 25.8 Mha to 4.8 Mha. Priority rankings are also unstable across sampling designs, and the biophysical-only model cannot distinguish ecological limitation from ongoing human pressure.

A pre-submission version of the Matters Arising was evaluated by [The Unjournal](https://unjournal.pubpub.org/pub/evalsumnaturalregeneration/).

**Author:** Cannon Cloud, WZB Berlin Social Science Center / Goethe University Frankfurt

---

## Repository Structure

```
matters-arising-williams-2024/
│
├── manuscript/
│   ├── matters_arising.tex          # Main manuscript (LaTeX)
│   ├── matters_arising_main.bib     # Main text references
│   └── matters_arising_si.bib       # Supplementary references
│
├── python/
│   ├── pantropical/                 # Pantropical PNV sampling pipeline (GEE)
│   │   ├── step1a_submit_batches.py
│   │   ├── step1b_merge_batches.py
│   │   ├── step2_sample_layers.py
│   │   ├── step3_grid_distances.py
│   │   └── step4_merge_all.py
│   │
│   └── brazil/                      # Brazil reanalysis pipeline (GEE)
│       ├── brazil_step0_export_masks.py
│       ├── brazil_step0b_export_annuli.py
│       ├── brazil_step1_sample_points.py
│       ├── brazil_step2_sample_layers.py
│       ├── brazil_step3_forest_variables.py
│       └── brazil_step4_merge.py
│
├── R/
│   ├── 01_train_models.R            # Train random forest models
│   ├── 01_correlations_figure.R     # Correlation heatmaps
│   ├── 02_figures.R                 # Main manuscript figures
│   ├── 03_model_performance.R       # OOB accuracy and AUC
│   └── 04_socioeconomic_models.R    # Socioeconomic extension models
│
├── figures/                         # Output figures (PDF)
│
├── SEQUENTIAL_WORKFLOW_PNV_SAMPLING.md   # Pantropical pipeline guide
├── brazil_pipeline_workflow.md           # Brazil pipeline guide
├── R_WORKFLOW.md                         # R analysis guide
└── README.md
```

---

## Reproducing the Analysis

### 1. Python / GEE pipelines

Two separate Google Earth Engine pipelines build the training datasets. See their dedicated guides:

- **Pantropical correlation analysis:** [`SEQUENTIAL_WORKFLOW_PNV_SAMPLING.md`](SEQUENTIAL_WORKFLOW_PNV_SAMPLING.md)
- **Brazil reanalysis:** [`brazil_pipeline_workflow.md`](brazil_pipeline_workflow.md)

Both require a GEE account.

### 2. R analysis

See [`R_WORKFLOW.md`](R_WORKFLOW.md). Run scripts in numbered order, starting from `01_train_models.R`. Requires R ≥ 4.2 and the packages listed in that file.

---

## Data

The Williams et al. prediction rasters and PNV masks are available from their [Zenodo repository](https://doi.org/10.5281/zenodo.13168079). The Fagan et al. regrowth classification is available from [Global Forest Watch](https://www.globalforestwatch.org/). All other data sources are listed in the manuscript Supplementary Methods (SM1).

This repository does not include large raster files. The `figures/` directory contains the PDF outputs referenced in the manuscript.

---

## Citation

> Cloud, C. (2025). Training prevalence and spatial sampling design drive aggregate estimates of tropical forest regeneration potential. *Matters Arising*, responding to Williams et al. (2024) *Nature* 636, 131–137.

---

## License

Code: [MIT License](LICENSE)
