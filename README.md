# LDAM-Ferroptosis-TLE

Analysis code for the manuscript:

**"Lipid-droplet-accumulating microglia adopt an iron-storing, ferroptosis-resistant transcriptional state in temporal lobe epilepsy: a multi-omics and spatial transcriptomic study"**

Zhenglin Liu, Jian Yin*
Department of Neurosurgery, The Second Affiliated Hospital of Dalian Medical University, Dalian, China
*Corresponding author: jyin@dmu.edu.cn | ORCID: 0000-0002-5167-9821

## Overview

This repository contains all analysis code used in the study, covering:

- Bulk transcriptome differential expression, WGCNA, and machine-learning signature construction (13 GEO cohorts)
- Single-cell / single-nucleus RNA-seq analysis of microglial states (LDAM scoring)
- Direction-resolved ferroptosis scoring (driver / suppressor / net, FerrDb-annotated)
- Donor-level pseudobulk statistics
- Genetic causal inference: MAGMA gene-set analysis (Python reimplementation), scDRS, gsMap, SMR/HEIDI, colocalization, cis-pQTL Mendelian randomization, TWAS
- Spatial transcriptomics: human hippocampus (GSE288888) and rat epileptogenesis time series (GSE301260), GraphST, cell2location, COMMOT
- Extracellular vesicle spatial co-localization
- Contamination-gated diagnostic signature (S100A10/CD9/SYNGR2) and six-gene blood panel with external validation

## Repository structure

```
R/        Standardized R scripts (45 analysis scripts + _config.R + README_code.md)
python/   Python scripts for genetic, spatial, and single-cell analyses (82 scripts)
```

## Requirements

- **R** >= 4.3 (see `R/README_code.md` for per-script package dependencies)
- **Python** >= 3.10 with: scanpy, pandas, numpy, scipy, scikit-learn, matplotlib, torch (CPU), cell2location, commot, xgboost, shap

## Data availability

All input data are publicly available from GEO, the ILAE Consortium on Complex Epilepsies, IEU OpenGWAS, GTEx, MetaBrain, eQTLGen, OneK1K, and deCODE. Accession numbers are listed in the manuscript Methods. Raw data are not redistributed here.

Set the environment variable `LDAM_BASE` to your local data directory before running the R scripts (see `R/_config.R`).

## Citation

If you use this code, please cite the manuscript (citation details will be added upon publication).[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23480511.svg)](https://doi.org/10.5281/zenodo.23480511)

## License

MIT License (see LICENSE).
