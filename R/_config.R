# ============================================================================
# Script      : _config.R
# Purpose     : Central path configuration shared by all analysis scripts.
# Manuscript  : Methods 2.1 (reproducibility configuration)
# Input       : none (environment variable LDAM_BASE, optional)
# Output      : defines BASE_DIR / DATA_DIR / OUT_DIR / FIG_DIR / GS_DIR
# Dependencies: base R only
# Date        : standardized 2026-10
# Usage       : source("_config.R") at the top of a session, or simply run any
#               script from the project root - every script self-defines
#               BASE_DIR <- Sys.getenv("LDAM_BASE", ".") if not yet defined.
# ============================================================================

# Project root: set the environment variable LDAM_BASE, or run from the project root.
BASE_DIR <- Sys.getenv("LDAM_BASE", ".")

# Conventional sub-directories (created on demand by individual scripts)
DATA_DIR <- BASE_DIR                          # raw data files live in the project root
OUT_DIR  <- file.path(BASE_DIR, "results")    # result tables
FIG_DIR  <- file.path(BASE_DIR, "figures")    # figures
GS_DIR   <- file.path(BASE_DIR, "genesets")   # frozen LDAM / ferroptosis gene sets

# Legacy external directories referenced by two legacy scripts (07, fix_04):
# parameterised via environment variables; data are NOT part of this repository.
GLIOMA_GEO_DIR  <- Sys.getenv("GLIOMA_GEO_DIR",  file.path(BASE_DIR, "data/external/glioma/GEO"))
GLIOMA_TCGA_DIR <- Sys.getenv("GLIOMA_TCGA_DIR", file.path(BASE_DIR, "data/external/glioma/TCGA"))
