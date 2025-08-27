# global.R — CellPrint (scRNA-seq): globals, libs, options, helpers

# ==== App Settings ====
APP_NAME       <- "CellPrint"
MAX_UPLOAD_MB  <- getOption("cellprint.max_upload_mb", 1024)   # override via options()
set.seed(1234)

# Allow large uploads (set once here)
options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)

# ==== Libraries ====
library(shiny)
library(shinydashboard)
library(Seurat)
library(ggplot2)
library(dplyr)
library(DT)
library(tibble)
library(shinybusy)
library(monocle3)
# --- Added for CellEntropy module ---
library(SingleR)                     # Cell type annotation
library(celldex)                     # Reference datasets
library(SingleCellExperiment)        # SCE structure
library(SummarizedExperiment)        # colData, assays, etc.
library(scater)                      # runPCA, runUMAP, plotting helpers
library(ggridges)                    # Ridge plots
library(ggrepel)                     # Label repelling for UMAP
library(viridis)                     # Color gradients
library(gridExtra)                   # Arranging multiple plots
library(Matrix)
library(SeuratObject)
library(reshape2)
# --- Added for RTI module ---
library(igraph)                      # Graph-based pseudotime on SNN
# --- Added for BSI module ---
library(tidyr)                       # pivot_longer / data reshaping for CSVs
library(hexbin)                      # hexbin backend for ggplot2 stat_summary_hex
library(RANN)

# ==== Palette / Theme (lightweight defaults you can reuse) ====
CELLPRINT_BLUE <- "#3c8dbc"

# ==== Small Utils (used across app) ====
`%||%` <- function(a, b) if (!is.null(a)) a else b

is_seurat <- function(obj) inherits(obj, "Seurat")

safe_read_rds <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("File not found or unreadable: ", path %||% "<empty path>")
  }
  readRDS(path)
}

infer_dataset_type <- function(filepath, obj) {
  fname <- basename(filepath %||% "")
  if (nzchar(fname) && grepl("^GSE\\d+", fname, ignore.case = TRUE)) return("GEO Dataset")
  if (!is.null(obj@misc) && !is.null(obj@misc$dataset_type)) return(as.character(obj@misc$dataset_type))
  "User Upload"
}

safe_seurat_summary <- function(obj) {
  if (is.null(obj)) return("No object loaded.")
  paste(capture.output(print(obj)), collapse = "\n")
}
