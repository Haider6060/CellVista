# UI module for Cellular Regulatory Turbulence Index (RTI) + Early-warning + Cell types
rti_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        width = 3,
        wellPanel(
          h4("RTI Analysis"),
          actionButton(ns("run"), "Run RTI"),
          br(), br(),
          numericInput(ns("npcs"), "PCA dimensions (same as clustering)", value = 20, min = 10, step = 1),
          helpText("Uses the existing SNN from clustering. If missing, neighbors will be recomputed with these dims."),
          tags$hr(),
          h4("Early-warning (branch commitment)"),
          numericInput(ns("ew_bins"),   "Bins along pseudotime", value = 50, min = 20, step = 5),
          numericInput(ns("ew_tau"),    "Label-entropy threshold (τ)", value = 0.20, min = 0, max = 1, step = 0.01),
          numericInput(ns("ew_win"),    "Consecutive low-entropy bins", value = 5, min = 3, step = 1),
          numericInput(ns("ew_min_n"),  "Min. cells per bin", value = 50, min = 1, step = 1),
          actionButton(ns("run_ew"), "Compute early-warning"),
          helpText("Commitment = first bin where neighbor label-entropy ≤ τ for the given number of consecutive bins, with ≥ min cells/bin. Peak RTI is searched before commitment; lead = commitment − peak."),
          tags$hr(),
          h4("Cell type annotation (SingleR)"),
          actionButton(ns("run_types"), "Annotate cell types"),
          helpText("Automatically labels cells with SingleR using celldex references (HumanPrimaryCellAtlas/MouseRNAseq; chosen by gene overlap). Adds 'cell_type' and 'major_type' to metadata and enables type-level plots and CSVs.")
        )
      ),
      column(
        width = 9,
        tabsetPanel(
          tabPanel(
            "RTI UMAP",
            plotOutput(ns("umap_rti"), height = 500),
            downloadButton(ns("dl_umap_rti"), "Download UMAP")
          ),
          tabPanel(
            "RTI by Cluster",
            plotOutput(ns("vln_rti"), height = 500),
            downloadButton(ns("dl_vln_rti"), "Download Violin plot")
          ),
          tabPanel(
            "Cluster Summary",
            tableOutput(ns("tbl_summary")),
            downloadButton(ns("dl_summary_csv"), "Download CSV")
          ),
          tabPanel(
            "Early-warning",
            plotOutput(ns("ew_plot"), height = 520),
            fluidRow(
              column(6, downloadButton(ns("dl_ew_plot"), "Download Early-warning plot")),
              column(6, downloadButton(ns("dl_ew_csv"),  "Download Early-warning summary"))
            ),
            br(),
            tableOutput(ns("ew_summary"))
          ),
          tabPanel(
            "Cell types",
            br(),
            strong("UMAP colored by cell type"),
            plotOutput(ns("umap_celltype"), height = 480),
            downloadButton(ns("dl_umap_celltype"), "Download UMAP by cell type"),
            tags$hr(),
            fluidRow(
              column(
                width = 6,
                strong("RTI by cell type"),
                plotOutput(ns("vln_rti_celltype"), height = 420),
                downloadButton(ns("dl_vln_rti_celltype"), "Download violin (cell type)")
              ),
              column(
                width = 6,
                strong("RTI by major family"),
                plotOutput(ns("vln_rti_majortype"), height = 420),
                downloadButton(ns("dl_vln_rti_majortype"), "Download violin (major type)")
              )
            ),
            tags$hr(),
            fluidRow(
              column(
                width = 6,
                strong("Type-level summary (n, median, IQR)"),
                tableOutput(ns("tbl_type_summary")),
                downloadButton(ns("dl_type_summary_csv"), "Download type summary CSV")
              ),
              column(
                width = 6,
                strong("Cluster ↔ type mapping (max-overlap)"),
                tableOutput(ns("tbl_cluster_type_map")),
                downloadButton(ns("dl_cluster_type_map_csv"), "Download cluster↔type CSV")
              )
            )
          )
        )
      )
    )
  )
}
