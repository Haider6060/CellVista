# ui_cellentropy.R

cellentropy_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("CellEntropy"),
    fluidRow(
      column(
        12,
        actionButton(ns("run_entropy"), "Run Entropy Analysis"),
        span(id = ns("done_msg"), style = "margin-left:12px; font-weight:600;")
      )
    ),
    tags$hr(),
    
    h4("Bar Plot: Mean Entropy per Cell Type"),
    fluidRow(
      column(
        12,
        actionButton(ns("plot_bar"), "Plot Bar"),
        downloadButton(ns("download_bar"), "Download Bar plot"),
        br(), br(),
        plotOutput(ns("bar_plot"), height = "460px")
      )
    ),
    tags$hr(),
    
    h4("UMAP: Entropy Landscape with Cell-Type Labels"),
    fluidRow(
      column(
        12,
        actionButton(ns("plot_umap"), "Plot UMAP"),
        downloadButton(ns("download_umap"), "Download UMAP plot"),
        br(), br(),
        plotOutput(ns("umap_plot"), height = "520px")
      )
    ),
    tags$hr(),
    
    h4("Gene Entropy Tables"),
    fluidRow(
      column(
        12,
        # changed to downloads so the browser asks where to save
        downloadButton(ns("download_csv_all"), "Download CSV: All Genes (per type)"),
        downloadButton(ns("download_csv_top"), "Download CSV: Top-50 High-Entropy (per type)")
      )
    ),
    tags$hr(),
    
    h4("Heatmap: Top-30 High-Entropy Genes"),
    fluidRow(
      column(
        12,
        actionButton(ns("plot_heatmap"), "Plot Heatmap"),
        downloadButton(ns("download_heatmap"), "Download Heatmap "),
        br(), br(),
        plotOutput(ns("heatmap_plot"), height = "700px")
      )
    )
  )
}
