# ui_qc.R — QC module UI (no library calls here)

qc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        width = 12, title = "QC Controls", status = "primary", solidHeader = TRUE,
        div(
          style = "display:flex; gap:10px; flex-wrap:wrap;",
          actionButton(ns("run_qc"),   "Run QC analysis",  icon = icon("play")),
          actionButton(ns("show_vln"), "Plot QC (violin)", icon = icon("chart-area")),
          downloadButton(ns("dl_vln"), "Download QC violin"),
          actionButton(ns("show_sct"), "Plot QC scatter",  icon = icon("dot-circle")),
          downloadButton(ns("dl_sct"), "Download QC scatter")
        ),
        br(),
        uiOutput(ns("qc_status"))
      )
    ),
    fluidRow(
      box(
        width = 12, title = "QC Violin (nFeature, nCount, % mito, % ribo)",
        status = "primary", solidHeader = TRUE,
        plotOutput(ns("vln_plot"), height = 480)
      )
    ),
    fluidRow(
      box(
        width = 12, title = "QC Scatter (UMIs vs % mitochondrial)",
        status = "primary", solidHeader = TRUE,
        plotOutput(ns("sct_plot"), height = 520)
      )
    )
  )
}
