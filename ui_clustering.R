# ui_clustering.R
# Minimal, universal clustering UI (uses uploaded object from Data tab — no re-upload)

clustering_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        width = 12, title = "Clustering & UMAP", status = "primary", solidHeader = TRUE,
        div(
          class = "d-flex",
          sliderInput(ns("resolution"), "Clustering Resolution:",
                      min = 0.2, max = 1.2, value = 0.5, step = 0.1, width = "400px"),
          # Button changed to default (white/gray) — no blue
          actionButton(ns("run"), "Run Clustering and UMAP",
                       class = "btn btn-default", style = "margin-left:16px;"),
          downloadButton(ns("dl_umap"), "Download UMAP",
                         style = "margin-left:16px;")
        ),
        tags$hr(),
        plotOutput(ns("umap"), height = "600px")
      )
    )
  )
}
