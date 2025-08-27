# UI module for Boundary Sharpness Index (BSI)
BSI_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 3,
        shiny::wellPanel(
          shiny::h4("Boundary Sharpness Index (BSI)"),
          shiny::helpText(
            "BSI quantifies the sharpness of transcriptional boundaries in the embedding.",
            "Unlabeled = geometry-only; Label-aware = neighbor label disagreement."
          ),
          shiny::actionButton(ns("run"), "Run BSI analysis"),
          shiny::br(), shiny::br(),
          shiny::strong("Status:"),
          shiny::textOutput(ns("status"), inline = TRUE)
        )
      ),
      shiny::column(
        width = 9,
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Label-aware UMAP",
            shiny::actionButton(ns("plot_labelaware"), "Plot label-aware UMAP"),
            shiny::downloadButton(ns("dl_labelaware"), "Download UMAP"),
            shiny::br(), shiny::br(),
            shiny::plotOutput(ns("p_labelaware"), height = 520)
          ),
          shiny::tabPanel(
            "Unlabeled UMAP",
            shiny::actionButton(ns("plot_unlabeled"), "Plot unlabeled UMAP"),
            shiny::downloadButton(ns("dl_unlabeled"), "Download UMAP"),
            shiny::br(), shiny::br(),
            shiny::plotOutput(ns("p_unlabeled"), height = 520)
          ),
          shiny::tabPanel(
            "Violin (by cluster)",
            shiny::actionButton(ns("plot_violin"), "Plot violin"),
            shiny::div(
              style = "display:flex; gap:12px; align-items:center; margin-top:8px;",
              shiny::downloadButton(ns("dl_violin"),   "Download PLOT"),
              shiny::downloadButton(ns("dl_bsi_csv"), "Download BSI×Cluster CSV")
            ),
            shiny::br(),
            shiny::plotOutput(ns("p_violin"), height = 520)
          )
        )
      )
    )
  )
}
