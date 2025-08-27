# ui_pseudotime.R
pseudotime_ui <- function(id) {
  ns <- NS(id)
  
  # CSS to force a dropdown caret on the multi-select "normal_vals"
  caret_css <- sprintf("
    /* make the selectize box show a clear ▼ caret */
    #%s + .selectize-control .selectize-input { position: relative; }
    #%s + .selectize-control .selectize-input::after {
      content: '\\25BE';              /* ▼ */
      position: absolute;
      right: 10px;
      top: 50%%;
      transform: translateY(-50%%);
      color: #555;
      pointer-events: none;
      font-size: 14px;
    }
  ", ns('normal_vals'), ns('normal_vals'))
  
  tagList(
    tags$head(tags$style(HTML(caret_css))),
    box(width = 12, title = "Pseudotime Trajectory Analysis", status = "primary", solidHeader = TRUE,
        fluidRow(
          column(4, actionButton(ns("run"), "Run Pseudotime Analysis")),
          column(8, div(style="margin-top:8px;color:#666;",
                        "Tip: run clustering first; this module reuses the uploaded Seurat object."))
        ),
        tags$hr(),
        
        # ---- Optional grouping override ----
        div(id = ns("group_panel"),
            h4("Optional: Fix Normal vs Tumor mapping"),
            div(style="color:#666;margin-bottom:6px;",
                "If auto-detection shows only one class, choose a column and mark which values are Normal; others will be Tumor."),
            fluidRow(
              column(4, selectInput(ns("group_col"), "Metadata column:", choices = NULL)),
              column(8, uiOutput(ns("group_values_ui")))   # rendered in server
            ),
            actionButton(ns("apply_group"), "Apply Mapping")
        ),
        
        tags$hr(),
        
        h4("UMAP (pseudotime)"),
        fluidRow(
          column(4, actionButton(ns("show_umap"), "Show UMAP")),
          column(4, downloadButton(ns("dl_umap"), "Download UMAP"))
        ),
        plotOutput(ns("plot_umap"), height = "520px"),
        
        tags$hr(),
        
        h4("Trajectory — States"),
        fluidRow(
          column(4, actionButton(ns("show_states"), "Show States")),
          column(4, downloadButton(ns("dl_states"), "Download States"))
        ),
        plotOutput(ns("plot_states"), height = "520px"),
        
        tags$hr(),
        
        h4("Trajectory — Normal vs Tumor"),
        fluidRow(
          column(4, actionButton(ns("show_group"), "Show Normal vs Tumor")),
          column(4, downloadButton(ns("dl_group"), "Download Group"))
        ),
        plotOutput(ns("plot_group"), height = "520px")
    )
  )
}
