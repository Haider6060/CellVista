# ui.R — CellVista: Data Upload + QC + Clustering + Pseudotime + CellEntropy + RTI + BSI — no library calls here

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CellVista"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload",           tabName = "data",        icon = icon("database")),
      menuItem("QC & Normalization",    tabName = "qc",          icon = icon("heartbeat")),
      menuItem("Clustering & UMAP",     tabName = "clustering",  icon = icon("project-diagram")),
      # NEW: Pseudotime (inserted after Clustering & UMAP)
      menuItem("Pseudotime",            tabName = "pseudotime",  icon = icon("hourglass-half")),
      menuItem("CellEntropy",           tabName = "cellentropy", icon = icon("bolt")),
      menuItem("RTI (Regulatory Turbulence)", tabName = "rti",   icon = icon("wind")),
      # ---------- NEW MODULE: BSI (Boundary Sharpness Index) ----------
      menuItem("BSI (Boundary Sharpness)", tabName = "bsi", icon = icon("draw-polygon"))
    )
  ),
  dashboardBody(
    # Uniform box header color + minor UI polish
    tags$head(
      tags$style(HTML("
        .box>.box-header { background-color: #3c8dbc !important; color: #fff !important; }
        .upload-note { color:#555; margin-top:6px; }
        .inline-row .shiny-text-output { display:inline; margin-right:16px; }
        .help-note { color:#444; font-size: 13px; line-height: 1.4; }
      "))
    ),
    
    tabItems(
      # --- Data Upload & Details ---
      tabItem(
        tabName = "data",
        fluidRow(
          box(
            width = 12, title = "Upload Seurat Object (.rds)", status = "primary", solidHeader = TRUE,
            fileInput("rds_file", label = NULL, accept = ".rds", buttonLabel = "Browse…"),
            div(class = "upload-note",
                sprintf("Maximum upload size: %s MB", format(MAX_UPLOAD_MB, big.mark=","))),
            uiOutput("upload_status")
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Seurat Object Information", status = "primary", solidHeader = TRUE,
            tags$pre(style = "max-height: 260px; overflow:auto; white-space: pre-wrap;",
                     textOutput("seurat_summary"))
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Dataset Summary", status = "primary", solidHeader = TRUE,
            fluidRow(
              column(6, div(class = "inline-row", strong("Type: "), textOutput("dataset_type", container = span))),
              column(6, div(class = "inline-row", strong("File: "), textOutput("dataset_name", container = span)))
            )
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Metadata Preview", status = "primary", solidHeader = TRUE,
            sliderInput("n_preview", "Rows to show:", min = 5, max = 100, value = 20, step = 5),
            DTOutput("metadata_table")
          )
        )
      ),
      
      # --- QC & Normalization (module UI; appears BEFORE clustering) ---
      tabItem(
        tabName = "qc",
        qc_ui("qc")   # defined in ui_qc.R
      ),
      
      # --- Clustering & UMAP (module UI) ---
      tabItem(
        tabName = "clustering",
        clustering_ui("clustering")   # defined in ui_clustering.R
      ),
      
      # --- Pseudotime (inserted after clustering, before CellEntropy) ---
      tabItem(
        tabName = "pseudotime",
        fluidRow(
          box(
            width = 12, title = "Pseudotime Trajectories", status = "primary", solidHeader = TRUE,
            div(class = "help-note",
                tags$ul(
                  tags$li("Uses the already loaded Seurat object (post-clustering/UMAP)."),
                  tags$li("Computes per-cell pseudotime ordering and (if supported by your module) branch/lineage assignments."),
                  tags$li("Visualizations typically include: UMAP colored by pseudotime, lineage/branch overlays, density along pseudotime, and top dynamic genes."),
                  tags$li("Provides PNG (1000 dpi) and CSV exports for pseudotime values and assignments."),
                  tags$li("This tab is wired to your module in ui_pseudotime.R / server_pseudotime.R.")
                )
            )
          )
        ),
        # Module UI (defined in ui_pseudotime.R)
        pseudotime_ui("pseudotime")
      ),
      
      # --- CellEntropy (analysis module) ---
      tabItem(
        tabName = "cellentropy",
        fluidRow(
          box(
            width = 12, title = "CellEntropy — Information-Theoretic Profiling", status = "primary", solidHeader = TRUE,
            div(class = "help-note",
                tags$ul(
                  tags$li("Uses the already uploaded dataset (no re-upload)."),
                  tags$li("Runs per-cell entropy, cell-type annotation, mean entropy per type, and gene-level entropy per type."),
                  tags$li("Generates bar plot (mean entropy by cell type), UMAP entropy map with labels,"),
                  tags$li("Saves CSVs: all genes with entropy per type, and top-50 high-entropy genes per type."),
                  tags$li("Draws heatmap of top-30 high-entropy genes (PNG downloads at 1000 dpi).")
                )
            )
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Run & Status", status = "primary", solidHeader = TRUE,
            cellentropy_ui("cellentropy")   # defined in ui_cellentropy.R; module id = "cellentropy"
          )
        )
      ),
      
      # --- RTI (Regulatory Turbulence Index) module UI ---
      tabItem(
        tabName = "rti",
        fluidRow(
          box(
            width = 12, title = "Cellular Regulatory Turbulence Index (RTI)", status = "primary", solidHeader = TRUE,
            div(class = "help-note",
                tags$ul(
                  tags$li("Runs RTI using the SAME PCA dims and SNN neighbors as clustering."),
                  tags$li("Computes TI, DD, IM, H components and aggregate RTI in [0,1]."),
                  tags$li("Shows UMAP colored by RTI, violin by cluster, and per-cluster CSV export.")
                )
            )
          )
        ),
        rti_ui("rti")  # from rti_ui.R
      ),
      
      # --- BSI (Boundary Sharpness Index) module UI ---
      tabItem(
        tabName = "bsi",
        fluidRow(
          box(
            width = 12, title = "Boundary Sharpness Index (BSI)", status = "primary", solidHeader = TRUE,
            div(class = "help-note",
                tags$ul(
                  tags$li("Click “Run BSI analysis” to compute per-cell BSI on the current embedding (prefers UMAP; falls back to PCA)."),
                  tags$li("Computes two scores: Unlabeled BSI (geometry-only kNN distance gap, k = 30) and Label-aware BSI (neighbor label disagreement with existing clusters)."),
                  tags$li("Plots (rendered only when clicked): UMAP hexbin colored by BSI (label-aware / unlabeled) and violin plot by cluster."),
                  tags$li("Downloads: Each plot is available as PNG at 1000 dpi (individual download buttons in each tab)."),
                  tags$li("This tab is wired to the BSI module files: BSI_ui.R / BSI_server.R.")
                )
            )
          )
        ),
        # Module UI (defined in BSI_ui.R)
        BSI_ui("bsi")
      )
    )
  )
)
