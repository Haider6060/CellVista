# App.R — CellPrint: Data Upload + QC + Clustering + Pseudotime + CellEntropy + RTI + BSI
# (All libraries are loaded in global.R)

# ==== Sources ====
source("global.R")                  # libraries, options, helpers

# Core modules
source("ui_qc.R")                   # QC module UI
source("server_qc.R")               # QC module server
source("ui_clustering.R")           # Clustering module UI
source("server_clustering.R")       # Clustering module server

# Pseudotime
source("ui_pseudotime.R")           # Pseudotime module UI
source("server_pseudotime.R")       # Pseudotime module server

# CellEntropy
source("ui_cellentropy.R")          # CellEntropy module UI
source("server_cellentropy.R")      # CellEntropy module server

# RTI
source("rti_ui.R")                  # RTI module UI
source("rti_server.R")              # RTI module server

# BSI (Boundary Sharpness Index) — NEW
source("BSI_ui.R")                  # BSI module UI
source("BSI_server.R")              # BSI module server


# Main UI (fallback to ui.R if main_ui.R not found)
if (file.exists("main_ui.R")) {
  source("main_ui.R")               # defines `ui`
} else if (file.exists("ui.R")) {
  source("ui.R")                    # defines `ui`
} else {
  stop("Neither main_ui.R nor ui.R found. Please add your main UI file.")
}

# ==== Options ====
if (!exists("MAX_UPLOAD_MB", inherits = FALSE)) MAX_UPLOAD_MB <- 1024
options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)
options(shiny.fullstacktrace = TRUE)

# ==== Helpers (define only if absent) ====
if (!exists("%||%", inherits = FALSE)) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
if (!exists("infer_dataset_type", inherits = FALSE)) {
  infer_dataset_type <- function(filepath, obj) {
    fname <- basename(filepath %||% "")
    if (nzchar(fname) && grepl("^GSE\\d+", fname, ignore.case = TRUE)) return("GEO Dataset")
    if (!is.null(obj@misc) && !is.null(obj@misc$dataset_type)) return(as.character(obj@misc$dataset_type))
    "User Upload"
  }
}
if (!exists("safe_seurat_summary", inherits = FALSE)) {
  safe_seurat_summary <- function(obj) {
    if (is.null(obj)) return("No object loaded.")
    paste(capture.output(print(obj)), collapse = "\n")
  }
}

# ==== Server ====
server <- function(input, output, session) {
  seurat_obj   <- reactiveVal(NULL)
  dataset_type <- reactiveVal("—")
  dataset_name <- reactiveVal("—")
  
  observeEvent(input$rds_file, {
    req(input$rds_file)
    
    if (!is.null(input$rds_file$size) &&
        input$rds_file$size > (MAX_UPLOAD_MB * 1024^2)) {
      showNotification(
        sprintf("File exceeds %s MB limit. Please upload a smaller file or increase MAX_UPLOAD_MB.", MAX_UPLOAD_MB),
        type = "error", duration = NULL
      )
      return(invisible(NULL))
    }
    
    fpath <- input$rds_file$datapath
    dataset_name(basename(input$rds_file$name))
    
    obj <- tryCatch(readRDS(fpath),
                    error = function(e) {
                      showNotification(paste("Failed to read RDS:", e$message),
                                       type = "error", duration = NULL)
                      NULL
                    })
    
    if (is.null(obj) || !inherits(obj, "Seurat")) {
      showNotification("Uploaded file is not a valid Seurat object (.rds).", type = "error", duration = NULL)
      return(invisible(NULL))
    }
    
    seurat_obj(obj)
    dataset_type(infer_dataset_type(input$rds_file$name, obj))
  }, ignoreInit = TRUE, priority = 10)
  
  output$upload_status <- renderUI({
    if (is.null(input$rds_file)) {
      tags$span("No file uploaded yet.")
    } else {
      sz_mb <- round((input$rds_file$size %||% 0) / 1024^2, 2)
      tags$span(style = "color:#2ca02c;",
                sprintf("Upload complete — %s (%.2f MB)", input$rds_file$name, sz_mb))
    }
  })
  output$seurat_summary <- renderText(safe_seurat_summary(seurat_obj()))
  output$dataset_type   <- renderText(dataset_type())
  output$dataset_name   <- renderText(dataset_name())
  
  output$metadata_table <- renderDT({
    obj <- seurat_obj(); req(obj)
    meta <- tryCatch({
      md <- obj@meta.data
      if (!is.data.frame(md)) md <- as.data.frame(md)
      if (nrow(md) == 0) data.frame(message = "No metadata available in object.") else md
    }, error = function(e) data.frame(message = "No metadata available in object."))
    n <- input$n_preview %||% 20
    DT::datatable(
      head(meta, n),
      rownames = TRUE,
      options = list(pageLength = n, lengthChange = FALSE, searching = TRUE, scrollX = TRUE)
    )
  })
  
  # ==== Modules ====
  qc_server("qc", seurat_obj_reactive = seurat_obj)
  clustering_server("clustering", seurat_obj_reactive = seurat_obj)
  
  # Pseudotime (placed after clustering, before CellEntropy)
  pseudotime_server("pseudotime", seurat_obj_reactive = seurat_obj)
  
  cellentropy_server("cellentropy", seurat_obj_reactive = seurat_obj)
  rti_server("rti", seurat_obj_reactive = seurat_obj, npcs = 20)  # RTI module (same dims as clustering by default)
  
  # BSI module server hookup (NOTE: BSI_server expects obj_r = reactive Seurat object)
  BSI_server("bsi", obj_r = seurat_obj)
}

# ==== Launch App ====
shinyApp(ui = ui, server = server)
