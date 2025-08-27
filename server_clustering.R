# server_clustering.R
# Minimal, universal clustering server (consumes existing Seurat object from Data tab)
# Usage from app.R: clustering_server("clustering", seurat_obj_reactive = seurat_obj)

clustering_server <- function(id, seurat_obj_reactive, nfeatures = 2000, npcs = 20) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(proc = NULL)
    
    # Reset processed object when source object changes
    observeEvent(seurat_obj_reactive(), {
      rv$proc <- NULL
    }, ignoreInit = TRUE)
    
    # --- Helpers ---
    set_active_assay <- function(obj) {
      if ("RNA" %in% names(obj@assays)) DefaultAssay(obj) <- "RNA"
      else DefaultAssay(obj) <- names(obj@assays)[1]
      obj
    }
    add_percent_mt <- function(obj) {
      assay_now <- DefaultAssay(obj)
      try_patterns <- c("^(MT-|mt-|Mt-)", "^mt\\.", "^MT\\.")
      added <- FALSE
      for (pat in try_patterns) {
        ok <- tryCatch(any(grepl(pat, rownames(obj[[assay_now]]))), error = function(e) FALSE)
        if (isTRUE(ok)) {
          obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = pat, assay = assay_now)
          added <- TRUE
          break
        }
      }
      if (!added) obj[["percent.mt"]] <- 0
      obj
    }
    
    # --- Run universal pipeline on click ---
    observeEvent(input$run, {
      obj <- seurat_obj_reactive()
      validate(need(!is.null(obj) && inherits(obj, "Seurat"),
                    "Please upload a valid Seurat object in the Data tab first."))
      
      withProgress(message = "Clustering…", value = 0, {
        obj <- set_active_assay(obj);                incProgress(0.1)
        obj <- add_percent_mt(obj);                  incProgress(0.2)
        obj <- NormalizeData(obj, verbose = FALSE);  incProgress(0.35)
        obj <- FindVariableFeatures(obj, selection.method = "vst",
                                    nfeatures = max(200, min(nfeatures, nrow(obj))),
                                    verbose = FALSE);                         incProgress(0.5)
        vars <- intersect("percent.mt", colnames(obj@meta.data))
        obj <- ScaleData(obj, features = VariableFeatures(obj),
                         vars.to.regress = vars, verbose = FALSE);           incProgress(0.65)
        obj <- RunPCA(obj, features = VariableFeatures(obj),
                      npcs = max(10, npcs), verbose = FALSE);                incProgress(0.8)
        use_dims <- 1:max(10, npcs)
        obj <- FindNeighbors(obj, dims = use_dims, verbose = FALSE)
        obj <- RunUMAP(obj, dims = use_dims, verbose = FALSE)
        obj <- FindClusters(obj, resolution = as.numeric(input$resolution), verbose = FALSE)
        rv$proc <- obj;                                                     incProgress(1)
      })
    })
    
    # --- Plot ---
    umap_plot <- reactive({
      req(rv$proc)
      DimPlot(rv$proc, reduction = "umap", label = TRUE, repel = TRUE) +
        theme(legend.position = "right") +
        guides(color = guide_legend(title = "Cluster"))
    })
    
    output$umap <- renderPlot({
      umap_plot()
    }, res = 120)
    
    # --- Download (PNG, 1000 dpi) ---
    output$dl_umap <- downloadHandler(
      filename = function() {
        paste0("umap_clusters_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
      },
      content = function(file) {
        g <- umap_plot()
        ggplot2::ggsave(filename = file, plot = g, width = 6, height = 5, dpi = 1000)
      }
    )
    
    # Return processed object if the app wants it
    invisible(list(seurat = reactive(rv$proc)))
  })
}
