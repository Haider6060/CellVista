# Server module for Boundary Sharpness Index (BSI)
# Usage (works with either arg name):
#   BSI_server("bsi", seurat_obj_reactive = seurat_obj)
#   BSI_server("bsi", obj_r = seurat_obj)
#
# Notes:
# - Expects Seurat, ggplot2, RANN, (optional) hexbin to be loaded in global.R
# - No library() calls here, per your style.

BSI_server <- function(id, seurat_obj_reactive = NULL, obj_r = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ---- pick reactive getter (support both argument names) ----
    get_obj <- if (!is.null(seurat_obj_reactive)) seurat_obj_reactive else obj_r
    if (is.null(get_obj)) stop("Provide a reactive Seurat object via `seurat_obj_reactive` or `obj_r`.")
    
    # ---- helpers ----
    choose_assay <- function(o) {
      assays <- Seurat::Assays(o)
      if ("RNA" %in% assays) return("RNA")
      assays[1]
    }
    is_normalized <- function(o, assay) {
      x <- try(Seurat::GetAssayData(o, assay = assay, slot = "data"), silent = TRUE)
      inherits(x, "dgCMatrix") && nrow(x) > 0 && length(x@x) > 0
    }
    has_scaled <- function(o, assay) {
      x <- try(o[[assay]]@scale.data, silent = TRUE)
      is.matrix(x) && nrow(x) > 0
    }
    ensure_preprocessing <- function(o) {
      # 1) set default assay (prefer RNA)
      a <- choose_assay(o)
      Seurat::DefaultAssay(o) <- a
      
      # 2) Normalize if needed
      if (!is_normalized(o, a)) {
        o <- Seurat::NormalizeData(o, normalization.method = "LogNormalize",
                                   scale.factor = 1e4, verbose = FALSE)
      }
      
      # 3) HVGs
      hvg <- try(Seurat::VariableFeatures(o), silent = TRUE)
      if (!is.character(hvg) || length(hvg) < 200) {
        o <- Seurat::FindVariableFeatures(o, selection.method = "vst",
                                          nfeatures = 2000, verbose = FALSE)
        hvg <- Seurat::VariableFeatures(o)
      }
      
      # 4) ScaleData (ensure assay has scale.data)
      if (!has_scaled(o, a)) {
        o <- Seurat::ScaleData(o, features = hvg, verbose = FALSE)
      }
      
      # 5) PCA
      if (!("pca" %in% Seurat::Reductions(o))) {
        npcs <- max(20, min(50, length(hvg)))
        o <- Seurat::RunPCA(o, features = hvg, npcs = npcs, verbose = FALSE)
      }
      
      # 6) SNN + UMAP (for visualization; prefer UMAP later)
      pcs_available <- ncol(Seurat::Embeddings(o, "pca"))
      dims_use <- 1:min(30, pcs_available)
      if (!"RNA_snn" %in% names(o@graphs)) {
        o <- Seurat::FindNeighbors(o, reduction = "pca", dims = dims_use, verbose = FALSE)
      }
      if (!("umap" %in% Seurat::Reductions(o))) {
        o <- Seurat::RunUMAP(o, reduction = "pca", dims = dims_use, verbose = FALSE)
      }
      o
    }
    pick_reduction <- function(o) {
      r <- Seurat::Reductions(o)
      if ("umap" %in% r) "umap" else if ("pca" %in% r) "pca" else stop("No UMAP/PCA found.")
    }
    make_hex_plot <- function(df, zcol, title) {
      if (!base::requireNamespace("hexbin", quietly = TRUE)) {
        # Fallback to 2D bins if hexbin missing
        return(
          ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y)) +
            ggplot2::stat_summary_2d(
              ggplot2::aes(z = .data[[zcol]], fill = ggplot2::after_stat(..value..)),
              fun = base::mean, bins = 120
            ) +
            ggplot2::scale_fill_viridis_c(name = title) +
            ggplot2::coord_equal() +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(
              panel.grid = ggplot2::element_blank(),
              axis.title = ggplot2::element_blank(),
              axis.text  = ggplot2::element_blank(),
              axis.ticks = ggplot2::element_blank()
            )
        )
      }
      ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y)) +
        ggplot2::stat_summary_hex(
          ggplot2::aes(z = .data[[zcol]], fill = ggplot2::after_stat(..value..)),
          fun = base::mean, bins = 80
        ) +
        ggplot2::scale_fill_viridis_c(name = title) +
        ggplot2::coord_equal() +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          axis.text  = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()
        )
    }
    make_violin <- function(o, col_unl, col_lab) {
      Seurat::Idents(o) <- "seurat_clusters"
      Seurat::VlnPlot(o,
                      features = c(col_unl, col_lab),
                      pt.size = 0.1, flip = TRUE
      )
    }
    
    # ---- reactive state ----
    rv <- shiny::reactiveValues(
      ready = FALSE,
      o = NULL,
      df = NULL,
      red = NULL,
      col_unl = NULL,
      col_lab = NULL
    )
    output$status <- shiny::renderText("Idle.")
    
    # (1) Run BSI analysis
    shiny::observeEvent(input$run, ignoreInit = TRUE, {
      output$status <- shiny::renderText("Analysis is running ...")
      
      shiny::withProgress(message = "Computing BSI", value = 0, {
        o <- get_obj(); req(o)
        if (!inherits(o, "Seurat")) stop("Reactive object is not a Seurat object.")
        
        # Ensure preprocessing EXACTLY as in the interactive run
        shiny::incProgress(0.25, detail = "Preprocessing (normalize → HVGs → scale → PCA → UMAP) ...")
        o <- ensure_preprocessing(o)
        
        # Choose embedding (prefer UMAP)
        red <- pick_reduction(o)
        emb <- Seurat::Embeddings(o, red)
        k <- 30
        
        # Unlabeled BSI: distance-gap between k-th and (k+1)-th neighbor
        shiny::incProgress(0.55, detail = "Unlabeled BSI (kNN distance gap) ...")
        k_ext <- k + 5
        nn <- RANN::nn2(data = emb, query = emb, k = k_ext)
        dst <- nn$nn.dists
        kth  <- dst[, k + 1]  # skip self
        k1th <- dst[, k + 2]
        bsi_raw <- base::pmax(k1th - kth, 0)
        rng <- stats::quantile(bsi_raw, probs = c(0.01, 0.99), na.rm = TRUE)
        bsi_unl <- (bsi_raw - rng[1]) / base::max(rng[2] - rng[1], .Machine$double.eps)
        bsi_unl <- base::pmin(base::pmax(bsi_unl, 0), 1)
        col_unl <- base::sprintf("BSI_%s_k%d", red, k)
        o[[col_unl]] <- bsi_unl
        
        # Label-aware BSI: fraction of neighbors with different cluster labels
        shiny::incProgress(0.8, detail = "Label-aware BSI (neighbor label disagreement) ...")
        if (!"seurat_clusters" %in% base::colnames(o@meta.data)) {
          # use same dims as PCA neighbors for consistency
          pcs_available <- ncol(Seurat::Embeddings(o, "pca"))
          dims_use <- 1:min(30, pcs_available)
          if (!"RNA_snn" %in% names(o@graphs)) {
            o <- Seurat::FindNeighbors(o, reduction = "pca", dims = dims_use, verbose = FALSE)
          }
          o <- Seurat::FindClusters(o, resolution = 0.6, verbose = FALSE)
        }
        nn2 <- RANN::nn2(data = emb, query = emb, k = k + 1) # includes self
        idx <- nn2$nn.idx
        labs <- base::as.integer(base::factor(o$seurat_clusters))
        labs_nn <- base::matrix(labs[idx], nrow = nrow(idx))[, -1, drop = FALSE]
        diff_frac <- base::rowMeans(labs_nn != labs)
        col_lab <- base::sprintf("BSI_labelaware_%s_k%d", red, k)
        o[[col_lab]] <- diff_frac
        
        # Cache plotting frame
        shiny::incProgress(0.95, detail = "Preparing plot data ...")
        df <- base::data.frame(
          x = emb[, 1],
          y = emb[, 2],
          BSI_unlabeled  = o[[col_unl]][, 1],
          BSI_labelaware = o[[col_lab]][, 1]
        )
        
        rv$o <- o
        rv$df <- df
        rv$red <- red
        rv$col_unl <- col_unl
        rv$col_lab <- col_lab
        rv$ready <- TRUE
      })
      
      output$status <- shiny::renderText("Analysis complete.")
    })
    
    # (2) Plot label-aware UMAP (on click only)
    shiny::observeEvent(input$plot_labelaware, ignoreInit = TRUE, {
      output$p_labelaware <- shiny::renderPlot({
        if (!isTRUE(rv$ready)) {
          return(ggplot2::ggplot() + ggplot2::theme_void() +
                   ggplot2::geom_text(ggplot2::aes(0, 0, label = "Run BSI analysis first")))
        }
        make_hex_plot(rv$df, "BSI_labelaware", "BSI (label-aware)")
      })
    })
    
    # Download: label-aware UMAP (PNG, 1000 dpi)
    output$dl_labelaware <- shiny::downloadHandler(
      filename = function() "BSI_labelaware_umap.png",
      content = function(file) {
        if (!isTRUE(rv$ready)) stop("Run BSI analysis first.")
        p <- make_hex_plot(rv$df, "BSI_labelaware", "BSI (label-aware)")
        ggplot2::ggsave(file, plot = p, width = 6, height = 5, units = "in", dpi = 1000)
      }
    )
    
    # (4) Plot unlabeled UMAP (on click only)
    shiny::observeEvent(input$plot_unlabeled, ignoreInit = TRUE, {
      output$p_unlabeled <- shiny::renderPlot({
        if (!isTRUE(rv$ready)) {
          return(ggplot2::ggplot() + ggplot2::theme_void() +
                   ggplot2::geom_text(ggplot2::aes(0, 0, label = "Run BSI analysis first")))
        }
        make_hex_plot(rv$df, "BSI_unlabeled", "BSI (unlabeled)")
      })
    })
    
    # Download: unlabeled UMAP (PNG, 1000 dpi)
    output$dl_unlabeled <- shiny::downloadHandler(
      filename = function() "BSI_unlabeled_umap.png",
      content = function(file) {
        if (!isTRUE(rv$ready)) stop("Run BSI analysis first.")
        p <- make_hex_plot(rv$df, "BSI_unlabeled", "BSI (unlabeled)")
        ggplot2::ggsave(file, plot = p, width = 6, height = 5, units = "in", dpi = 1000)
      }
    )
    
    # (5) Violin (on click only)
    shiny::observeEvent(input$plot_violin, ignoreInit = TRUE, {
      output$p_violin <- shiny::renderPlot({
        if (!isTRUE(rv$ready)) {
          return(ggplot2::ggplot() + ggplot2::theme_void() +
                   ggplot2::geom_text(ggplot2::aes(0, 0, label = "Run BSI analysis first")))
        }
        make_violin(rv$o, rv$col_unl, rv$col_lab)
      })
    })
    
    # Download: violin (PNG, 1000 dpi)
    output$dl_violin <- shiny::downloadHandler(
      filename = function() "BSI_violin.png",
      content = function(file) {
        if (!isTRUE(rv$ready)) stop("Run BSI analysis first.")
        p <- make_violin(rv$o, rv$col_unl, rv$col_lab)
        ggplot2::ggsave(file, plot = p, width = 7, height = 5, units = "in", dpi = 1000)
      }
    )
    
    # ---------- NEW: Download cluster-specific BSI CSV ----------
    output$dl_bsi_csv <- shiny::downloadHandler(
      filename = function() "BSI_clusters.csv",
      content = function(file) {
        if (!isTRUE(rv$ready)) stop("Run BSI analysis first.")
        o <- rv$o
        col_unl <- rv$col_unl
        col_lab <- rv$col_lab
        df <- base::data.frame(
          Cell_ID        = base::colnames(o),
          Cluster        = base::as.character(o$seurat_clusters),
          BSI_unlabeled  = base::as.numeric(o[[col_unl]][, 1]),
          BSI_labelaware = base::as.numeric(o[[col_lab]][, 1]),
          stringsAsFactors = FALSE
        )
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    # ------------------------------------------------------------
  })
}
