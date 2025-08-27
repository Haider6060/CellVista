# server_qc.R — QC module server (robust QC metrics + clean PNG downloads)

qc_server <- function(id, seurat_obj_reactive) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(
      ready  = FALSE,
      obj    = NULL,   # QC-prepared Seurat object
      vln    = NULL,   # on-screen violin grid
      sct    = NULL,   # on-screen scatter
      thr    = NULL,   # list(umi_min, gene_min, mito_max)
      counts = NULL    # list(n_pass, n_fail)
    )
    
    # ---------- Helpers ----------
    .detect_mito_pattern <- function(obj) {
      feats <- rownames(obj)
      if (any(grepl("^MT-", feats))) return("^MT-")  # human
      if (any(grepl("^mt-", feats))) return("^mt-")  # mouse lower
      if (any(grepl("^Mt-", feats))) return("^Mt-")  # mixed
      "^MT-"  # fallback
    }
    .detect_ribo_pattern <- function(obj) {
      feats <- rownames(obj)
      if (any(grepl("^RPS", feats)) || any(grepl("^RPL", feats))) return("^RPS|^RPL") # human
      if (any(grepl("^Rps", feats)) || any(grepl("^Rpl", feats))) return("^Rps|^Rpl") # mouse
      "^RPS|^RPL"
    }
    .ensure_qc_cols <- function(obj) {
      meta_names <- colnames(obj@meta.data)
      if (!("percent.mt" %in% meta_names)) {
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = .detect_mito_pattern(obj))
      }
      if (!("percent.rb" %in% meta_names)) {
        obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = .detect_ribo_pattern(obj))
      }
      obj
    }
    
    # ---------- RUN QC ----------
    observeEvent(input$run_qc, {
      obj <- seurat_obj_reactive()
      validate(need(!is.null(obj), "Please upload a Seurat .rds first."))
      
      # Ensure QC columns
      obj <- .ensure_qc_cols(obj)
      
      # Thresholds (data-driven; change later if you add UI controls)
      df_qc <- FetchData(obj, vars = c("nCount_RNA","nFeature_RNA","percent.mt"))
      umi_min  <- as.numeric(quantile(df_qc$nCount_RNA,   0.05, na.rm = TRUE))
      gene_min <- as.numeric(quantile(df_qc$nFeature_RNA, 0.05, na.rm = TRUE))
      mito_max <- as.numeric(quantile(df_qc$percent.mt,   0.95, na.rm = TRUE))
      rv$thr <- list(umi_min = umi_min, gene_min = gene_min, mito_max = mito_max)
      
      # On-screen VIOLIN (multi-feature grid)
      n_groups <- if (!is.null(obj$orig.ident)) length(unique(obj$orig.ident)) else 1L
      cols <- if (n_groups <= 1) "#3c8dbc" else grDevices::hcl.colors(n_groups, palette = "Set3")
      
      rv$vln <- VlnPlot(
        obj,
        features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),
        group.by = if (!is.null(obj$orig.ident)) "orig.ident" else NULL,
        pt.size = 0,
        cols = cols,
        ncol = 2
      ) + theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        plot.title  = element_text(size = 14, face = "bold"),
        axis.title  = element_text(size = 12)
      )
      
      # On-screen SCATTER (UMIs vs % mito) built with ggplot (no FeatureScatter dependency)
      df_qc$QC_status <- with(df_qc, ifelse(
        nCount_RNA < umi_min | nFeature_RNA < gene_min | percent.mt > mito_max,
        "FAIL (low UMIs / low genes / high mito)",
        "PASS (meets all QC thresholds)"
      ))
      df_qc$QC_status <- factor(
        df_qc$QC_status,
        levels = c("FAIL (low UMIs / low genes / high mito)",
                   "PASS (meets all QC thresholds)")
      )
      n_pass <- sum(df_qc$QC_status == "PASS (meets all QC thresholds)", na.rm = TRUE)
      n_fail <- sum(df_qc$QC_status == "FAIL (low UMIs / low genes / high mito)", na.rm = TRUE)
      rv$counts <- list(n_pass = n_pass, n_fail = n_fail)
      
      rv$sct <- ggplot2::ggplot(df_qc, ggplot2::aes(x = nCount_RNA, y = percent.mt, color = QC_status)) +
        ggplot2::annotate("rect", xmin = -Inf, xmax = umi_min, ymin = -Inf, ymax = Inf,
                          alpha = 0.08, fill = "red") +
        ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = mito_max, ymax = Inf,
                          alpha = 0.08, fill = "red") +
        ggplot2::geom_point(alpha = 0.55, size = 0.8) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_color_manual(values = c(
          "FAIL (low UMIs / low genes / high mito)" = "#d62728",
          "PASS (meets all QC thresholds)"          = "#2ca02c"
        ), name = "Cell QC classification") +
        ggplot2::geom_vline(xintercept = umi_min,  linetype = "dashed") +
        ggplot2::geom_hline(yintercept = mito_max, linetype = "dashed") +
        ggplot2::labs(
          title = "QC scatter: UMIs vs % mitochondrial",
          subtitle = sprintf("PASS = UMIs ≥ %.0f, genes ≥ %.0f, mito ≤ %.1f%% | PASS: %d  FAIL: %d",
                             umi_min, gene_min, mito_max, n_pass, n_fail),
          x = "UMIs per cell (log10)",
          y = "% mitochondrial reads"
        ) +
        ggplot2::annotate("text", x = umi_min, y = -Inf, vjust = -0.7, hjust = -0.05,
                          label = sprintf("UMI threshold = %.0f", umi_min), size = 3.6) +
        ggplot2::annotate("text", x = min(df_qc$nCount_RNA, na.rm = TRUE), y = mito_max,
                          vjust = -0.7, hjust = 0,
                          label = sprintf("Mito threshold = %.1f%%", mito_max), size = 3.6) +
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = "top")
      
      rv$obj <- obj
      rv$ready <- TRUE
      showNotification("QC completed successfully.", type = "message")
    }, priority = 10)
    
    # ---------- STATUS ----------
    output$qc_status <- renderUI({
      if (!rv$ready) return(tags$span("Ready: click \"Run QC analysis\"."))
      thr <- rv$thr; cnt <- rv$counts %||% list(n_pass = NA, n_fail = NA)
      tags$span(
        sprintf("QC ready. Thresholds — UMIs ≥ %.0f, genes ≥ %.0f, mito ≤ %.1f%% | PASS: %s  FAIL: %s",
                thr$umi_min, thr$gene_min, thr$mito_max, cnt$n_pass, cnt$n_fail),
        style = "color:#2ca02c; font-weight:600;"
      )
    })
    
    # ---------- SHOW PLOTS ON CLICK ----------
    observeEvent(input$show_vln, { req(rv$ready, rv$vln); output$vln_plot <- renderPlot(rv$vln, res = 120) })
    observeEvent(input$show_sct, { req(rv$ready, rv$sct); output$sct_plot <- renderPlot(rv$sct, res = 120) })
    
    # ---------- DOWNLOADS (true PNG @ 1000 dpi) ----------
    # Violin: build 4 separate panels, stack vertically, very wide canvas; use short labels
    output$dl_vln <- downloadHandler(
      filename    = function() sprintf("cellprint_qc_violin_%s.png", format(Sys.time(), "%Y%m%d_%H%M%S")),
      contentType = "image/png",
      content     = function(file) {
        req(rv$ready, rv$obj)
        obj <- rv$obj
        
        # Short, unique labels for samples; fallback to single group if orig.ident missing
        if (is.null(obj$orig.ident)) {
          obj$`.__cp_short` <- factor(rep("All", ncol(obj)))
        } else {
          grp_fac    <- factor(obj$orig.ident)
          long_lvls  <- levels(grp_fac)
          short_lvls <- abbreviate(long_lvls, minlength = 8, strict = TRUE)
          obj$`.__cp_short` <- factor(grp_fac, labels = short_lvls)
        }
        
        feats <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb")
        n_groups <- length(levels(obj$`.__cp_short`))
        cols <- if (n_groups <= 1) "#3c8dbc" else grDevices::hcl.colors(n_groups, palette = "Set3")
        
        vln_list <- lapply(feats, function(f) {
          VlnPlot(
            obj, features = f, group.by = ".__cp_short",
            pt.size = 0, cols = cols
          ) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6, lineheight = 0.8),
              plot.title  = ggplot2::element_text(size = 14, face = "bold"),
              axis.title  = ggplot2::element_text(size = 12),
              plot.margin = ggplot2::margin(t = 6, r = 6, b = 26, l = 6)
            ) +
            ggplot2::labs(title = f, x = "Samples (abbreviated)")
        })
        
        # Wide canvas per group
        width_in  <- max(28, min(60, 0.75 * n_groups))
        height_in <- 3.6 * length(vln_list)
        
        # Save as real PNG with base device (prevents HTML fallbacks)
        png(filename = file, width = width_in, height = height_in, units = "in", res = 1000)
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = length(vln_list), ncol = 1)))
        for (i in seq_along(vln_list)) {
          print(vln_list[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
        }
        dev.off()
      }
    )
    
    # Scatter: save ggplot directly as PNG
    output$dl_sct <- downloadHandler(
      filename    = function() sprintf("cellprint_qc_scatter_%s.png", format(Sys.time(), "%Y%m%d_%H%M%S")),
      contentType = "image/png",
      content     = function(file) {
        req(rv$ready, rv$sct)
        ggplot2::ggsave(
          filename  = file,
          plot      = rv$sct,
          device    = "png",
          dpi       = 1000,
          width     = 12,
          height    = 8,
          units     = "in",
          limitsize = FALSE
        )
      }
    )
  })
}
