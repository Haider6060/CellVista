# server_cellentropy.R

cellentropy_server <- function(id, seurat_obj_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ---------------- State ----------------
    rv <- reactiveValues(
      sce = NULL,
      barPlot = NULL,
      umapPlot = NULL,
      heatPlot = NULL,
      geneEntropyAll = NULL,
      geneEntropyTop = NULL
    )
    
    # ---------------- Reset on dataset change ----------------
    observeEvent(seurat_obj_reactive(), {
      rv$sce <- NULL
      rv$barPlot <- NULL
      rv$umapPlot <- NULL
      rv$heatPlot <- NULL
      rv$geneEntropyAll <- NULL
      rv$geneEntropyTop <- NULL
      
      output$bar_plot <- renderPlot(NULL)
      output$umap_plot <- renderPlot(NULL)
      output$heatmap_plot <- renderPlot(NULL)
      output$done_msg <- renderText("")
    }, ignoreInit = FALSE)
    
    # ---------------- Helpers (no library() calls here) ----------------
    safe_coerce_dgc <- function(m) {
      if (inherits(m, "dgCMatrix")) return(m)
      if (is.matrix(m)) return(Matrix::Matrix(m, sparse = TRUE))
      tryCatch(as(m, "dgCMatrix"), error = function(...) Matrix::Matrix(as.matrix(m), sparse = TRUE))
    }
    safelayer <- function(assay, prefer = c("counts","raw","umi","data","lognorm","logcounts")) {
      av <- SeuratObject::Layers(assay); for (p in prefer) if (p %in% av) return(p); av[1]
    }
    as_sce_universal <- function(x) {
      if (inherits(x, "SingleCellExperiment")) return(x)
      if (inherits(x, "Seurat")) {
        da <- Seurat::DefaultAssay(x); assay_obj <- x[[da]]
        lyr_counts <- safelayer(assay_obj, c("counts","raw","umi"))
        lyr_data   <- safelayer(assay_obj, c("data","lognorm","logcounts","counts"))
        cnt <- SeuratObject::GetAssayData(assay_obj, layer = lyr_counts)
        dat <- SeuratObject::GetAssayData(assay_obj, layer = lyr_data)
        cell_order <- colnames(cnt)
        md <- x@meta.data[match(cell_order, rownames(x@meta.data)), , drop = FALSE]
        sce <- SingleCellExperiment::SingleCellExperiment(
          assays = list(counts = safe_coerce_dgc(cnt), logcounts = safe_coerce_dgc(dat))
        )
        SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(md, row.names = cell_order)
        return(sce)
      }
      stop("Unsupported object type for CellEntropy.")
    }
    choose_reference <- function(genes) {
      ref_h <- celldex::HumanPrimaryCellAtlasData()
      ref_m <- celldex::MouseRNAseqData()
      if (length(intersect(genes, rownames(ref_h))) >= length(intersect(genes, rownames(ref_m)))) ref_h else ref_m
    }
    
    # Sparse-safe per-cell entropy on a dgCMatrix (no dense ops)
    col_entropy_sparse <- function(P) {
      # P must contain probabilities (columns sum to 1); dgCMatrix (CSC)
      ent <- numeric(ncol(P))
      xr <- P@x
      pp <- P@p
      for (j in seq_len(ncol(P))) {
        if (pp[j] < pp[j+1]) {
          idx <- (pp[j] + 1L):pp[j+1L]
          x <- xr[idx]
          ent[j] <- -sum(x * log2(x))
        } else {
          ent[j] <- 0
        }
      }
      ent
    }
    
    # Sparse-safe per-gene entropy across cells for a subset (no dense row extraction)
    row_entropy_by_type <- function(sub_mat) {
      # sub_mat: dgCMatrix genes x cells
      totals <- Matrix::rowSums(sub_mat)                # fast, stays sparse
      keep   <- which(totals > 0)
      if (!length(keep)) return(rep(NA_real_, nrow(sub_mat)))
      
      # triplet form for per-row aggregation without densifying
      trip <- Matrix::summary(sub_mat[keep, , drop = FALSE])  # i (row in 'keep'), j, x
      # probability per (gene,row) across cells
      p <- trip$x / totals[keep][trip$i]
      # entropy contribution per nonzero
      contrib <- p * log2(p)
      # sum by row index in 'keep'
      s <- rowsum(contrib, group = trip$i, reorder = FALSE)
      out <- rep(NA_real_, nrow(sub_mat))
      out[keep] <- -as.numeric(s)
      out
    }
    
    # ---------------- Run Analysis (button) ----------------
    observeEvent(input$run_entropy, {
      obj <- seurat_obj_reactive()
      
      if (is.null(obj)) {
        showNotification("No Seurat object loaded. Please upload a .rds first.", type = "error", duration = 7)
        return(invisible(NULL))
      }
      
      withProgress(message = "CellEntropy: running analysis…", value = 0, {
        shinybusy::show_modal_spinner(text = "Computing…")
        
        tryCatch({
          # Convert to SCE
          incProgress(0.05, detail = "Converting to SingleCellExperiment")
          sce <- as_sce_universal(obj)
          
          # SingleR annotation
          incProgress(0.20, detail = "Annotating cell types (SingleR)")
          ref <- choose_reference(rownames(sce))
          common <- intersect(rownames(sce), rownames(ref))
          pred <- SingleR::SingleR(test = sce[common, ], ref = ref[common, ], labels = ref$label.main)
          pruned <- SingleR::pruneScores(pred)
          labels <- pred$labels; labels[is.na(pruned)] <- "Unknown"
          SummarizedExperiment::colData(sce)$cell_type <- labels
          
          # Per-cell entropy (sparse-safe)
          incProgress(0.40, detail = "Computing per-cell entropy")
          X <- SummarizedExperiment::assay(sce, "logcounts")
          X <- safe_coerce_dgc(X)
          X@x[X@x < 0] <- 0
          # normalize columns to probability
          cs <- Matrix::colSums(X); cs[cs == 0] <- 1
          P <- X
          P@x <- P@x / rep.int(cs, diff(P@p))
          ent_cell <- col_entropy_sparse(P)
          SummarizedExperiment::colData(sce)$entropy <- as.numeric(ent_cell)
          
          # Bar plot (single-color gradient: light-to-dark)
          incProgress(0.55, detail = "Preparing bar plot")
          df_bar <- data.frame(
            cell_type = SummarizedExperiment::colData(sce)$cell_type,
            entropy   = SummarizedExperiment::colData(sce)$entropy
          )
          df_bar <- dplyr::summarise(dplyr::group_by(df_bar, cell_type),
                                     mean_entropy = mean(entropy), .groups = "drop")
          
          rv$barPlot <- ggplot2::ggplot(df_bar, ggplot2::aes(x = stats::reorder(cell_type, mean_entropy),
                                                             y = mean_entropy, fill = mean_entropy)) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::geom_text(ggplot2::aes(label = round(mean_entropy, 2)), vjust = -0.5, size = 4) +
            ggplot2::scale_fill_gradient(low = "#dbe9f6", high = "#084594", guide = "none") +  # single hue
            ggplot2::theme_classic(base_size = 14) +
            ggplot2::labs(x = "Cell type", y = "Mean entropy")
          
          # UMAP + entropy landscape
          incProgress(0.70, detail = "Embedding (PCA/UMAP)")
          if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(sce)) {
            sce <- scater::runPCA(sce, ncomponents = 30, exprs_values = "logcounts")
            sce <- scater::runUMAP(sce, dimred = "PCA", n_neighbors = 30, min_dist = 0.3)
          }
          df_um <- data.frame(
            UMAP1 = SingleCellExperiment::reducedDim(sce, "UMAP")[,1],
            UMAP2 = SingleCellExperiment::reducedDim(sce, "UMAP")[,2],
            entropy = SummarizedExperiment::colData(sce)$entropy,
            cell_type = SummarizedExperiment::colData(sce)$cell_type
          )
          cents <- dplyr::summarise(dplyr::group_by(df_um, cell_type),
                                    UMAP1 = stats::median(UMAP1), UMAP2 = stats::median(UMAP2), .groups = "drop")
          rv$umapPlot <- ggplot2::ggplot(df_um, ggplot2::aes(UMAP1, UMAP2)) +
            ggplot2::stat_summary_2d(ggplot2::aes(z = entropy, fill = ggplot2::after_stat(value)), bins = 120) +
            viridis::scale_fill_viridis(name = "Entropy", option = "magma") +
            ggplot2::geom_point(data = cents, ggplot2::aes(UMAP1, UMAP2),
                                inherit.aes = FALSE, size = 1.2, color = "white") +
            ggrepel::geom_label_repel(data = cents,
                                      ggplot2::aes(UMAP1, UMAP2, label = cell_type),
                                      inherit.aes = FALSE, size = 3.0, label.size = 0.2,
                                      fill = "white", alpha = 0.95, max.overlaps = Inf,
                                      box.padding = 0.35, point.padding = 0.25, seed = 42) +
            ggplot2::coord_equal() +
            ggplot2::theme_void(base_size = 14) +
            ggplot2::theme(legend.position = "right")
          
          # Gene-level entropy per type (sparse-safe; no dense rows)
          incProgress(0.85, detail = "Computing gene-level entropy")
          mat2 <- SummarizedExperiment::assay(sce, "logcounts")
          mat2 <- safe_coerce_dgc(mat2); mat2@x[mat2@x < 0] <- 0
          ct <- factor(SummarizedExperiment::colData(sce)$cell_type)
          cts <- levels(ct)
          idx_by_type <- lapply(cts, function(k) which(ct == k))
          
          ent_list <- lapply(seq_along(cts), function(i) {
            cols <- idx_by_type[[i]]
            if (length(cols) < 2) return(NULL)
            sub <- mat2[, cols, drop = FALSE]                    # genes x cells (dgCMatrix)
            ent_rows <- row_entropy_by_type(sub)                 # vector length nrow(mat2)
            data.frame(cell_type = cts[i], gene = rownames(mat2),
                       gene_entropy = ent_rows, n_cells = length(cols),
                       stringsAsFactors = FALSE)
          })
          all_ge <- do.call(rbind, ent_list)
          all_ge <- all_ge[!is.na(all_ge$gene_entropy), ]
          rv$geneEntropyAll <- all_ge
          rv$geneEntropyTop <- dplyr::ungroup(
            dplyr::slice_max(dplyr::group_by(all_ge, cell_type),
                             order_by = gene_entropy, n = 50, with_ties = FALSE)
          )
          
          # Heatmap (top-30 per type; per-type mean expr; z-score by gene)
          incProgress(0.95, detail = "Building heatmap")
          top30 <- dplyr::ungroup(
            dplyr::slice_max(dplyr::group_by(all_ge, cell_type),
                             order_by = gene_entropy, n = 30, with_ties = FALSE)
          )
          genes_use <- unique(top30$gene)
          means_by_type <- lapply(cts, function(k) {
            cols <- which(ct == k); if (!length(cols)) return(NULL)
            data.frame(gene = genes_use,
                       cell_type = k,
                       mean_expr = as.numeric(Matrix::rowMeans(mat2[genes_use, cols, drop = FALSE])),
                       stringsAsFactors = FALSE)
          })
          hm_df <- do.call(rbind, means_by_type)
          zwide <- reshape2::dcast(hm_df, gene ~ cell_type, value.var = "mean_expr")
          mat_hm <- as.matrix(zwide[, -1, drop = FALSE]); rownames(mat_hm) <- zwide$gene
          mat_hm <- scale(mat_hm)  # z-score by gene
          hm_long <- reshape2::melt(mat_hm, varnames = c("gene","cell_type"), value.name = "z")
          
          rv$heatPlot <- ggplot2::ggplot(hm_long, ggplot2::aes(cell_type, gene, fill = z)) +
            ggplot2::geom_tile() +
            viridis::scale_fill_viridis(option = "mako", name = "Z-score") +
            ggplot2::theme_minimal(base_size = 12) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                           panel.grid = ggplot2::element_blank()) +
            ggplot2::labs(x = "Cell type", y = "Top-30 high-entropy genes")
          
          rv$sce <- sce
          output$done_msg <- renderText("Analysis completed.")
          showNotification("CellEntropy: analysis finished.", type = "message", duration = 6)
          
        }, error = function(e) {
          showNotification(paste("CellEntropy failed:", e$message), type = "error", duration = NULL)
        }, finally = {
          shinybusy::remove_modal_spinner()
        })
      })
    })
    
    # ---------------- Renderers ----------------
    observeEvent(input$plot_bar,   { req(rv$barPlot);    output$bar_plot    <- renderPlot(rv$barPlot, res = 120) })
    observeEvent(input$plot_umap,  { req(rv$umapPlot);   output$umap_plot   <- renderPlot(rv$umapPlot, res = 120) })
    observeEvent(input$plot_heatmap, {
      req(rv$heatPlot); output$heatmap_plot <- renderPlot(rv$heatPlot, res = 120)
    })
    
    # ---------------- Downloads (PNG & CSV prompts) ----------------
    output$download_bar <- downloadHandler(
      filename = function() "mean_entropy_barplot.png",
      content  = function(file) ggplot2::ggsave(file, plot = rv$barPlot, width = 8, height = 6, dpi = 1000)
    )
    output$download_umap <- downloadHandler(
      filename = function() "umap_entropy.png",
      content  = function(file) ggplot2::ggsave(file, plot = rv$umapPlot, width = 7, height = 5, dpi = 1000)
    )
    output$download_heatmap <- downloadHandler(
      filename = function() "top30_entropy_heatmap.png",
      content  = function(file) ggplot2::ggsave(file, plot = rv$heatPlot, width = 9, height = 12, dpi = 1000)
    )
    
    output$download_csv_all <- downloadHandler(
      filename = function() "gene_entropy_by_type_all.csv",
      content = function(file) {
        req(rv$geneEntropyAll)
        utils::write.csv(rv$geneEntropyAll, file, row.names = FALSE)
      }
    )
    output$download_csv_top <- downloadHandler(
      filename = function() "gene_entropy_top50_by_type.csv",
      content = function(file) {
        req(rv$geneEntropyTop)
        utils::write.csv(rv$geneEntropyTop, file, row.names = FALSE)
      }
    )
  })
}



















