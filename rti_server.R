# rti_server.R — Cellular Regulatory Turbulence Index (RTI) module server
# NOTE: No library() calls here. All packages are loaded from global.R.

rti_server <- function(id, seurat_obj_reactive, npcs = 20) {
  moduleServer(id, function(input, output, session) {
    
    rv <- reactiveValues(obj = NULL, summary = NULL, p_umap = NULL, p_vln = NULL)
    
    # ---------- Small helpers (self-contained) ----------
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    
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
    
    has_scaled <- function(obj, assay) {
      # Seurat v5 prefers 'layer'; fall back to 'slot' if needed
      ok <- FALSE
      try({
        m <- GetAssayData(obj, assay = assay, layer = "scale.data")
        ok <- is.matrix(m) || inherits(m, "dgCMatrix")
        if (ok) ok <- (nrow(m) > 0 && ncol(m) > 0)
      }, silent = TRUE)
      if (!ok) {
        try({
          m <- GetAssayData(obj, assay = assay, slot = "scale.data")
          ok <- (is.matrix(m) || inherits(m, "dgCMatrix")) && (nrow(m) > 0 && ncol(m) > 0)
        }, silent = TRUE)
      }
      ok
    }
    
    build_neighbors_lists <- function(g) {
      A <- as(g, "dgCMatrix")
      idx_list <- vector("list", ncol(A))
      wts_list <- vector("list", ncol(A))
      for (i in seq_len(ncol(A))) {
        nz <- which(A[i, ] != 0)
        nbrs <- c(nz, i)                           # neighbors + self
        wts  <- c(as.numeric(A[i, nz]), 1e-6)      # tiny self-weight
        if (length(unique(nbrs)) != length(nbrs)) {
          df <- data.frame(nbrs = nbrs, wts = wts)
          df <- stats::aggregate(wts ~ nbrs, data = df, FUN = sum)
          nbrs <- as.integer(df$nbrs)
          wts  <- as.numeric(df$wts)
        }
        idx_list[[i]] <- as.integer(nbrs)
        wts_list[[i]] <- as.numeric(wts)
      }
      list(idx = idx_list, wts = wts_list)
    }
    
    graph_pseudotime <- function(g, cellnames, nCount) {
      # Work on sparse graph; build igraph edges without densifying
      A <- as(g, "dgCMatrix")
      A_sym <- as(A + t(A), "dgCMatrix")
      S <- summary(A_sym)  # i, j, x (1-based)
      df_edges <- data.frame(
        from   = cellnames[S$i],
        to     = cellnames[S$j],
        weight = 1 / S$x,   # stronger edge -> smaller distance
        stringsAsFactors = FALSE
      )
      g_ig <- igraph::graph_from_data_frame(
        df_edges, directed = FALSE,
        vertices = data.frame(name = cellnames, stringsAsFactors = FALSE)
      )
      root_idx  <- which.min(nCount)
      root_name <- cellnames[root_idx]
      d <- as.numeric(igraph::distances(g_ig, v = root_name, to = igraph::V(g_ig),
                                        weights = igraph::E(g_ig)$weight)[1, ])
      d[!is.finite(d)] <- max(d[is.finite(d)], na.rm = TRUE)
      (d - min(d)) / (max(d) - min(d) + 1e-6)
    }
    
    ensure_umap <- function(obj, use_dims) {
      if (!"umap" %in% names(obj@reductions)) {
        obj <- RunUMAP(obj, dims = use_dims, verbose = FALSE)
      }
      obj
    }
    
    velocity_from_neighbors <- function(X, tvec, idx_list, wts_list) {
      n <- nrow(X); p <- ncol(X)
      vel <- matrix(NA_real_, nrow = n, ncol = p,
                    dimnames = list(rownames(X), colnames(X)))
      for (i in seq_len(n)) {
        nbrs <- idx_list[[i]]
        w    <- wts_list[[i]]
        tt   <- tvec[nbrs]
        W    <- sum(w)
        tbar <- sum(w * tt) / W
        tcen <- tt - tbar
        denom <- sum(w * tcen * tcen) + 1e-8
        Xi <- X[nbrs, , drop = FALSE]
        num <- as.numeric(crossprod(tcen * w, Xi))
        vel[i, ] <- num / denom
      }
      vel
    }
    
    rti_components <- function(vel, idx_list, tvec) {
      eps <- 1e-8
      # Local mean velocity
      mu_mat <- do.call(rbind, lapply(idx_list, function(ix) colMeans(vel[ix, , drop = FALSE])))
      rownames(mu_mat) <- rownames(vel)
      mu_norm  <- sqrt(rowSums(mu_mat^2))
      mu_norm2 <- mu_norm * mu_norm
      v_sq     <- rowSums(vel^2)
      
      # Turbulence Intensity (TI)
      n_i     <- vapply(idx_list, length, integer(1))
      sum_vsq <- vapply(idx_list, function(ix) sum(v_sq[ix]), numeric(1))
      k_i     <- (sum_vsq - n_i * mu_norm2) / (2 * n_i + eps)
      TI      <- sqrt(pmax(0, 2 * k_i)) / (mu_norm + eps)
      
      # Directional Disorder (DD)
      DD <- numeric(nrow(vel))
      for (i in seq_len(nrow(vel))) {
        ix  <- idx_list[[i]]
        Vi  <- vel[ix, , drop = FALSE]
        mui <- t(mu_mat[i, , drop = FALSE])
        dots  <- as.numeric(Vi %*% mui)
        vnorm <- sqrt(rowSums(Vi^2))
        DD[i] <- 1 - mean(dots / (vnorm * (mu_norm[i] + eps) + eps), na.rm = TRUE)
      }
      
      # Intermittency (IM)
      IM <- numeric(nrow(vel))
      for (i in seq_len(nrow(vel))) {
        ix  <- idx_list[[i]]
        ufl <- sweep(vel[ix, , drop = FALSE], 2, mu_mat[i, ], FUN = "-")
        s   <- sqrt(rowSums(ufl^2))
        med <- stats::median(s)
        IM[i] <- stats::mad(s, constant = 1) / (med + eps)
      }
      
      # Spectral Entropy (H)
      H <- numeric(nrow(vel))
      for (i in seq_len(nrow(vel))) {
        ix  <- idx_list[[i]]
        ord <- ix[order(tvec[ix])]
        s   <- sqrt(rowSums(vel[ord, , drop = FALSE]^2))
        m   <- length(s)
        if (m < 4) { H[i] <- 0; next }
        s <- s - mean(s)
        S <- stats::fft(s)
        P <- Mod(S)^2 / m
        K <- floor(m / 2)
        if (K < 2) { H[i] <- 0; next }
        p <- P[2:(K+1)]
        ps <- sum(p)
        if (ps <= 0) { H[i] <- 0; next }
        p <- p / ps
        H[i] <- -sum(p * log(p)) / log(length(p))
      }
      
      list(TI = TI, DD = DD, IM = IM, H = H)
    }
    
    rti_aggregate <- function(TI, DD, IM, H) {
      eps <- 1e-8
      components <- cbind(TI = TI, DD = DD, IM = IM, H = H)
      robust_z <- function(x) (x - stats::median(x)) / (stats::mad(x, constant = 1) + eps)
      Z <- apply(components, 2, robust_z)
      raw <- rowMeans(Z)
      (raw - min(raw)) / (max(raw) - min(raw) + eps)
    }
    
    make_cluster_summary <- function(meta) {
      if (!("seurat_clusters" %in% colnames(meta)) || !("RTI" %in% colnames(meta))) {
        return(data.frame(message = "Clusters or RTI not available.", stringsAsFactors = FALSE))
      }
      cl <- meta$seurat_clusters
      r  <- meta$RTI
      u  <- sort(unique(cl))
      n  <- integer(length(u))
      med <- q25 <- q75 <- numeric(length(u))
      for (i in seq_along(u)) {
        ix <- which(cl == u[i])
        vals <- r[ix]
        n[i]   <- length(ix)
        med[i] <- stats::median(vals, na.rm = TRUE)
        q25[i] <- stats::quantile(vals, 0.25, na.rm = TRUE, names = FALSE)
        q75[i] <- stats::quantile(vals, 0.75, na.rm = TRUE, names = FALSE)
      }
      df <- data.frame(
        cluster = u,
        n_cells = n,
        RTI_median = med,
        RTI_q25 = q25,
        RTI_q75 = q75,
        stringsAsFactors = FALSE
      )
      df[order(-df$RTI_median), , drop = FALSE]
    }
    
    # ---------- Run RTI on click (with progress + spinner) ----------
    observeEvent(input$run, {
      obj <- seurat_obj_reactive()
      validate(need(!is.null(obj) && inherits(obj, "Seurat"),
                    "Please upload a valid Seurat object in the Data tab first."))
      
      use_npcs <- max(10, as.numeric(input$npcs %||% npcs))
      use_dims <- 1:use_npcs
      
      shinybusy::show_modal_spinner(text = "Running RTI analysis…", color = "#3c8dbc", spin = "fading-circle")
      on.exit(shinybusy::remove_modal_spinner(), add = TRUE)
      
      withProgress(message = "RTI analysis", value = 0, {
        incProgress(0.02, detail = "Preparing assay / QC covariates")
        obj <- set_active_assay(obj)
        obj <- add_percent_mt(obj)
        assay_now <- DefaultAssay(obj)
        
        incProgress(0.12, detail = "Normalization & variable features")
        if (length(VariableFeatures(obj)) == 0) {
          obj <- NormalizeData(obj, verbose = FALSE)
          obj <- FindVariableFeatures(
            obj, selection.method = "vst",
            nfeatures = max(200, min(2000, nrow(obj))),
            verbose = FALSE
          )
        }
        
        incProgress(0.28, detail = "Scaling (regressing percent.mt)")
        if (!has_scaled(obj, assay_now)) {
          vars <- intersect("percent.mt", colnames(obj@meta.data))
          obj <- ScaleData(
            obj,
            features = VariableFeatures(obj),
            vars.to.regress = vars,
            verbose = FALSE
          )
        }
        
        incProgress(0.40, detail = "PCA")
        if (!"pca" %in% names(obj@reductions)) {
          obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = use_npcs, verbose = FALSE)
        }
        
        incProgress(0.52, detail = "Neighbors / graph")
        g <- obj@graphs$RNA_snn
        if (is.null(g)) {
          obj <- FindNeighbors(obj, dims = use_dims, verbose = FALSE)
          g <- obj@graphs$RNA_snn %||% obj@graphs[[1]]
        }
        
        incProgress(0.58, detail = "UMAP (for visualization)")
        obj <- ensure_umap(obj, use_dims)
        
        incProgress(0.66, detail = "Neighbor lists")
        nb <- build_neighbors_lists(g)
        idx_list <- nb$idx; wts_list <- nb$wts
        
        incProgress(0.74, detail = "Graph pseudotime")
        tvec <- graph_pseudotime(g, colnames(obj), obj$nCount_RNA)
        obj$pseudotime <- tvec
        
        incProgress(0.82, detail = "Local velocity")
        X <- Embeddings(obj, "pca")[, use_dims, drop = FALSE]
        vel <- velocity_from_neighbors(X, tvec, idx_list, wts_list)
        obj@reductions$vel <- CreateDimReducObject(embeddings = vel, key = "VEL_", assay = DefaultAssay(obj))
        
        incProgress(0.90, detail = "RTI components")
        comps <- rti_components(vel, idx_list, tvec)
        obj$RTI_TI <- comps$TI
        obj$RTI_DD <- comps$DD
        obj$RTI_IM <- comps$IM
        obj$RTI_H  <- comps$H
        
        incProgress(0.96, detail = "Aggregate RTI & plots")
        obj$RTI <- rti_aggregate(obj$RTI_TI, obj$RTI_DD, obj$RTI_IM, obj$RTI_H)
        
        # If clusters are absent, create them (so violin/summary work)
        if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
          # use existing neighbor graph/dims; resolution = 0.8 like your outside analysis
          obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
        }
        
        rv$p_umap <- FeaturePlot(obj, features = "RTI", reduction = "umap", cols = c("navy", "gold")) +
          ggplot2::ggtitle("RTI per cell") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14))
        
        rv$p_vln <- VlnPlot(obj, features = "RTI", group.by = "seurat_clusters", pt.size = 0) +
          ggplot2::ggtitle("RTI by cluster") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14))
        
        rv$summary <- make_cluster_summary(obj@meta.data)
        
        rv$obj <- obj
        incProgress(1, detail = "Done")
      })
    })
    
    # ---------- Outputs ----------
    output$umap_rti <- renderPlot({ req(rv$p_umap); rv$p_umap }, res = 120)
    output$vln_rti  <- renderPlot({ req(rv$p_vln);  rv$p_vln  }, res = 120)
    output$tbl_summary <- renderTable({ req(rv$summary); rv$summary }, striped = TRUE, bordered = TRUE, digits = 3)
    
    # ---------- Downloads ----------
    output$dl_umap_rti <- downloadHandler(
      filename = function() paste0("RTI_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) ggplot2::ggsave(file, plot = rv$p_umap, width = 6, height = 5, dpi = 1000)
    )
    output$dl_vln_rti <- downloadHandler(
      filename = function() paste0("RTI_violin_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) ggplot2::ggsave(file, plot = rv$p_vln, width = 6, height = 5, dpi = 1000)
    )
    output$dl_summary_csv <- downloadHandler(
      filename = function() paste0("RTI_cluster_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) utils::write.csv(rv$summary, file = file, row.names = FALSE)
    )
    
    # =========================
    # NEW: Early-warning add-on
    # =========================
    # (Relies on RTI having been computed above; uses rv$obj)
    
    # small helper (plot smoothing only)
    ew_moving_average <- function(x, k = 5L) {
      out <- as.numeric(stats::filter(x, rep(1/k, k), sides = 2))
      out[is.na(out)] <- x[is.na(out)]
      out
    }
    
    observeEvent(input$run_ew, {
      validate(need(!is.null(rv$obj), "Run RTI first (the early-warning uses its pseudotime/graph)."))
      obj <- rv$obj
      
      # parameters (must exist in UI; if not present, defaults are used harmlessly)
      nbin  <- max(20, as.integer(input$ew_bins %||% 50))
      tau   <- as.numeric(input$ew_tau %||% 0.20)
      win   <- max(3, as.integer(input$ew_win %||% 5))
      min_n <- max(1, as.integer(input$ew_min_n %||% 50))
      
      shinybusy::show_modal_spinner(text = "Computing early-warning…", color = "#3c8dbc", spin = "fading-circle")
      on.exit(shinybusy::remove_modal_spinner(), add = TRUE)
      
      withProgress(message = "Early-warning", value = 0, {
        incProgress(0.1, detail = "Neighbor label-entropy")
        # neighbors (SNN already built above)
        g <- obj@graphs$RNA_snn %||% obj@graphs[[1]]
        validate(need(!is.null(g), "Neighbor graph not found. Please run clustering / neighbors first."))
        nb <- build_neighbors_lists(g)
        idx_list <- nb$idx
        
        validate(need("seurat_clusters" %in% colnames(obj@meta.data), "seurat_clusters missing."))
        lab <- as.character(obj$seurat_clusters)
        H_label <- numeric(length(idx_list))
        for (i in seq_along(idx_list)) {
          labs <- lab[idx_list[[i]]]
          tt   <- table(labs)
          p    <- as.numeric(tt) / sum(tt)
          K    <- length(p)
          if (K <= 1) { H_label[i] <- 0 } else {
            H_label[i] <- -sum(p * log(p)) / log(K)
          }
        }
        obj$H_label <- H_label
        
        incProgress(0.35, detail = "Binning pseudotime")
        validate(need("pseudotime" %in% colnames(obj@meta.data), "pseudotime missing."))
        pt  <- obj$pseudotime
        rti <- obj$RTI
        brks    <- seq(0, 1, length.out = nbin + 1)
        centers <- head(brks, -1) + diff(brks)/2
        bin_id  <- cut(pt, breaks = brks, include.lowest = TRUE, labels = FALSE)
        counts  <- as.numeric(table(factor(bin_id, levels = 1:nbin)))
        
        prof <- stats::aggregate(
          data.frame(H = H_label, RTI = rti),
          by = list(bin = bin_id),
          FUN = function(x) if (length(x)) stats::median(x, na.rm = TRUE) else NA_real_
        )
        prof$pseudotime <- centers[prof$bin]
        prof$n <- counts[prof$bin]
        prof <- prof[is.finite(prof$H) & is.finite(prof$RTI), ]
        prof <- prof[order(prof$pseudotime), ]
        
        incProgress(0.6, detail = "Commitment / peak / lead")
        low  <- prof$H <= tau
        ok_n <- prof$n >= min_n
        persist <- logical(nrow(prof))
        for (i in seq_len(nrow(prof))) {
          j <- i:min(nrow(prof), i + win - 1)
          persist[i] <- all(low[j] & ok_n[j], na.rm = TRUE)
        }
        commit_idx <- which(persist)[1]
        if (length(commit_idx) == 0 || is.na(commit_idx)) {
          showNotification("No commitment point detected with current thresholds.", type = "warning")
          rv$ew_plot <- NULL
          rv$ew_summary_df <- data.frame(
            commitment_pt = NA_real_, peak_pt = NA_real_, lead_time = NA_real_,
            peak_RTI = NA_real_, p_lead = NA_real_, p_peak = NA_real_,
            tau = tau, persist_bins = win, min_cells_per_bin = min_n,
            n_bins_used = nrow(prof)
          )
          rv$ew_profile <- prof
          return(invisible(NULL))
        }
        commit_pt <- prof$pseudotime[commit_idx]
        pre_idx <- which(prof$pseudotime < commit_pt)
        peak_idx <- pre_idx[ which.max(prof$RTI[pre_idx]) ]
        peak_pt  <- prof$pseudotime[peak_idx]
        peak_rti <- prof$RTI[peak_idx]
        lead_time <- commit_pt - peak_pt
        
        incProgress(0.8, detail = "Permutation p-values")
        set.seed(123)
        B <- 200
        lead_perm <- numeric(B); peak_perm <- numeric(B)
        bin_profile <- function(pt_vec, rti_vec, nbin) {
          br <- seq(0, 1, length.out = nbin + 1)
          id <- cut(pt_vec, breaks = br, include.lowest = TRUE, labels = FALSE)
          ctr <- head(br, -1) + diff(br)/2
          ct  <- as.numeric(table(factor(id, levels = 1:nbin)))
          pr <- stats::aggregate(
            data.frame(RTI = rti_vec),
            by = list(bin = id),
            FUN = function(x) if (length(x)) stats::median(x, na.rm = TRUE) else NA_real_
          )
          pr$pseudotime <- ctr[pr$bin]
          pr$n <- ct[pr$bin]
          pr <- pr[is.finite(pr$RTI), ]
          pr[order(pr$pseudotime), ]
        }
        for (b in seq_len(B)) {
          pt_perm <- sample(pt)
          pr <- bin_profile(pt_perm, rti, nbin)
          pre <- which(pr$pseudotime < commit_pt)
          if (length(pre) == 0) { lead_perm[b] <- NA_real_; peak_perm[b] <- NA_real_; next }
          k <- pre[ which.max(pr$RTI[pre]) ]
          peak_perm[b] <- pr$RTI[k]
          lead_perm[b] <- commit_pt - pr$pseudotime[k]
        }
        lead_perm <- lead_perm[is.finite(lead_perm)]
        peak_perm <- peak_perm[is.finite(peak_perm)]
        p_lead <- (1 + sum(lead_perm >= lead_time)) / (length(lead_perm) + 1)
        p_peak <- (1 + sum(peak_perm >= peak_rti)) / (length(peak_perm) + 1)
        
        incProgress(0.95, detail = "Publication plot")
        plotdf <- data.frame(pseudotime = prof$pseudotime, RTI = prof$RTI, H = prof$H)
        plotdf$RTI_s <- ew_moving_average(plotdf$RTI, 5)
        plotdf$H_s   <- ew_moving_average(plotdf$H,   5)
        longdf <- rbind(
          data.frame(pseudotime = plotdf$pseudotime, value = plotdf$RTI_s, series = "RTI"),
          data.frame(pseudotime = plotdf$pseudotime, value = plotdf$H_s,   series = "Neighbor label-entropy")
        )
        peak_df <- data.frame(pseudotime = peak_pt, value = peak_rti)
        cols <- c("RTI" = "#2b8cbe", "Neighbor label-entropy" = "#d95f0e")
        
        rv$ew_plot <- ggplot2::ggplot(longdf, ggplot2::aes(pseudotime, value, color = series, linetype = series)) +
          ggplot2::geom_ribbon(
            data = subset(longdf, series == "RTI" & pseudotime <= commit_pt),
            ggplot2::aes(x = pseudotime, ymin = 0, ymax = value),
            inherit.aes = FALSE, fill = cols["RTI"], alpha = 0.12
          ) +
          ggplot2::geom_line(linewidth = 1.1) +
          ggplot2::geom_vline(xintercept = commit_pt, color = "#333333", linewidth = 1) +
          ggplot2::geom_point(data = peak_df, ggplot2::aes(x = pseudotime, y = value),
                              color = cols["RTI"], size = 2.6, inherit.aes = FALSE) +
          ggplot2::annotate("segment", x = peak_pt, xend = commit_pt, y = peak_rti, yend = peak_rti,
                            colour = cols["RTI"], arrow = grid::arrow(length = grid::unit(0.18, "cm"))) +
          ggplot2::annotate("text", x = (peak_pt + commit_pt)/2, y = peak_rti + 0.035,
                            label = sprintf("lead = %.2f (p = %.3f)", lead_time, p_lead),
                            colour = cols["RTI"]) +
          ggplot2::scale_color_manual(values = cols, name = "") +
          ggplot2::scale_linetype_manual(values = c("RTI" = "solid", "Neighbor label-entropy" = "longdash"), name = "") +
          ggplot2::labs(
            title = "RTI early-warning",
            subtitle = sprintf("Commitment = %.2f, Peak RTI = %.2f (p_peak = %.3f)", commit_pt, peak_rti, p_peak),
            x = "Pseudotime (binned)", y = "Normalized value (0–1)"
          ) +
          ggplot2::theme_classic(base_size = 14) +
          ggplot2::theme(
            legend.position = "top",
            legend.key.width = grid::unit(2, "lines"),
            plot.title = ggplot2::element_text(face = "bold", size = 16),
            plot.subtitle = ggplot2::element_text(size = 12),
            panel.border = ggplot2::element_rect(color = "#333333", fill = NA, linewidth = 0.5)
          )
        
        rv$ew_profile <- prof
        rv$ew_summary_df <- data.frame(
          commitment_pt = commit_pt,
          peak_pt = peak_pt,
          lead_time = lead_time,
          peak_RTI = peak_rti,
          p_lead = p_lead,
          p_peak = p_peak,
          tau = tau,
          persist_bins = win,
          min_cells_per_bin = min_n,
          n_bins_used = nrow(prof),
          stringsAsFactors = FALSE
        )
        
        incProgress(1, detail = "Done")
      })
    })
    
    output$ew_plot <- renderPlot({ req(rv$ew_plot); rv$ew_plot }, res = 120)
    output$ew_summary <- renderTable({ req(rv$ew_summary_df); rv$ew_summary_df }, bordered = TRUE, striped = TRUE, digits = 3)
    
    output$dl_ew_plot <- downloadHandler(
      filename = function() paste0("RTI_early_warning_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) ggplot2::ggsave(file, plot = rv$ew_plot, width = 6, height = 5, dpi = 1000)
    )
    output$dl_ew_csv <- downloadHandler(
      filename = function() paste0("RTI_early_warning_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) utils::write.csv(rv$ew_summary_df, file = file, row.names = FALSE)
    )
    
    # =========================================
    # NEW: Cell-type annotation (SingleR + celldex)
    # =========================================
    # helpers for SCE conversion + reference choice + major-type mapping
    safe_coerce_dgc <- function(m) {
      if (inherits(m, "dgCMatrix")) return(m)
      if (is.matrix(m)) return(Matrix::Matrix(m, sparse = TRUE))
      tryCatch(as(m, "dgCMatrix"), error = function(...) Matrix::Matrix(as.matrix(m), sparse = TRUE))
    }
    safelayer <- function(assay, prefer = c("counts","raw","umi","data","lognorm","logcounts")) {
      av <- SeuratObject::Layers(assay)
      for (p in prefer) if (p %in% av) return(p)
      av[1]
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
      stop("Unsupported object type for SingleR conversion.")
    }
    choose_reference <- function(genes) {
      ref_h <- celldex::HumanPrimaryCellAtlasData()
      ref_m <- celldex::MouseRNAseqData()
      if (length(intersect(genes, rownames(ref_h))) >= length(intersect(genes, rownames(ref_m)))) ref_h else ref_m
    }
    map_major_type <- function(labels) {
      x <- tolower(gsub("[^A-Za-z0-9_ ]", " ", labels))
      major <- rep("Other", length(x))
      set <- function(idx, val) { major[idx] <<- val }
      set(grepl("t[_ ]?cell|cd4|cd8|th", x), "T_cells")
      set(grepl("\\bb[_ ]?cell|plasma", x), "B_cells")
      set(grepl("\\bnk\\b|natural[_ ]killer", x), "NK_cells")
      set(grepl("dendrit|\\bdc\\b", x), "Dendritic")
      set(grepl("mono|macro|neutro|granulo|myeloid|baso|eos", x), "Myeloid")
      set(grepl("endothel", x), "Endothelial")
      set(grepl("fibro|stromal", x), "Fibroblasts")
      set(grepl("epithel", x), "Epithelial")
      set(grepl("eryth|reticulocyte", x), "Erythroid")
      set(grepl("platelet|megakaryo", x), "Platelets")
      set(grepl("smooth[_ ]?muscle", x), "Smooth_muscle")
      major
    }
    
    observeEvent(input$run_types, {
      # Need a Seurat object (can run before or after RTI)
      base_obj <- rv$obj %||% seurat_obj_reactive()
      validate(need(!is.null(base_obj) && inherits(base_obj, "Seurat"),
                    "Please upload a valid Seurat object in the Data tab first."))
      
      shinybusy::show_modal_spinner(text = "Annotating cell types (SingleR)…", color = "#3c8dbc", spin = "fading-circle")
      on.exit(shinybusy::remove_modal_spinner(), add = TRUE)
      
      withProgress(message = "Cell-type annotation", value = 0, {
        incProgress(0.10, detail = "Converting to SingleCellExperiment")
        sce <- as_sce_universal(base_obj)
        
        incProgress(0.35, detail = "Choosing celldex reference")
        ref <- choose_reference(rownames(sce))
        common <- intersect(rownames(sce), rownames(ref))
        validate(need(length(common) > 500, "Too few overlapping genes with reference for SingleR."))
        
        incProgress(0.60, detail = "Running SingleR")
        pred <- SingleR::SingleR(test = sce[common, ], ref = ref[common, ], labels = ref$label.main)
        pruned <- SingleR::pruneScores(pred)
        labels <- pred$labels
        labels[is.na(pruned)] <- "Unknown"
        
        # attach back to Seurat
        base_obj$cell_type <- labels
        base_obj$major_type <- map_major_type(labels)
        
        # Ensure UMAP exists for plots
        use_npcs <- max(10, as.numeric(input$npcs %||% npcs))
        use_dims <- 1:use_npcs
        if (!"umap" %in% names(base_obj@reductions)) {
          if (!"pca" %in% names(base_obj@reductions)) {
            if (length(VariableFeatures(base_obj)) == 0) {
              base_obj <- NormalizeData(base_obj, verbose = FALSE)
              base_obj <- FindVariableFeatures(base_obj, selection.method = "vst",
                                               nfeatures = max(200, min(2000, nrow(base_obj))), verbose = FALSE)
            }
            base_obj <- RunPCA(base_obj, features = VariableFeatures(base_obj), npcs = use_npcs, verbose = FALSE)
          }
          base_obj <- RunUMAP(base_obj, dims = use_dims, verbose = FALSE)
        }
        
        # Build plots
        rv$p_umap_celltype <- DimPlot(base_obj, reduction = "umap", group.by = "cell_type",
                                      label = TRUE, repel = TRUE) +
          ggplot2::ggtitle("UMAP colored by cell type")
        
        rv$p_vln_rti_celltype <- tryCatch({
          VlnPlot(base_obj, features = "RTI", group.by = "cell_type", pt.size = 0) +
            ggplot2::ggtitle("RTI by cell type")
        }, error = function(e) ggplot2::ggplot() + ggplot2::ggtitle("RTI by cell type (RTI not computed yet)"))
        
        rv$p_vln_rti_majortype <- tryCatch({
          VlnPlot(base_obj, features = "RTI", group.by = "major_type", pt.size = 0) +
            ggplot2::ggtitle("RTI by major type")
        }, error = function(e) ggplot2::ggplot() + ggplot2::ggtitle("RTI by major type (RTI not computed yet)"))
        
        # Summaries
        df_meta <- base_obj@meta.data
        if (!is.null(df_meta$cell_type)) {
          rv$type_summary_df <- df_meta |>
            dplyr::group_by(cell_type) |>
            dplyr::summarise(
              n = dplyr::n(),
              RTI_median = stats::median(RTI, na.rm = TRUE),
              RTI_q25    = stats::quantile(RTI, 0.25, na.rm = TRUE),
              RTI_q75    = stats::quantile(RTI, 0.75, na.rm = TRUE),
              .groups = "drop"
            ) |>
            dplyr::arrange(dplyr::desc(RTI_median))
        } else {
          rv$type_summary_df <- data.frame(message = "cell_type metadata unavailable")
        }
        
        if (!is.null(df_meta$cell_type) && !is.null(df_meta$seurat_clusters)) {
          rv$cluster_type_map <- df_meta |>
            dplyr::mutate(seurat_clusters = as.character(seurat_clusters),
                          cell_type = as.character(cell_type)) |>
            dplyr::count(seurat_clusters, cell_type, name = "n") |>
            dplyr::group_by(seurat_clusters) |>
            dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) |>
            dplyr::ungroup() |>
            dplyr::arrange(as.numeric(seurat_clusters))
        } else {
          rv$cluster_type_map <- data.frame(message = "seurat_clusters or cell_type missing")
        }
        
        # keep object with added annotations
        rv$obj <- base_obj
        incProgress(1, detail = "Done")
      })
    })
    
    # Renderers for Cell types tab
    output$umap_celltype <- renderPlot({ req(rv$p_umap_celltype); rv$p_umap_celltype }, res = 120)
    output$vln_rti_celltype <- renderPlot({ req(rv$p_vln_rti_celltype); rv$p_vln_rti_celltype }, res = 120)
    output$vln_rti_majortype <- renderPlot({ req(rv$p_vln_rti_majortype); rv$p_vln_rti_majortype }, res = 120)
    output$tbl_type_summary <- renderTable({ req(rv$type_summary_df); rv$type_summary_df }, striped = TRUE, bordered = TRUE, digits = 3)
    output$tbl_cluster_type_map <- renderTable({ req(rv$cluster_type_map); rv$cluster_type_map }, striped = TRUE, bordered = TRUE, digits = 0)
    
    # Downloads for Cell types tab
    output$dl_umap_celltype <- downloadHandler(
      filename = function() paste0("UMAP_cell_types_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) ggplot2::ggsave(file, plot = rv$p_umap_celltype, width = 7, height = 5, dpi = 1000)
    )
    output$dl_vln_rti_celltype <- downloadHandler(
      filename = function() paste0("RTI_violin_cell_type_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) ggplot2::ggsave(file, plot = rv$p_vln_rti_celltype, width = 7, height = 5, dpi = 1000)
    )
    output$dl_vln_rti_majortype <- downloadHandler(
      filename = function() paste0("RTI_violin_major_type_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) ggplot2::ggsave(file, plot = rv$p_vln_rti_majortype, width = 7, height = 5, dpi = 1000)
    )
    output$dl_type_summary_csv <- downloadHandler(
      filename = function() paste0("RTI_type_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) utils::write.csv(rv$type_summary_df, file = file, row.names = FALSE)
    )
    output$dl_cluster_type_map_csv <- downloadHandler(
      filename = function() paste0("RTI_cluster_type_map_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content  = function(file) utils::write.csv(rv$cluster_type_map, file = file, row.names = FALSE)
    )
    
    # Expose updated Seurat object if the app wants it
    invisible(list(seurat = reactive(rv$obj)))
  })
}
