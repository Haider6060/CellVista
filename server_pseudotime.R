# server_pseudotime.R — FINAL with placeholder in values selector

pseudotime_server <- function(id, seurat_obj_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(cds = NULL, plot_umap = NULL, plot_states = NULL, plot_group = NULL)
    
    # ---------------- helpers ----------------
    build_cds_from_seurat <- function(obj) {
      assay_now <- if ("RNA" %in% names(obj@assays)) "RNA" else names(obj@assays)[1]
      SeuratObject::DefaultAssay(obj) <- assay_now
      
      layers_avail  <- tryCatch(SeuratObject::Layers(obj[[assay_now]]), error = function(e) character(0))
      counts_layers <- grep("^counts", layers_avail, value = TRUE)
      if (length(counts_layers) == 0) {
        def_layer <- tryCatch(SeuratObject::DefaultLayer(obj[[assay_now]]), error = function(e) NULL)
        if (is.null(def_layer) && length(layers_avail) > 0) def_layer <- layers_avail[1]
        counts_layers <- def_layer
      }
      
      mats  <- lapply(counts_layers, function(ly) Seurat::GetAssayData(obj, assay = assay_now, layer = ly))
      genes <- Reduce(intersect, lapply(mats, rownames))
      mats  <- lapply(mats, function(m) m[genes, , drop = FALSE])
      cts   <- do.call(cbind, mats)
      cts   <- cts[, !duplicated(colnames(cts)), drop = FALSE]
      
      cmeta  <- obj@meta.data
      common <- intersect(colnames(cts), rownames(cmeta))
      cts    <- cts[, common, drop = FALSE]
      cmeta  <- cmeta[common, , drop = FALSE]
      
      gmeta <- data.frame(gene_short_name = rownames(cts), row.names = rownames(cts))
      monocle3::new_cell_data_set(cts, cell_metadata = cmeta, gene_metadata = gmeta)
    }
    
    has_norm  <- function(x) grepl("(^|[^a-z])(normal|adjacent|adj|para.?cancer|control|ctrl|healthy|benign|nat|nml|non[- ]?tumou?r)([^a-z]|$)", x, perl = TRUE)
    has_tumor <- function(x) grepl("(^|[^a-z])(tumou?r|tumoural|tumoral|cancer|carcinoma|malignan|primary|metasta|lesion|neoplasm|adenocarcinoma|squamous|gbm|hcc|crc|melanoma|sarcoma|leukemia|lymphoma|myeloma|case)([^a-z]|$)", x, perl = TRUE)
    
    auto_group <- function(md) {
      L <- tolower
      cand_cell <- grep("(^sample$)|orig.ident|sample|group|condition|status|disease|type|source|title|tissue|phenotype|diagnos|histolog|project",
                        colnames(md), ignore.case = TRUE, value = TRUE)
      cand_cell <- unique(c(colnames(md)[tolower(colnames(md)) == "sample"], cand_cell))
      
      for (col in cand_cell) {
        s <- L(as.character(md[[col]]))
        if (any(has_norm(s), na.rm = TRUE))
          return(factor(ifelse(has_norm(s), "Normal", "Tumor"), levels = c("Normal","Tumor")))
      }
      
      cand_sample <- grep("(^sample$)|orig.ident|sample|patient|donor|run|library|batch|dataset",
                          colnames(md), ignore.case = TRUE, value = TRUE)
      if (length(cand_sample) == 0) {
        txt <- names(which(sapply(md, function(x) is.character(x) || is.factor(x))))
        small <- txt[sapply(txt, function(x) length(unique(md[[x]])) <= max(3, round(nrow(md)*0.1)))]
        cand_sample <- if (length(small)) small else txt[1]
      }
      sCol <- cand_sample[1]
      sid  <- as.character(md[[sCol]])
      uL   <- L(unique(sid))
      if (any(has_norm(uL))) {
        lab <- ifelse(has_norm(uL), "Normal", "Tumor")
        return(factor(setNames(lab, unique(sid))[sid], levels = c("Normal","Tumor")))
      }
      
      if (any(has_tumor(uL))) return(factor(rep("Tumor", nrow(md)), levels = c("Normal","Tumor")))
      factor(rep("Tumor", nrow(md)), levels = c("Normal","Tumor"))
    }
    
    apply_manual_group <- function(cds, col, normals) {
      md <- as.data.frame(SummarizedExperiment::colData(cds))
      vals <- as.character(md[[col]])
      g <- ifelse(vals %in% normals, "Normal", "Tumor")
      factor(g, levels = c("Normal","Tumor"))
    }
    
    choose_root_cells <- function(cds, max_n = 200) {
      md <- as.data.frame(SummarizedExperiment::colData(cds))
      roots <- rownames(md)[md$group %in% "Normal"]
      if (length(roots) == 0) {
        libs <- Matrix::colSums(SummarizedExperiment::assay(cds, "counts"))
        roots <- names(libs)[order(libs)][1L]
      } else {
        roots <- roots[seq_len(min(max_n, length(roots)))]
      }
      roots
    }
    
    # ---- Plots (no root/branch/leaf labels) ----
    make_plot_umap <- function(cds) {
      monocle3::plot_cells(
        cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE,
        label_cell_groups = FALSE, label_roots = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
        cell_size = 0.3, trajectory_graph_segment_size = 0.7, trajectory_graph_color = "black"
      ) +
        ggplot2::scale_color_gradientn(colors = c("#08306B","#2171B5","#6BAED6","#C6DBEF"), name = "Pseudotime") +
        ggplot2::theme_classic(base_size = 12) + ggplot2::theme(legend.position = "right")
    }
    
    make_plot_states <- function(cds) {
      SummarizedExperiment::colData(cds)$State <- factor(as.integer(as.factor(monocle3::partitions(cds))))
      monocle3::plot_cells(
        cds, color_cells_by = "State", show_trajectory_graph = TRUE,
        label_cell_groups = FALSE, label_roots = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
        cell_size = 0.3, trajectory_graph_segment_size = 0.7, trajectory_graph_color = "black"
      ) +
        ggplot2::scale_color_brewer(palette = "Set2", name = "State") +
        ggplot2::theme_classic(base_size = 12)
    }
    
    make_plot_group <- function(cds) {
      monocle3::plot_cells(
        cds, color_cells_by = "group", show_trajectory_graph = TRUE,
        label_cell_groups = FALSE, label_roots = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
        cell_size = 0.3, trajectory_graph_segment_size = 0.7, trajectory_graph_color = "black"
      ) +
        ggplot2::scale_color_manual(values = c("Normal"="#F8766D","Tumor"="#00BFC4"), name = "Group") +
        ggplot2::theme_classic(base_size = 12)
    }
    
    # ---------------- pipeline ----------------
    observeEvent(input$run, {
      obj <- seurat_obj_reactive()
      validate(need(!is.null(obj) && inherits(obj, "Seurat"),
                    "Please upload a valid Seurat object in the Data tab first."))
      
      withProgress(message = "Running pseudotime analysis…", value = 0, {
        incProgress(0.10, detail = "Building cell_data_set")
        cds <- build_cds_from_seurat(obj)
        
        incProgress(0.30, detail = "Preprocess & UMAP")
        cds <- monocle3::preprocess_cds(cds, num_dim = 50)
        cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
        
        incProgress(0.55, detail = "Cluster & Learn graph")
        set.seed(1)
        cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
        cds <- monocle3::learn_graph(cds, use_partition = FALSE,
                                     learn_graph_control = list(ncenter = 200))
        
        incProgress(0.75, detail = "Order cells (rooted)")
        md <- as.data.frame(SummarizedExperiment::colData(cds))
        auto_g <- auto_group(md)
        SummarizedExperiment::colData(cds)$group <- auto_g
        roots <- choose_root_cells(cds, max_n = 200)
        cds <- monocle3::order_cells(cds, root_cells = roots, reduction_method = "UMAP")
        SummarizedExperiment::colData(cds)$State <- factor(as.integer(as.factor(monocle3::partitions(cds))))
        
        # manual mapping UI (with placeholder)
        cols <- names(md)[sapply(md, function(x) is.character(x) || is.factor(x))]
        updateSelectInput(session, "group_col", choices = cols, selected = cols[1])
        output$group_values_ui <- renderUI({
          req(input$group_col)
          vals <- unique(as.character(md[[input$group_col]]))
          selectizeInput(
            ns("normal_vals"),
            "Values treated as Normal:",
            choices = sort(vals),
            multiple = TRUE,
            options = list(placeholder = "Click to choose one or more values…")
          )
        })
        
        rv$cds <- cds
        rv$plot_umap <- rv$plot_states <- rv$plot_group <- NULL
        incProgress(1)
        cl <- table(SummarizedExperiment::colData(cds)$group)
        if (length(cl[cl > 0]) < 2) {
          showNotification("Only one class found automatically. Use the panel above to mark Normal values; others become Tumor.",
                           type = "warning", duration = NULL)
        } else {
          showNotification("Pseudotime analysis completed.", type = "message", duration = 5)
        }
      })
    }, ignoreInit = TRUE)
    
    observeEvent(input$apply_group, {
      req(rv$cds, input$group_col)
      g <- apply_manual_group(rv$cds, input$group_col, input$normal_vals %||% character(0))
      SummarizedExperiment::colData(rv$cds)$group <- g
      if (!is.null(rv$plot_group)) rv$plot_group <- make_plot_group(rv$cds)
      showNotification("Group mapping applied.", type = "message", duration = 4)
    })
    
    # on-demand plots
    observeEvent(input$show_umap,   { req(rv$cds); rv$plot_umap   <- make_plot_umap(rv$cds) })
    observeEvent(input$show_states, { req(rv$cds); rv$plot_states <- make_plot_states(rv$cds) })
    observeEvent(input$show_group,  { req(rv$cds); rv$plot_group  <- make_plot_group(rv$cds) })
    
    output$plot_umap   <- renderPlot({ req(rv$plot_umap);   rv$plot_umap },   res = 120)
    output$plot_states <- renderPlot({ req(rv$plot_states); rv$plot_states }, res = 120)
    output$plot_group  <- renderPlot({ req(rv$plot_group);  rv$plot_group },  res = 120)
    
    # downloads (1000 dpi)
    output$dl_umap <- downloadHandler(
      filename = function() paste0("pseudotime_umap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) {
        req(rv$cds)
        p <- if (!is.null(rv$plot_umap)) rv$plot_umap else make_plot_umap(rv$cds)
        ggplot2::ggsave(file, plot = p, width = 6, height = 5, dpi = 1000)
      }
    )
    output$dl_states <- downloadHandler(
      filename = function() paste0("trajectory_states_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) {
        req(rv$cds)
        p <- if (!is.null(rv$plot_states)) rv$plot_states else make_plot_states(rv$cds)
        ggplot2::ggsave(file, plot = p, width = 6, height = 5, dpi = 1000)
      }
    )
    output$dl_group <- downloadHandler(
      filename = function() paste0("trajectory_group_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content  = function(file) {
        req(rv$cds)
        p <- if (!is.null(rv$plot_group)) rv$plot_group else make_plot_group(rv$cds)
        ggplot2::ggsave(file, plot = p, width = 6, height = 5, dpi = 1000)
      }
    )
  })
}
