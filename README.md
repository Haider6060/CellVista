# CellVista: Single-Cell Transcriptome Analysis with Novel Indices  

## Overview  

**CellVista** is a user-friendly **Shiny application** designed for researchers and bioinformaticians working on **single-cell RNA-seq (scRNA-seq)** analysis.  
The app integrates **standard scRNA-seq workflows** (QC, clustering, pseudotime) with **three novel analysis modules**:  

- **CellEntropy** → quantifies transcriptional diversity across cells  
- **RTI (Regulatory Turbulence Index)** → measures transcriptional instability and regulatory noise  
- **BSI (Boundary Sharpness Index)** → a new metric that quantifies the clarity of transcriptional boundaries between clusters or states  

With **CellVista**, users can easily explore single-cell landscapes, identify stable vs. transitional cell states, and gain deeper insights into tumor heterogeneity and microenvironment dynamics.  

---

## ✅ Key Features  

- Data input via `.rds` files (Seurat v5 objects)  
- Automated **QC & normalization**  
- **Clustering and UMAP** visualization  
- **Pseudotime trajectory inference**  
- **Entropy landscape transcriptome** analysis  
- **RTI (Regulatory Turbulence Index)** for instability quantification  
- **BSI (Boundary Sharpness Index)** for cluster boundary sharpness  
- High-resolution plots and CSV export  

---

## 📂 Input Requirements  

- `.rds` file (Seurat v5 format)  

### ➤ Preparing Your Data  

```r
# Example (from GEO matrix):
seurat_obj <- CreateSeuratObject(counts = your_matrix)
saveRDS(seurat_obj, file = "your_dataset.rds")
```  

---

## 🧪 Datasets Tested  

CellVista has been successfully tested on multiple real-world datasets:  

| Dataset Type | Source | Description |
|--------------|--------|-------------|
| Lung cancer (GEO) | GEO | Primary development dataset |
| Breast cancer (GEO) | GEO | Validation dataset |
| Pancreatic cancer (GEO) | GEO | Tumor microenvironment |
| Brain tumor (GEO) | GEO | Glioblastoma test dataset |  
|Hypothalamus development (GEO) | GEO | Non-cancer developmental dataset (mouse, 11 timepoints)

All datasets processed successfully, confirming **robustness and generalizability**.  

---

## 🎯 Quick Start  

### Option 1: Run in RStudio  

1. Open `app.R`  
2. Click **Run App**  

### Option 2: R Console  

```r
shiny::runApp("your_app_folder_path")
```  

---

## 📦 Demo Files Included  

| File Name | Description |
|-----------|-------------|
| `demo_lung_seurat.rds` | GEO lung cancer subsample |
| `demo_breast_seurat.rds` | Breast cancer subsample |
| `demo_pancreas_seurat.rds` | Pancreatic cancer subsample |
| `demo_brain_seurat.rds` | Glioblastoma subsample |  
| demo_hypothalamus_seurat.rds | Mouse hypothalamus developmental subsample (non-cancer)
> ⚠️ These are 200-cell subsamples to meet GitHub limits. Full datasets used for testing are available on request.

---


## 📄 License  

MIT License  

---

## 📬 Contact  

**Developer:** Haider  
📧 haider@emails.bjut.edu.cn  

---

## 📚 Citation  

If **CellVista** contributes to your research, please cite this tool in your publication.  


