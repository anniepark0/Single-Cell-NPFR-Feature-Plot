library(Seurat)
library(Rutils)

# Download Thirst2 dataset from GEO (GSE207799)
f <- "GSE207799_Thirst2_SCT_trimPlus.Robj"
if (!file.exists(f)) {
  options(timeout = 3600)
  download.file(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207799/suppl/", f, ".gz"),
                paste0(f, ".gz"), mode = "wb")
  R.utils::gunzip(paste0(f, ".gz"))
}
thirst <- UpdateSeuratObject(get(load(f)))


# Subset Vmat+ DAT+ (dopaminergic) cells and re-cluster
da.cells <- subset(thirst, subset = Vmat > 4 & DAT > 3)
da <- CreateSeuratObject(da.cells@assays$RNA@counts, meta.data = da.cells@meta.data["orig.ident"],
                         min.cells = 1, min.features = 200, project = "VmatDAT")
da <- NormalizeData(da)
da <- ScaleData(da, vars.to.regress = "orig.ident")
da <- FindVariableFeatures(da)
da <- RunPCA(da)
da <- RunUMAP(da, dims = 1:50)
da <- FindNeighbors(da)
da <- FindClusters(da, resolution = 4)

# Plots
DimPlot(da, label = TRUE) + NoLegend()
FeaturePlot(da, "NPFR", cols = c("grey80", "#5FAFA4"), order = TRUE, alpha = 0.6) + NoAxes()


