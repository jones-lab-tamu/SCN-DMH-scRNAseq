# Copyright (c) Mustafa Ozturgut
# Analysis done for the Jones Lab, Texas A&M University

# Whole-droplet EGFP analysis of the Jones Lab SCN 10x Flex data (submitted QC). Makes the paper figures
# (Supplementary Figure 3a-b, Figure 4b-f), two EGFP+ vs EGFP- volcano plots, the EGFP+ neurons per Ward group at three
# EGFP cutoffs and the Morris 2021 / Wen 2020 gene-list scores on the putative SCN UMAP. All figures are saved as PDF.

# =============================================================================
# 1. Settings and libraries
# =============================================================================
# Working folder: holds the Cell Ranger 9 raw matrix (Jones Lab data/) and the Morris/Wen gene lists (reference_data/)
setwd("D:/Mustafa Active analysis/Mustafa/Collaboration/Jones Lab")

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(fmsb)
library(ggrepel)
library(readxl)
library(ggrastr)

# Output folders, one per figure set
dir.create(file.path("whole droplet paper figures", "Supplementary Figure 3"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("whole droplet paper figures", "Figure 4"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("whole droplet paper figures", "Volcano Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("whole droplet paper figures", "EGFP per Ward Group"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("whole droplet paper figures", "Morris and Wen Gene List Scores"), recursive = TRUE, showWarnings = FALSE)

# Plot titles: bold, size 12
bold_title <- theme(plot.title = element_text(face = "bold", size = 12))

# Probe barcode letters 17-24 (10x CG000621 Rev D)
barcode_to_portion <- c(ACTTTAGG = "BC001", AACGGGAA = "BC002", AGTAGGCT = "BC003", ATGTTGAC = "BC004")

# =============================================================================
# 2. Load counts and call EGFP per droplet
# =============================================================================
# Unfiltered matrix; barcodes with 200+ genes, genes in 3+ barcodes
jones_counts <- Read10X(data.dir = "Jones Lab data")
jones_seurat <- CreateSeuratObject(counts = jones_counts, project = "JonesLab", min.cells = 3, min.features = 200)
# Number of genes and of barcodes with 200+ genes (16,888 and 46,960)
print(dim(jones_seurat))

# A Flex barcode is 24 letters: letters 1-16 name the droplet, letters 17-24 the probe barcode (BC001-BC004).
# EGFP is summed over the whole droplet, including barcodes below the 200-gene cut, so the droplet name is taken
# for every raw barcode (to sum EGFP) and again for the cells kept in jones_seurat (to look the sum up)
egfp_rows <- c("EGFP-1", "EGFP-2", "EGFP-3")
raw_egfp_umi <- Matrix::colSums(jones_counts[egfp_rows, ])
raw_droplet <- substr(colnames(jones_counts), 1, 16)
egfp_per_droplet <- tapply(raw_egfp_umi, raw_droplet, sum)
jones_seurat$droplet <- substr(colnames(jones_seurat), 1, 16)
jones_seurat$probe_barcode <- substr(colnames(jones_seurat), 17, 24)

# Number of cells in each droplet; TRUE if the cell is the only one in its droplet
cells_per_droplet <- table(jones_seurat$droplet)
cells_in_droplet <- as.numeric(cells_per_droplet[jones_seurat$droplet])
jones_seurat$no_other_cell <- cells_in_droplet == 1
# EGFP UMIs summed over the cell's whole droplet
jones_seurat$egfp_droplet_umi <- as.numeric(egfp_per_droplet[jones_seurat$droplet])
# Whole-droplet EGFP+ call: alone in the droplet and at least 1 EGFP UMI in the droplet
jones_seurat$gfp_positive <- jones_seurat$no_other_cell & jones_seurat$egfp_droplet_umi >= 1

# Mitochondrial percent (genes starting with mt-)
jones_seurat[["percent.mt"]] <- PercentageFeatureSet(jones_seurat, pattern = "^mt-")

# Probe barcode letters to portion names BC001-BC004
portion_names <- as.character(barcode_to_portion[jones_seurat$probe_barcode])
jones_seurat$portion <- factor(portion_names, levels = barcode_to_portion)
print(table(jones_seurat$portion, useNA = "ifany"))

# =============================================================================
# 3. Quality control filter
# =============================================================================
# QC filter: 1,000 to 8,000 genes and mitochondrial RNA below 10%
jones_seurat <- subset(jones_seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 10)
# Cells kept after QC (40,854)
print(ncol(jones_seurat))
print(table(jones_seurat$portion))

# =============================================================================
# 4. SingleR cell-type labels (MouseRNAseqData)
# =============================================================================
jones_seurat[["RNA"]] <- JoinLayers(jones_seurat[["RNA"]])
jones_seurat <- NormalizeData(jones_seurat, verbose = FALSE)
ref_mouse <- MouseRNAseqData()
sce_jones <- as.SingleCellExperiment(jones_seurat, assay = "RNA")
pred_jones <- SingleR(test = sce_jones, ref = ref_mouse, labels = ref_mouse$label.main)
jones_seurat$pruned_labels <- pred_jones$pruned.labels
print(table(jones_seurat$pruned_labels, useNA = "ifany"))

# =============================================================================
# 5. All cells: PCA, clustering and UMAP
# =============================================================================
jones_seurat <- FindVariableFeatures(jones_seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
jones_seurat <- ScaleData(jones_seurat, verbose = FALSE)
jones_seurat <- RunPCA(jones_seurat, npcs = 30, verbose = FALSE)
jones_seurat <- FindNeighbors(jones_seurat, dims = 1:30, verbose = FALSE)
jones_seurat <- FindClusters(jones_seurat, resolution = 0.5, verbose = FALSE)
jones_seurat <- RunUMAP(jones_seurat, dims = 1:30, verbose = FALSE)
jones_seurat$seurat_clusters <- Idents(jones_seurat)
# Number of all-cell clusters (31) and whole-droplet EGFP+ cells after QC (410)
print(length(unique(jones_seurat$seurat_clusters)))
print(sum(jones_seurat$gfp_positive))

# =============================================================================
# 6. Supplementary Figure 3a: SingleR labels
# =============================================================================
# One colour per SingleR label
cell_type_colors <- c(Neurons = "#08519c", Macrophages = "#A63603", "Epithelial cells" = "#54278F", Astrocytes = "#F16913",
                                Microglia = "#FEC44F", Oligodendrocytes = "#CB181D", "Endothelial cells" = "#238B45", Fibroblasts = "#FFD92F",
                                "T cells" = "#4292C6", "NK cells" = "#DD3497", Cardiomyocytes = "#F768A1", Monocytes = "#67001F")

# One row per cell: UMAP coordinates (stored by RunUMAP), SingleR label and EGFP call
cell_umap <- FetchData(jones_seurat, vars = c("umap_1", "umap_2", "pruned_labels", "gfp_positive"))
# Labels in the order they first appear; cells without a label are left out of the legend
cell_type_levels <- unique(cell_umap$pruned_labels[!is.na(cell_umap$pruned_labels)])
cell_umap$cell_type <- factor(cell_umap$pruned_labels, levels = cell_type_levels)
# EGFP+ neurons, drawn in black on top
egfp_neurons <- cell_umap[cell_umap$pruned_labels %in% "Neurons" & cell_umap$gfp_positive, ]
egfp_neurons$legend <- "EGFP+ neurons"

p_3a <- ggplot(cell_umap, aes(x = umap_1, y = umap_2, color = cell_type)) +
  geom_point(size = 0.8, alpha = 0.82, stroke = 0, shape = 16) +
  scale_color_manual(values = cell_type_colors[cell_type_levels], limits = cell_type_levels, breaks = cell_type_levels, name = "", drop = FALSE) +
  geom_point(data = egfp_neurons, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 3.5, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(size = 2.45, color = "black"))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Cell Types and EGFP+ Neurons", x = "UMAP 1", y = "UMAP 2") +
  bold_title
print(p_3a)

# =============================================================================
# 7. Neuron subset
# =============================================================================
# The 16 all-cell clusters that are neurons
neuron_clusters <- c("0", "2", "3", "4", "5", "7", "8", "9", "10", "11", "12", "15", "16", "18", "20", "22")
neuron_seurat <- subset(jones_seurat, idents = neuron_clusters)

# Neuron recluster (PCs 1-20, default resolution)
neuron_seurat[["RNA"]] <- JoinLayers(neuron_seurat[["RNA"]])
neuron_seurat <- NormalizeData(neuron_seurat, verbose = FALSE)
neuron_seurat <- FindVariableFeatures(neuron_seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
neuron_seurat <- ScaleData(neuron_seurat, features = VariableFeatures(neuron_seurat), verbose = FALSE)
neuron_seurat <- RunPCA(neuron_seurat, features = VariableFeatures(neuron_seurat), verbose = FALSE)
neuron_seurat <- FindNeighbors(neuron_seurat, dims = 1:20, verbose = FALSE)
neuron_seurat <- FindClusters(neuron_seurat, verbose = FALSE)
neuron_seurat <- RunUMAP(neuron_seurat, dims = 1:20, verbose = FALSE)
# Number of neurons (28,894)
print(ncol(neuron_seurat))

# =============================================================================
# 8. Neuron re-clustering (PCs 1-15, resolution 0.8)
# =============================================================================
neuron_seurat[["RNA"]] <- JoinLayers(neuron_seurat[["RNA"]])
neuron_seurat <- NormalizeData(neuron_seurat, verbose = FALSE)
neuron_seurat <- FindVariableFeatures(neuron_seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
neuron_seurat <- ScaleData(neuron_seurat, features = rownames(neuron_seurat), verbose = FALSE)
neuron_seurat <- RunPCA(neuron_seurat, features = VariableFeatures(neuron_seurat), verbose = FALSE)
neuron_seurat <- FindNeighbors(neuron_seurat, dims = 1:15, verbose = FALSE)
set.seed(42)
neuron_seurat <- RunUMAP(neuron_seurat, dims = 1:15, verbose = FALSE)
set.seed(42)
neuron_seurat <- FindClusters(neuron_seurat, resolution = 0.8, verbose = FALSE)
neuron_seurat$clusters_res_0_8 <- Idents(neuron_seurat)
# Number of neuron clusters (27) and whole-droplet EGFP+ neurons (379)
print(length(unique(neuron_seurat$clusters_res_0_8)))
print(sum(neuron_seurat$gfp_positive))

# =============================================================================
# 9. Supplementary Figure 3b: Sst+ and Trh+ neurons on the neuron map
# =============================================================================
# One row per neuron: UMAP coordinates, Sst and Trh RNA (log-normalised) and EGFP call
neuron_umap <- FetchData(neuron_seurat, vars = c("umap_1", "umap_2", "Sst", "Trh", "gfp_positive"))
# Sst+ or Trh+ if any RNA; a neuron with both is shown as Sst+
neuron_umap$highlight <- "Other neurons"
neuron_umap$highlight[neuron_umap$Sst > 0] <- "Sst+ (including Sst+ Trh+)"
neuron_umap$highlight[neuron_umap$Sst == 0 & neuron_umap$Trh > 0] <- "Trh+ only"
neuron_umap$highlight <- factor(neuron_umap$highlight, levels = c("Other neurons", "Sst+ (including Sst+ Trh+)", "Trh+ only"))
# EGFP+ neurons, drawn in black on top
egfp_neuron_umap <- neuron_umap[neuron_umap$gfp_positive, ]
egfp_neuron_umap$legend <- "EGFP+ neurons"

p_3b <- ggplot(neuron_umap, aes(x = umap_1, y = umap_2, color = highlight)) +
  geom_point(size = 0.35, alpha = 0.88) +
  scale_color_manual(values = c("Other neurons" = "grey92", "Sst+ (including Sst+ Trh+)" = "#C77CFF", "Trh+ only" = "#66C2A5"),
                     limits = levels(neuron_umap$highlight), breaks = c("Sst+ (including Sst+ Trh+)", "Trh+ only"), drop = FALSE, name = NULL) +
  geom_point(data = egfp_neuron_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 3.5, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(shape = 16, size = 2.45, color = "black", alpha = 1))) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Neuron Clusters Screened for Removal") +
  bold_title
print(p_3b)

# =============================================================================
# 10. Putative SCN neurons (PCs 1-15, resolution 0.5)
# =============================================================================
# Remove the six neuron clusters screened out in Supplementary Figure 3b
removed_clusters <- c("4", "26", "23", "12", "17", "20")
scn_seurat <- subset(neuron_seurat, idents = removed_clusters, invert = TRUE)
scn_seurat[["RNA"]] <- JoinLayers(scn_seurat[["RNA"]])
scn_seurat <- NormalizeData(scn_seurat, verbose = FALSE)
scn_seurat <- FindVariableFeatures(scn_seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
scn_seurat <- ScaleData(scn_seurat, features = rownames(scn_seurat), verbose = FALSE)
scn_seurat <- RunPCA(scn_seurat, features = VariableFeatures(scn_seurat), verbose = FALSE)
scn_seurat <- FindNeighbors(scn_seurat, dims = 1:15, verbose = FALSE)
set.seed(42)
scn_seurat <- FindClusters(scn_seurat, resolution = 0.5, verbose = FALSE)
set.seed(42)
scn_seurat <- RunUMAP(scn_seurat, dims = 1:15, verbose = FALSE)
scn_seurat$seurat_clusters <- Idents(scn_seurat)
# Number of putative SCN neurons (24,179), neurons alone in their droplet (19,150) and whole-droplet EGFP+ neurons (355)
print(ncol(scn_seurat))
print(sum(scn_seurat$no_other_cell))
print(sum(scn_seurat$gfp_positive))

# One row per putative SCN neuron: UMAP coordinates and EGFP call, for the EGFP+ overlay on the UMAP panels below
scn_umap <- FetchData(scn_seurat, vars = c("umap_1", "umap_2", "gfp_positive"))
egfp_scn_umap <- scn_umap[scn_umap$gfp_positive, ]
egfp_scn_umap$legend <- "EGFP+ neurons"

# =============================================================================
# 11. Figure 4b: seven SCN genes on the UMAP
# =============================================================================
# One row per neuron: UMAP coordinates, EGFP call and log-normalised RNA of the seven genes
scn_genes <- c("Avp", "Vip", "Cck", "Nms", "Vipr2", "Grp", "Prokr2")
gene_umap <- FetchData(scn_seurat, vars = c("umap_1", "umap_2", "gfp_positive", scn_genes))

# Each gene scaled 0-1 across neurons: minus its minimum, divided by its range
gene_umap$Avp_01 <- (gene_umap$Avp - min(gene_umap$Avp)) / (max(gene_umap$Avp) - min(gene_umap$Avp))
gene_umap$Vip_01 <- (gene_umap$Vip - min(gene_umap$Vip)) / (max(gene_umap$Vip) - min(gene_umap$Vip))
gene_umap$Cck_01 <- (gene_umap$Cck - min(gene_umap$Cck)) / (max(gene_umap$Cck) - min(gene_umap$Cck))
gene_umap$Nms_01 <- (gene_umap$Nms - min(gene_umap$Nms)) / (max(gene_umap$Nms) - min(gene_umap$Nms))
gene_umap$Vipr2_01 <- (gene_umap$Vipr2 - min(gene_umap$Vipr2)) / (max(gene_umap$Vipr2) - min(gene_umap$Vipr2))
gene_umap$Grp_01 <- (gene_umap$Grp - min(gene_umap$Grp)) / (max(gene_umap$Grp) - min(gene_umap$Grp))
gene_umap$Prokr2_01 <- (gene_umap$Prokr2 - min(gene_umap$Prokr2)) / (max(gene_umap$Prokr2) - min(gene_umap$Prokr2))

# One table per gene, sorted from low to high so the neurons with the most RNA are drawn last (on top)
avp_umap <- gene_umap[order(gene_umap$Avp_01), ]
vip_umap <- gene_umap[order(gene_umap$Vip_01), ]
cck_umap <- gene_umap[order(gene_umap$Cck_01), ]
nms_umap <- gene_umap[order(gene_umap$Nms_01), ]
vipr2_umap <- gene_umap[order(gene_umap$Vipr2_01), ]
grp_umap <- gene_umap[order(gene_umap$Grp_01), ]
prokr2_umap <- gene_umap[order(gene_umap$Prokr2_01), ]

# One panel per gene, EGFP+ neurons in black (the EGFP+ legend is shown on the Avp panel only)
p_avp <- ggplot(avp_umap, aes(x = umap_1, y = umap_2, color = Avp_01)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(2, 2, 2, 2)) +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = unit(1.35, "inch"), barheight = unit(0.22, "inch"), title.position = "top", title.theme = element_text(size = 8))) +
  scale_color_gradientn(colours = c("#EDEDED", "#0571b0"), name = "Norm. RNA (0-1)", limits = c(0, 1), oob = scales::squish) +
  ggtitle("Avp") +
  theme(plot.title = element_text(face = "italic", color = "#0571b0", size = 11, hjust = 0.5))
p_vip <- ggplot(vip_umap, aes(x = umap_1, y = umap_2, color = Vip_01)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(2, 2, 2, 2)) +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = unit(1.35, "inch"), barheight = unit(0.22, "inch"), title.position = "top", title.theme = element_text(size = 8))) +
  guides(shape = "none") +
  scale_color_gradientn(colours = c("#EDEDED", "#00BFC4"), name = "Norm. RNA (0-1)", limits = c(0, 1), oob = scales::squish) +
  ggtitle("Vip") +
  theme(plot.title = element_text(face = "italic", color = "#00BFC4", size = 11, hjust = 0.5))
p_cck <- ggplot(cck_umap, aes(x = umap_1, y = umap_2, color = Cck_01)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(2, 2, 2, 2)) +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = unit(1.35, "inch"), barheight = unit(0.22, "inch"), title.position = "top", title.theme = element_text(size = 8))) +
  guides(shape = "none") +
  scale_color_gradientn(colours = c("#EDEDED", "#F8766D"), name = "Norm. RNA (0-1)", limits = c(0, 1), oob = scales::squish) +
  ggtitle("Cck") +
  theme(plot.title = element_text(face = "italic", color = "#F8766D", size = 11, hjust = 0.5))
p_nms <- ggplot(nms_umap, aes(x = umap_1, y = umap_2, color = Nms_01)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(2, 2, 2, 2)) +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = unit(1.35, "inch"), barheight = unit(0.22, "inch"), title.position = "top", title.theme = element_text(size = 8))) +
  guides(shape = "none") +
  scale_color_gradientn(colours = c("#EDEDED", "#C77CFF"), name = "Norm. RNA (0-1)", limits = c(0, 1), oob = scales::squish) +
  ggtitle("Nms") +
  theme(plot.title = element_text(face = "italic", color = "#C77CFF", size = 11, hjust = 0.5))
p_vipr2 <- ggplot(vipr2_umap, aes(x = umap_1, y = umap_2, color = Vipr2_01)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(2, 2, 2, 2)) +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = unit(1.35, "inch"), barheight = unit(0.22, "inch"), title.position = "top", title.theme = element_text(size = 8))) +
  guides(shape = "none") +
  scale_color_gradientn(colours = c("#EDEDED", "#E6AB02"), name = "Norm. RNA (0-1)", limits = c(0, 1), oob = scales::squish) +
  ggtitle("Vipr2") +
  theme(plot.title = element_text(face = "italic", color = "#E6AB02", size = 11, hjust = 0.5))
p_grp <- ggplot(grp_umap, aes(x = umap_1, y = umap_2, color = Grp_01)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(2, 2, 2, 2)) +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = unit(1.35, "inch"), barheight = unit(0.22, "inch"), title.position = "top", title.theme = element_text(size = 8))) +
  guides(shape = "none") +
  scale_color_gradientn(colours = c("#EDEDED", "#A6761D"), name = "Norm. RNA (0-1)", limits = c(0, 1), oob = scales::squish) +
  ggtitle("Grp") +
  theme(plot.title = element_text(face = "italic", color = "#A6761D", size = 11, hjust = 0.5))
p_prokr2 <- ggplot(prokr2_umap, aes(x = umap_1, y = umap_2, color = Prokr2_01)) +
  geom_point(size = 0.35, alpha = 0.9) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(2, 2, 2, 2)) +
  guides(color = guide_colorbar(direction = "horizontal", barwidth = unit(1.35, "inch"), barheight = unit(0.22, "inch"), title.position = "top", title.theme = element_text(size = 8))) +
  guides(shape = "none") +
  scale_color_gradientn(colours = c("#EDEDED", "#984EA3"), name = "Norm. RNA (0-1)", limits = c(0, 1), oob = scales::squish) +
  ggtitle("Prokr2") +
  theme(plot.title = element_text(face = "italic", color = "#984EA3", size = 11, hjust = 0.5))
p_4b <- wrap_plots(p_avp, p_vip, p_cck, p_nms, p_vipr2, p_grp, p_prokr2, ncol = 3) +
  plot_annotation(title = "Seven SCN Genes in Putative SCN Neurons", theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12)))
print(p_4b)

# =============================================================================
# 12. Figure 4c: neurons with detectable RNA for the seven genes
# =============================================================================
# Percent of neurons with any RNA of each gene, among all putative SCN neurons and among EGFP+ neurons
detected <- gene_umap[, scn_genes] > 0
pct_all <- 100 * colMeans(detected)
pct_egfp <- 100 * colMeans(detected[gene_umap$gfp_positive, ])
n_neurons_text <- format(nrow(gene_umap), big.mark = ",")
population_all <- paste0("All putative SCN neurons (n = ", n_neurons_text, ")")
population_egfp <- paste0("EGFP+ neurons (n = ", sum(gene_umap$gfp_positive), ")")
gene_column <- rep(scn_genes, 2)
population_column <- rep(c(population_all, population_egfp), each = 7)
pct_4c <- data.frame(gene = factor(gene_column, levels = scn_genes), population = population_column, percent = c(pct_all, pct_egfp))
pct_4c$percent_label <- formatC(pct_4c$percent, format = "f", digits = 1)
print(pct_4c)
fill_4c <- c("grey70", "black")
names(fill_4c) <- c(population_all, population_egfp)

p_4c <- ggplot(pct_4c, aes(x = gene, y = percent, fill = population)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = percent_label), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = fill_4c, name = NULL) +
  labs(x = "Gene", y = "% of neurons with detectable RNA", title = "Neurons with Detectable RNA for Seven SCN Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  bold_title
print(p_4c)

# =============================================================================
# 13. Figure 4d-f: Ward groups
# =============================================================================
# Mean RNA of the variable genes in each putative SCN cluster (genes in rows, clusters in columns)
cluster_avg_list <- AverageExpression(scn_seurat, features = VariableFeatures(scn_seurat), group.by = "seurat_clusters", layer = "data")
cluster_avg <- as.matrix(cluster_avg_list$RNA)
# z-score of each gene across clusters
clusters_by_genes_z <- scale(t(cluster_avg))
cluster_z <- t(clusters_by_genes_z)
cluster_z[is.na(cluster_z)] <- 0
# Drop the g prefix Seurat adds to numeric cluster names
colnames(cluster_z) <- sub("^g", "", colnames(cluster_z))

# Ward.D2 tree of the clusters (Euclidean distance between their z-score profiles), cut into nine groups
cluster_dist <- dist(t(cluster_z), method = "euclidean")
ward_tree <- hclust(cluster_dist, method = "ward.D2")
ward_group_of_cluster <- cutree(ward_tree, k = 9)

# Tree with its leaves reordered so each Ward group sits together
ward_dend <- as.dendrogram(ward_tree)
leaf_group <- ward_group_of_cluster[labels(ward_dend)]
ward_dend <- reorder(ward_dend, leaf_group, agglo.FUN = mean)
# First and last leaf of each Ward group along the reordered tree, for the brackets under it
leaf_group_reordered <- ward_group_of_cluster[labels(ward_dend)]
leaf_position <- 1:length(leaf_group_reordered)
group_left <- tapply(leaf_position, leaf_group_reordered, min) - 0.5
group_right <- tapply(leaf_position, leaf_group_reordered, max) + 0.5
group_middle <- (group_left + group_right) / 2
bracket_y <- -0.20 * max(ward_tree$height)
tick_y <- -0.16 * max(ward_tree$height)
bracket_left <- group_left + 0.15
bracket_right <- group_right - 0.15

# Figure 4d: Ward tree with the nine groups bracketed under it
ward_title <- paste0("Ward Clustering of ", ncol(cluster_z), " Putative SCN Clusters into 9 Groups")
par(mar = c(6, 4, 4, 2))
plot(ward_dend, main = ward_title, xlab = "", ylab = "Height (Ward distance)")
mtext("Putative SCN clusters (resolution 0.5)", side = 1, line = 4.5, cex = 0.8)
segments(bracket_left, bracket_y, bracket_right, bracket_y, lwd = 1.5, xpd = TRUE)
segments(bracket_left, bracket_y, bracket_left, tick_y, lwd = 1.5, xpd = TRUE)
segments(bracket_right, bracket_y, bracket_right, tick_y, lwd = 1.5, xpd = TRUE)
text(group_middle, bracket_y, labels = paste0("G", names(group_left)), pos = 1, cex = 0.8, xpd = TRUE)

# Each neuron's Ward group, looked up from its cluster
ward_group_per_cell <- ward_group_of_cluster[as.character(scn_seurat$seurat_clusters)]
group_text <- paste("Group", ward_group_per_cell)
scn_seurat$group_label <- factor(group_text, levels = paste("Group", 1:9))
print(table(scn_seurat$group_label))

# Markers of each Ward group (Wilcoxon, RNA, positive only); top five per group by adjusted p
Idents(scn_seurat) <- "group_label"
ward_markers <- FindAllMarkers(scn_seurat, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
Idents(scn_seurat) <- "seurat_clusters"
ward_top5 <- arrange(ward_markers, cluster, p_val_adj)
ward_top5 <- group_by(ward_top5, cluster)
ward_top5 <- slice_head(ward_top5, n = 5)
ward_top5 <- ungroup(ward_top5)
print(ward_top5[, c("cluster", "gene", "avg_log2FC", "p_val_adj")], n = nrow(ward_top5))

# Figure 4d UMAP: Ward groups, EGFP+ neurons in black
p_4d <- DimPlot(scn_seurat, reduction = "umap", group.by = "group_label", pt.size = 0.5) +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, color = "black", size = 2.45) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), name = NULL) +
  guides(color = guide_legend(order = 1, override.aes = list(size = 3.5)),
         shape = guide_legend(order = 2, override.aes = list(size = 2.45, color = "black"))) +
  labs(title = "Ward Groups of Putative SCN Neurons", x = "UMAP 1", y = "UMAP 2") +
  theme_bw() +
  bold_title
print(p_4d)

# Figure 4e: top five markers per group, in group order, with Prokr2 first and Kcnip1 last
ward_top5_sorted <- arrange(ward_top5, cluster, p_val_adj, desc(avg_log2FC))
top_genes <- unique(ward_top5_sorted$gene)
marker_genes <- c("Prokr2", setdiff(top_genes, c("Prokr2", "Kcnip1")), "Kcnip1")
p_4e <- DotPlot(scn_seurat, features = rev(marker_genes), group.by = "group_label", col.min = -2.5, col.max = 2.5, scale.min = 0, scale.max = 100, cols = c("lightgrey", "blue")) +
  coord_flip() +
  labs(x = "Gene", y = NULL, title = "Top Marker Genes of Each Ward Group") +
  guides(colour = guide_colorbar(title = "Scaled average\nexpression"), size = guide_legend(title = "Percent detected")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = "italic")) +
  bold_title
print(p_4e)

# EGFP+ neurons in each Ward group, and each group's share of all EGFP+ neurons
group_names <- levels(scn_seurat$group_label)
egfp_group_table <- table(scn_seurat$group_label[scn_seurat$gfp_positive])
egfp_per_group <- as.numeric(egfp_group_table)
egfp_share <- 100 * egfp_per_group / sum(egfp_per_group)
ward_egfp <- data.frame(group = factor(group_names, levels = group_names), egfp_positive = egfp_per_group, share = egfp_share)
print(ward_egfp)

# Figure 4f radar: share of all EGFP+ neurons in each Ward group (the first two rows set the 100% and 0% axis limits)
radar_rows <- rbind(rep(100, 9), rep(0, 9), egfp_share)
radar_table <- as.data.frame(radar_rows)
colnames(radar_table) <- group_names
share_text <- format(round(egfp_share, 1), nsmall = 1)
radar_legend <- paste0(group_names, ": ", share_text, "% (n = ", egfp_per_group, ")")
par(mar = c(1, 1, 4, 12))
radarchart(radar_table, axistype = 0, pcol = rgb(0.2, 0.5, 0.9, 0.9), pfcol = rgb(0.2, 0.5, 0.9, 0.3), plwd = 2.5, cglcol = "grey", cglty = 1,
           title = "Share of EGFP+ Neurons in Each Ward Group")
legend(x = 1.4, y = 0.9, legend = radar_legend, title = "% of all EGFP+", bty = "n", cex = 0.85, xpd = NA)

# Figure 4f bars: share of all EGFP+ neurons in each Ward group (number of EGFP+ neurons on each bar)
p_4f_bars <- ggplot(ward_egfp, aes(x = group, y = share, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = egfp_positive), vjust = -0.4, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Ward group", y = "% of all EGFP+ neurons", title = "Share of EGFP+ Neurons in Each Ward Group") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  bold_title
print(p_4f_bars)

# =============================================================================
# 14. EGFP+ vs EGFP- volcano plots (neurons alone in their droplet; EGFP genes out)
# =============================================================================
# Two comparisons: all EGFP+ vs all EGFP- putative SCN neurons; Group 1 EGFP+ vs every EGFP- neuron (EGFP+ neurons of the
# other groups left out). Wilcoxon (FindMarkers), logfc.threshold = 0, min.pct = 0.25; genes called at adjusted p < 0.05
# and |log2 fold change| > 0.25. Vipr2 and Prokr2 are also tested with min.pct = 0 so they always have a point
alone_cells <- colnames(scn_seurat)[scn_seurat$no_other_cell]
alone_seurat <- subset(scn_seurat, cells = alone_cells)
genes_without_egfp <- setdiff(rownames(scn_seurat), egfp_rows)
label_genes <- c("Vipr2", "Prokr2")

# EGFP+ or EGFP- for each comparison; "excluded" neurons are not tested
all_status <- rep("EGFP-", ncol(alone_seurat))
all_status[alone_seurat$gfp_positive] <- "EGFP+"
alone_seurat$all_status <- all_status
g1_status <- rep("EGFP-", ncol(alone_seurat))
g1_status[alone_seurat$gfp_positive] <- "excluded"
g1_status[alone_seurat$gfp_positive & alone_seurat$group_label == "Group 1"] <- "EGFP+"
alone_seurat$g1_status <- g1_status
print(table(alone_seurat$all_status))
print(table(alone_seurat$g1_status))

de_all <- FindMarkers(alone_seurat, assay = "RNA", group.by = "all_status", ident.1 = "EGFP+", ident.2 = "EGFP-", features = genes_without_egfp, logfc.threshold = 0, min.pct = 0.25)
de_g1 <- FindMarkers(alone_seurat, assay = "RNA", group.by = "g1_status", ident.1 = "EGFP+", ident.2 = "EGFP-", features = genes_without_egfp, logfc.threshold = 0, min.pct = 0.25)
de_all_label_genes <- FindMarkers(alone_seurat, assay = "RNA", group.by = "all_status", ident.1 = "EGFP+", ident.2 = "EGFP-", features = label_genes, logfc.threshold = 0, min.pct = 0)
de_g1_label_genes <- FindMarkers(alone_seurat, assay = "RNA", group.by = "g1_status", ident.1 = "EGFP+", ident.2 = "EGFP-", features = label_genes, logfc.threshold = 0, min.pct = 0)

# The Vipr2 and Prokr2 rows from the min.pct = 0 test replace any row from the main test
main_rows_all <- !(rownames(de_all) %in% label_genes)
de_all <- rbind(de_all[main_rows_all, ], de_all_label_genes)
main_rows_g1 <- !(rownames(de_g1) %in% label_genes)
de_g1 <- rbind(de_g1[main_rows_g1, ], de_g1_label_genes)
# Gene names (the row names) into a gene column
de_all$gene <- rownames(de_all)
de_g1$gene <- rownames(de_g1)

# Called genes: detected in at least 25% of either group, adjusted p < 0.05 and |log2 fold change| > 0.25
called_all <- pmax(de_all$pct.1, de_all$pct.2) >= 0.25 & de_all$p_val_adj < 0.05
de_all$direction <- "Not significant"
de_all$direction[called_all & de_all$avg_log2FC > 0.25] <- "Up in EGFP+"
de_all$direction[called_all & de_all$avg_log2FC < -0.25] <- "Down in EGFP+"
called_g1 <- pmax(de_g1$pct.1, de_g1$pct.2) >= 0.25 & de_g1$p_val_adj < 0.05
de_g1$direction <- "Not significant"
de_g1$direction[called_g1 & de_g1$avg_log2FC > 0.25] <- "Up in EGFP+"
de_g1$direction[called_g1 & de_g1$avg_log2FC < -0.25] <- "Down in EGFP+"
# -log10 adjusted p for the y axis (1e-300 keeps p = 0 finite)
de_all$neg_log10_padj <- -log10(de_all$p_val_adj + 1e-300)
de_g1$neg_log10_padj <- -log10(de_g1$p_val_adj + 1e-300)

# Labels: the 10 called genes with the smallest adjusted p on each side, plus Vipr2 and Prokr2
up_rows_all <- de_all[de_all$direction == "Up in EGFP+", ]
down_rows_all <- de_all[de_all$direction == "Down in EGFP+", ]
top_up_all <- head(up_rows_all$gene[order(up_rows_all$p_val_adj)], 10)
top_down_all <- head(down_rows_all$gene[order(down_rows_all$p_val_adj)], 10)
show_all <- de_all$gene %in% c(top_up_all, top_down_all, label_genes)
de_all$label <- ""
de_all$label[show_all] <- de_all$gene[show_all]
up_rows_g1 <- de_g1[de_g1$direction == "Up in EGFP+", ]
down_rows_g1 <- de_g1[de_g1$direction == "Down in EGFP+", ]
top_up_g1 <- head(up_rows_g1$gene[order(up_rows_g1$p_val_adj)], 10)
top_down_g1 <- head(down_rows_g1$gene[order(down_rows_g1$p_val_adj)], 10)
show_g1 <- de_g1$gene %in% c(top_up_g1, top_down_g1, label_genes)
de_g1$label <- ""
de_g1$label[show_g1] <- de_g1$gene[show_g1]

# Genes not called drawn small and first, called genes drawn larger on top
not_called_all <- de_all[de_all$direction == "Not significant", ]
called_rows_all <- de_all[de_all$direction != "Not significant", ]
not_called_all$point_size <- 0.7
called_rows_all$point_size <- 1.2
de_all <- rbind(not_called_all, called_rows_all)
not_called_g1 <- de_g1[de_g1$direction == "Not significant", ]
called_rows_g1 <- de_g1[de_g1$direction != "Not significant", ]
not_called_g1$point_size <- 0.7
called_rows_g1$point_size <- 1.2
de_g1 <- rbind(not_called_g1, called_rows_g1)

# Number of called genes per comparison
up_all <- sum(de_all$direction == "Up in EGFP+")
down_all <- sum(de_all$direction == "Down in EGFP+")
up_g1 <- sum(de_g1$direction == "Up in EGFP+")
down_g1 <- sum(de_g1$direction == "Down in EGFP+")
print(c(up_all = up_all, down_all = down_all, up_g1 = up_g1, down_g1 = down_g1))

# Symmetric x axis for each volcano
x_limit_all <- 1.05 * max(abs(de_all$avg_log2FC))
x_limit_g1 <- 1.05 * max(abs(de_g1$avg_log2FC))

# Volcanoes: dashed lines at |log2 fold change| = 0.25 and adjusted p = 0.05, called-gene counts in the top corners
p_volcano_all <- ggplot(de_all, aes(avg_log2FC, neg_log10_padj)) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = direction, size = point_size)) +
  geom_text_repel(aes(label = label), fontface = "italic", size = 3, max.overlaps = Inf,
                  min.segment.length = 0, segment.colour = "grey45", segment.size = 0.3, box.padding = 0.4, bg.colour = "white", bg.r = 0.15, seed = 42) +
  scale_colour_manual(values = c("Up in EGFP+" = "#B2182B", "Down in EGFP+" = "#2166AC", "Not significant" = "grey82"), guide = "none") +
  scale_size_identity() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(x = expression(log[2]~"fold change (EGFP+ / EGFP\u2212)"), y = expression(-log[10]~"adjusted p"), title = "EGFP+ vs EGFP\u2212 Putative SCN Neurons") +
  theme_classic() +
  bold_title +
  annotate("text", x = x_limit_all, y = Inf, label = paste0(up_all, " genes higher\nin EGFP+ neurons"), hjust = 1, vjust = 1.3, colour = "#B2182B", size = 3.6) +
  annotate("text", x = -x_limit_all, y = Inf, label = paste0(down_all, " genes lower\nin EGFP+ neurons"), hjust = 0, vjust = 1.3, colour = "#2166AC", size = 3.6) +
  scale_x_continuous(limits = c(-x_limit_all, x_limit_all))
print(p_volcano_all)
p_volcano_g1 <- ggplot(de_g1, aes(avg_log2FC, neg_log10_padj)) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = direction, size = point_size)) +
  geom_text_repel(aes(label = label), fontface = "italic", size = 3, max.overlaps = Inf,
                  min.segment.length = 0, segment.colour = "grey45", segment.size = 0.3, box.padding = 0.4, bg.colour = "white", bg.r = 0.15, seed = 42) +
  scale_colour_manual(values = c("Up in EGFP+" = "#B2182B", "Down in EGFP+" = "#2166AC", "Not significant" = "grey82"), guide = "none") +
  scale_size_identity() +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(x = expression(log[2]~"fold change (EGFP+ / EGFP\u2212)"), y = expression(-log[10]~"adjusted p"), title = "Group 1 EGFP+ vs All EGFP\u2212 Putative SCN Neurons") +
  theme_classic() +
  bold_title +
  annotate("text", x = x_limit_g1, y = Inf, label = paste0(up_g1, " genes higher\nin Group 1 EGFP+"), hjust = 1, vjust = 1.3, colour = "#B2182B", size = 3.6) +
  annotate("text", x = -x_limit_g1, y = Inf, label = paste0(down_g1, " genes lower\nin Group 1 EGFP+"), hjust = 0, vjust = 1.3, colour = "#2166AC", size = 3.6) +
  scale_x_continuous(limits = c(-x_limit_g1, x_limit_g1))
print(p_volcano_g1)

# =============================================================================
# 15. EGFP+ neurons per Ward group at 1+, 2+ and 3+ EGFP molecules
# =============================================================================
# Neurons alone in their droplet, and EGFP+ neurons at 1+, 2+ and 3+ EGFP molecules, counted per Ward group. The three
# EGFP+ sets are nested: every neuron with 3+ EGFP molecules is also in the 2+ and 1+ sets
table_all <- table(alone_seurat$group_label)
cells_all <- as.numeric(table_all)
table_1umi <- table(alone_seurat$group_label[alone_seurat$egfp_droplet_umi >= 1])
cells_1umi <- as.numeric(table_1umi)
table_2umi <- table(alone_seurat$group_label[alone_seurat$egfp_droplet_umi >= 2])
cells_2umi <- as.numeric(table_2umi)
table_3umi <- table(alone_seurat$group_label[alone_seurat$egfp_droplet_umi >= 3])
cells_3umi <- as.numeric(table_3umi)
n_alone_text <- format(sum(cells_all), big.mark = ",")
label_set_all <- paste0("All neurons with no other cell in the droplet (n = ", n_alone_text, ")")
label_set_1umi <- paste0("EGFP+, at least 1 EGFP molecule (n = ", sum(cells_1umi), ")")
label_set_2umi <- paste0("EGFP+, at least 2 EGFP molecules (n = ", sum(cells_2umi), ")")
label_set_3umi <- paste0("EGFP+, at least 3 EGFP molecules (n = ", sum(cells_3umi), ")")
# One table per set (percent of the set in each group, so each set sums to 100), stacked
set_pct_all <- 100 * cells_all / sum(cells_all)
bars_all <- data.frame(group = group_names, cell_set = label_set_all, cells = cells_all, pct = set_pct_all)
set_pct_1umi <- 100 * cells_1umi / sum(cells_1umi)
bars_1umi <- data.frame(group = group_names, cell_set = label_set_1umi, cells = cells_1umi, pct = set_pct_1umi)
set_pct_2umi <- 100 * cells_2umi / sum(cells_2umi)
bars_2umi <- data.frame(group = group_names, cell_set = label_set_2umi, cells = cells_2umi, pct = set_pct_2umi)
set_pct_3umi <- 100 * cells_3umi / sum(cells_3umi)
bars_3umi <- data.frame(group = group_names, cell_set = label_set_3umi, cells = cells_3umi, pct = set_pct_3umi)
cutoff_bars <- rbind(bars_all, bars_1umi, bars_2umi, bars_3umi)
cutoff_bars$group <- factor(cutoff_bars$group, levels = group_names)
cutoff_bars$cell_set <- factor(cutoff_bars$cell_set, levels = c(label_set_all, label_set_1umi, label_set_2umi, label_set_3umi))
print(cutoff_bars)

p_cutoffs <- ggplot(cutoff_bars, aes(x = group, y = pct, fill = cell_set)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = cells), position = position_dodge(width = 0.85), vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = c("grey80", "#a1d99b", "#41ab5d", "#006d2c"), name = NULL) +
  labs(x = "Ward group", y = "% of the set in each group (number on bar = neurons)", title = "EGFP+ Neurons per Ward Group at Three EGFP Cutoffs") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "vertical", axis.text.x = element_text(angle = 45, hjust = 1)) +
  bold_title
print(p_cutoffs)

# =============================================================================
# 16. Morris 2021 and Wen 2020 gene-list scores on the putative SCN UMAP
# =============================================================================
# Morris CT7.5 day-group lists A-H and the putative SCN and extra-SCN cell-class lists (GEO GSE167927), each sheet in full
morris_path <- file.path("reference_data", "GSE167927_CT7.5_SCN_Cell_Group_Genes.xlsx")
morris_classes_path <- file.path("reference_data", "GSE167927_CT7.5_Cell_Groups_All_Sequenced_Cells.xlsx")
morris_a <- read_excel(morris_path, sheet = "Genes Day Group A")
morris_b <- read_excel(morris_path, sheet = "Genes Day Group B")
morris_c <- read_excel(morris_path, sheet = "Genes Day Group C")
morris_d <- read_excel(morris_path, sheet = "Genes Day Group D")
morris_e <- read_excel(morris_path, sheet = "Genes Day Group E")
morris_f <- read_excel(morris_path, sheet = "Genes Day Group F")
morris_g <- read_excel(morris_path, sheet = "Genes Day Group G")
morris_h <- read_excel(morris_path, sheet = "Genes Day Group H")
morris_scn_cells <- read_excel(morris_classes_path, sheet = "Day Putative SCN cells")
morris_extra_scn <- read_excel(morris_classes_path, sheet = "Day Extra-SCN Hypothalamus cell")
# One table per list: Ensembl ID and list name
morris_a_genes <- data.frame(ensembl = morris_a$FeatureID, subtype = "Morris A (Gem/Ppp1r17/Anxa5)")
morris_b_genes <- data.frame(ensembl = morris_b$FeatureID, subtype = "Morris B (Vip/C1ql3/Calb2)")
morris_c_genes <- data.frame(ensembl = morris_c$FeatureID, subtype = "Morris C (Grp/Synpr/Cd24a)")
morris_d_genes <- data.frame(ensembl = morris_d$FeatureID, subtype = "Morris D (Prok2/Avp/Nms)")
morris_e_genes <- data.frame(ensembl = morris_e$FeatureID, subtype = "Morris E (Cck/Prokr2/Nsmf)")
morris_f_genes <- data.frame(ensembl = morris_f$FeatureID, subtype = "Morris F (Mef2c/Peg10/Pcdh11x)")
morris_g_genes <- data.frame(ensembl = morris_g$FeatureID, subtype = "Morris G (Gadd45a/Hspa1a/Cck)")
morris_h_genes <- data.frame(ensembl = morris_h$FeatureID, subtype = "Morris H (Vipr2/Lbh/Gas7)")
morris_scn_genes <- data.frame(ensembl = morris_scn_cells$FeatureID, subtype = "Morris putative SCN cells")
morris_extra_genes <- data.frame(ensembl = morris_extra_scn$FeatureID, subtype = "Morris extra-SCN hypothalamic cells")
morris_markers <- rbind(morris_a_genes, morris_b_genes, morris_c_genes, morris_d_genes, morris_e_genes,
                        morris_f_genes, morris_g_genes, morris_h_genes, morris_scn_genes, morris_extra_genes)
# Morris Ensembl IDs to gene names on the Flex probe panel (gene names as in the Seurat object: duplicates made unique, _ to -)
flex_features <- read.delim(file.path("Jones Lab data", "features.tsv.gz"), header = FALSE, col.names = c("ensembl", "symbol", "type"))
flex_unique_symbol <- make.unique(flex_features$symbol)
flex_features$object_name <- gsub("_", "-", flex_unique_symbol)
flex_feature_row <- match(morris_markers$ensembl, flex_features$ensembl)
morris_markers$gene <- flex_features$object_name[flex_feature_row]

# Wen 2020 Supplementary Table 1 SCN subtype markers with avg_logFC > 0 (matched to the Flex panel by symbol)
wen_table <- read_excel(file.path("reference_data", "Wen2020_NatNeurosci_Supplementary_Tables_1-6.xlsx"), sheet = "Supplementary Table 1", range = "U4:AA300")
wen_table <- wen_table[!is.na(wen_table$gene), ]
wen_table <- wen_table[wen_table$avg_logFC > 0, ]
wen_markers <- data.frame(gene = wen_table$gene, subtype = paste0("Wen ", wen_table$cluster))

# A gene enters the score only if it is in the object (genes detected in 3+ barcodes; EGFP rows excluded)
wen_levels <- paste0("Wen ", c("Avp/Nms", "Vip/Nms", "Grp/Vip", "Cck/C1ql3", "Cck/Bdnf"))
subtype_levels <- c(unique(morris_markers$subtype), wen_levels)
published_markers <- rbind(morris_markers[, c("gene", "subtype")], wen_markers[, c("gene", "subtype")])
published_markers$subtype <- factor(published_markers$subtype, levels = subtype_levels)
published_on_panel <- published_markers[published_markers$gene %in% genes_without_egfp, ]

# Genes of each list that entered the score
print(table(published_on_panel$subtype))

# AddModuleScore per list (Seurat defaults), z-scored across neurons
subtype_lists <- split(published_on_panel$gene, published_on_panel$subtype)
scn_seurat <- AddModuleScore(scn_seurat, features = subtype_lists, name = "subtype_score", assay = "RNA")
score_columns <- paste0("subtype_score", 1:length(subtype_lists))
subtype_scores <- as.matrix(scn_seurat@meta.data[, score_columns])
colnames(subtype_scores) <- names(subtype_lists)
subtype_z <- scale(subtype_scores)

# One row per neuron per list (UMAP coordinates and z-score), lowest scores drawn first
score_table <- data.frame(umap_1 = scn_umap$umap_1, umap_2 = scn_umap$umap_2, subtype_z, check.names = FALSE)
score_long <- pivot_longer(score_table, cols = -c(umap_1, umap_2), names_to = "subtype", values_to = "z")
score_long$subtype <- factor(score_long$subtype, levels = subtype_levels)
score_long <- score_long[order(score_long$z), ]

# Panel titles: list name on the first line, its marker genes on the second
subtype_panel_titles <- c("Morris A (Gem/Ppp1r17/Anxa5)" = "Morris group A\nGem/Ppp1r17/Anxa5",
                          "Morris B (Vip/C1ql3/Calb2)" = "Morris group B\nVip/C1ql3/Calb2",
                          "Morris C (Grp/Synpr/Cd24a)" = "Morris group C\nGrp/Synpr/Cd24a",
                          "Morris D (Prok2/Avp/Nms)" = "Morris group D\nProk2/Avp/Nms",
                          "Morris E (Cck/Prokr2/Nsmf)" = "Morris group E\nCck/Prokr2/Nsmf",
                          "Morris F (Mef2c/Peg10/Pcdh11x)" = "Morris group F\nMef2c/Peg10/Pcdh11x",
                          "Morris G (Gadd45a/Hspa1a/Cck)" = "Morris group G\nGadd45a/Hspa1a/Cck",
                          "Morris H (Vipr2/Lbh/Gas7)" = "Morris group H\nVipr2/Lbh/Gas7",
                          "Morris putative SCN cells" = "Morris\nputative SCN cells",
                          "Morris extra-SCN hypothalamic cells" = "Morris\nextra-SCN hypothalamic cells",
                          "Wen Avp/Nms" = "Wen\nAvp/Nms", "Wen Vip/Nms" = "Wen\nVip/Nms", "Wen Grp/Vip" = "Wen\nGrp/Vip",
                          "Wen Cck/C1ql3" = "Wen\nCck/C1ql3", "Wen Cck/Bdnf" = "Wen\nCck/Bdnf")
egfp_legend_label <- paste0("EGFP+ neurons (", nrow(egfp_scn_umap), ")")

# One UMAP per list (colour capped at z = 3), EGFP+ neurons in black
p_gene_lists <- ggplot(score_long, aes(x = umap_1, y = umap_2, colour = z)) +
  geom_point_rast(size = 0.15, raster.dpi = 300, dev = "ragg") +
  geom_point(data = egfp_scn_umap, aes(x = umap_1, y = umap_2, shape = legend), inherit.aes = FALSE, colour = "black", size = 0.6) +
  scale_shape_manual(values = c("EGFP+ neurons" = 16), labels = egfp_legend_label, name = NULL) +
  labs(title = "Published SCN Gene-List Scores on Putative SCN Neurons") +
  guides(shape = guide_legend(override.aes = list(size = 2.5))) +
  scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0, limits = c(-3, 3), oob = scales::squish, name = "Score (z)") +
  facet_wrap(~ subtype, ncol = 5, labeller = labeller(subtype = subtype_panel_titles)) +
  coord_fixed() +
  theme_void() +
  theme(strip.text = element_text(size = 9, margin = margin(2, 0, 4, 0)), legend.position = "bottom", legend.key.width = unit(8, "mm"),
        legend.key.height = unit(2.5, "mm"), panel.spacing = unit(4, "mm"), plot.margin = margin(8, 8, 8, 8)) +
  bold_title
print(p_gene_lists)

# =============================================================================
# 17. Print each figure again and save it as PDF
# =============================================================================
print(p_3a)
ggsave(file.path("whole droplet paper figures", "Supplementary Figure 3", "Supplementary_Figure_3a_SingleR_Cell_Types.pdf"), p_3a, width = 9, height = 6.5, device = grDevices::cairo_pdf)
print(p_3b)
ggsave(file.path("whole droplet paper figures", "Supplementary Figure 3", "Supplementary_Figure_3b_Neuron_Clusters_Screened_for_Removal.pdf"), p_3b, width = 8, height = 6.5, device = grDevices::cairo_pdf)
print(p_4b)
ggsave(file.path("whole droplet paper figures", "Figure 4", "Figure_4b_SCN_Genes_UMAP.pdf"), p_4b, width = 17, height = 12, device = grDevices::cairo_pdf)
print(p_4c)
ggsave(file.path("whole droplet paper figures", "Figure 4", "Figure_4c_Neurons_with_Detectable_RNA.pdf"), p_4c, width = 9, height = 5, device = grDevices::cairo_pdf)
# Ward tree: drawn on screen, then drawn again into its PDF file
par(mar = c(6, 4, 4, 2))
plot(ward_dend, main = ward_title, xlab = "", ylab = "Height (Ward distance)")
mtext("Putative SCN clusters (resolution 0.5)", side = 1, line = 4.5, cex = 0.8)
segments(bracket_left, bracket_y, bracket_right, bracket_y, lwd = 1.5, xpd = TRUE)
segments(bracket_left, bracket_y, bracket_left, tick_y, lwd = 1.5, xpd = TRUE)
segments(bracket_right, bracket_y, bracket_right, tick_y, lwd = 1.5, xpd = TRUE)
text(group_middle, bracket_y, labels = paste0("G", names(group_left)), pos = 1, cex = 0.8, xpd = TRUE)
pdf(file.path("whole droplet paper figures", "Figure 4", "Figure_4d_Ward_Dendrogram.pdf"), width = 7, height = 5.5)
par(mar = c(6, 4, 4, 2))
plot(ward_dend, main = ward_title, xlab = "", ylab = "Height (Ward distance)")
mtext("Putative SCN clusters (resolution 0.5)", side = 1, line = 4.5, cex = 0.8)
segments(bracket_left, bracket_y, bracket_right, bracket_y, lwd = 1.5, xpd = TRUE)
segments(bracket_left, bracket_y, bracket_left, tick_y, lwd = 1.5, xpd = TRUE)
segments(bracket_right, bracket_y, bracket_right, tick_y, lwd = 1.5, xpd = TRUE)
text(group_middle, bracket_y, labels = paste0("G", names(group_left)), pos = 1, cex = 0.8, xpd = TRUE)
dev.off()
print(p_4d)
ggsave(file.path("whole droplet paper figures", "Figure 4", "Figure_4d_Ward_Groups_UMAP.pdf"), p_4d, width = 8, height = 6, device = grDevices::cairo_pdf)
print(p_4e)
ggsave(file.path("whole droplet paper figures", "Figure 4", "Figure_4e_Ward_Group_Markers.pdf"), p_4e, width = 12, height = 10, device = grDevices::cairo_pdf)
# Radar: drawn on screen, then drawn again into its PDF file
par(mar = c(1, 1, 4, 12))
radarchart(radar_table, axistype = 0, pcol = rgb(0.2, 0.5, 0.9, 0.9), pfcol = rgb(0.2, 0.5, 0.9, 0.3), plwd = 2.5, cglcol = "grey", cglty = 1,
           title = "Share of EGFP+ Neurons in Each Ward Group")
legend(x = 1.4, y = 0.9, legend = radar_legend, title = "% of all EGFP+", bty = "n", cex = 0.85, xpd = NA)
pdf(file.path("whole droplet paper figures", "Figure 4", "Figure_4f_EGFP_Share_per_Ward_Group.pdf"), width = 10, height = 6.5)
par(mar = c(1, 1, 4, 12))
radarchart(radar_table, axistype = 0, pcol = rgb(0.2, 0.5, 0.9, 0.9), pfcol = rgb(0.2, 0.5, 0.9, 0.3), plwd = 2.5, cglcol = "grey", cglty = 1,
           title = "Share of EGFP+ Neurons in Each Ward Group")
legend(x = 1.4, y = 0.9, legend = radar_legend, title = "% of all EGFP+", bty = "n", cex = 0.85, xpd = NA)
dev.off()
print(p_4f_bars)
ggsave(file.path("whole droplet paper figures", "Figure 4", "Figure_4f_EGFP_Share_per_Ward_Group_Bars.pdf"), p_4f_bars, width = 8, height = 5, device = grDevices::cairo_pdf)
print(p_volcano_all)
ggsave(file.path("whole droplet paper figures", "Volcano Plots", "Volcano_EGFP_All_Putative_SCN_Neurons.pdf"), p_volcano_all, width = 8, height = 7.5, device = grDevices::cairo_pdf)
print(p_volcano_g1)
ggsave(file.path("whole droplet paper figures", "Volcano Plots", "Volcano_EGFP_Group_1_vs_All_EGFP_Negative.pdf"), p_volcano_g1, width = 8, height = 7.5, device = grDevices::cairo_pdf)
print(p_cutoffs)
ggsave(file.path("whole droplet paper figures", "EGFP per Ward Group", "EGFP_Neurons_per_Ward_Group_at_Three_Cutoffs.pdf"), p_cutoffs, width = 11, height = 6.5, device = grDevices::cairo_pdf)
print(p_gene_lists)
ggsave(file.path("whole droplet paper figures", "Morris and Wen Gene List Scores", "Morris_Wen_Gene_List_Scores_UMAP.pdf"), p_gene_lists, width = 11, height = 6.3, device = grDevices::cairo_pdf)

sessionInfo()
