library(Seurat)
library(dplyr)
library(stringr)

# Load Seurat object
obj <- readRDS("oneoverf_combined_submap_cells_snc_sct_complete.rds")

# Base folder for MapMyCells outputs
base_dir <- "MapMyCells/Output_Sample_level"
mapping_files <- list.files(path = base_dir, pattern = "*.csv$", recursive = TRUE, full.names = TRUE)

all_mappings <- list()

for (file in mapping_files) {
  # Extract sample ID
  sample_id <- str_extract(basename(file), "^SF[^_]+_[^_]+")
  
  
  # Try reading and processing the file
  mapping <- tryCatch({
    df <- read.csv(file, comment.char = "#")
    if (!"cell_id" %in% colnames(df)) stop("Missing cell_id")
    df$cell_id <- paste0(sample_id, "_", str_trim(df$cell_id))
    rownames(df) <- df$cell_id
    df[, c("class_name", "subclass_name"), drop = FALSE]
  }, error = function(e) {
    message(paste("Failed to load:", file))
    NULL
  })
  
  if (!is.null(mapping)) {
    all_mappings[[sample_id]] <- mapping
  }
}


# Combine and annotate
combined_mapping <- bind_rows(all_mappings)
common_cells <- intersect(rownames(combined_mapping), colnames(obj))
combined_mapping <- combined_mapping[common_cells, , drop = FALSE]
obj <- AddMetaData(obj, metadata = combined_mapping)

# Confirm result
cat("Annotated cells:", length(common_cells), "of", ncol(obj), "\n")


# Keep consistent colors
my_colors_2 <- c(
  "L2/3 IT"         = "#FEC89A",
  "L4 IT"           = "#FCD5CE",
  "L5 IT"           = "#FFDAB9",
  "L6 IT"           = "#F49595",
  "L6 IT Car3"      = "#F1CEBE",
  "L5 ET"           = "#F3D17C",
  "L6 CT"           = "#FCEFB4",
  "L5/6 NP"         = "#F9ED85",
  "L6b"             = "#FFCAE9",
  "Chandelier"      = "#A9E4DE",
  "Sst Chodl"       = "#83B8C6",
  "Lamp5"           = "#D0E444",    
  "Lamp5 Lhx6"      = "#D0EDEF",
  "Pax6"            = "#E1D3F8",
  "Vip"             = "#B0B5ED", 
  "Sncg"            = "#D6DFE8",
  "Pvalb"           = "#8FC0E3",    
  "Sst"             = "#7FB3D5",
  
  # Newly added distinct colors for unspecified clusters
  "OPC"             = "#8AC926",    # Vibrant green
  "Oligodendrocyte" = "#FF7F51",    # Bright orange-red
  "Microglia-PVM"   = "#6A4C93",    # Deep purple
  "Astrocyte"       = "#1982C4",    # Rich blue
  "VLMC"            = "#FF595E",    # Vibrant red-pink
  "Endothelial"     = "#52B69A"     # Teal green
)


# Visualize UMAP with new annotations
DimPlot(obj, label="T")
DimPlot(obj, reduction = "umap", group.by = "subclass_name", label = TRUE)
plot <- DimPlot(obj, label=F, group.by = "subclass_name", cols=my_colors_2)  +NoAxes()

#manual anottations using findallmarkers
cluster_ids <- c(
  "0" = "Tumor",
  "1" = "Oligodendrocyte",
  "2" = "Endothelial",
  "3" = "Oligodendrocyte",
  "4" = "Astrocyte",
  "5" = "Oligodendrocyte",
  "6" = "Neurons",
  "7" = "Myeloid/microglial",
  "8" = "Astrocyte",
  "9" = "Tumor",
  "10" = "Neurons",
  "11" = "Neurons",
  "12" = "Tumor",
  "13" = "Neurons",
  "14" = "Neurons",
  "15" = "Neurons",
  "16" = "Neurons",
  "17" = "Neurons",
  "18" = "Oligodendrocyte",
  "19" = "Tumor",
  "20" = "Oligodendrocyte",
  "21" = "Mesenchymal",
  "22" = "Erythroid",
  "23" = "Oligodendrocyte",
  "24" = "Neurons",
  "25" = "Tumor",
  "26" = "Myeloid",
  "27" = "Neurons",
  "28" = "Tumor",
  "29" = "Oligodendrocyte",
  "30" = "Tumor",
  "31" = "Neurons",
  "32" = "Neurons",
  "33" = "Astrocyte",
  "34" = "Tumor"
)

# Assign annotations to a new metadata column 'manual_annot'
# Convert cluster factor levels explicitly
cluster_labels <- as.character(obj$seurat_clusters)

# Map annotations directly without names
manual_annotations <- unname(cluster_ids[cluster_labels])

# Add to metadata directly
obj$manual_annot <- manual_annotations

my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
               "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

plot2 <- DimPlot(obj, label=F, group.by = "manual_annot", cols=my_colors)  +NoAxes()





#refined ids
# First copy existing annotations
obj$refined_anot <- obj$subclass_name

# Define your tumor cluster IDs (as numeric)
tumor_clusters <- c(0, 9, 12, 19, 25, 28, 30, 34)

# Update the refined_anot for cells in tumor clusters
obj$refined_anot[obj$seurat_clusters %in% tumor_clusters] <- "Tumor"


my_colors_3 <- c(
  "L2/3 IT"         = "#FEC89A",
  "L4 IT"           = "#FCD5CE",
  "L5 IT"           = "#FFDAB9",
  "L6 IT"           = "#F49595",
  "L6 IT Car3"      = "#F1CEBE",
  "L5 ET"           = "#F3D17C",
  "L6 CT"           = "#FCEFB4",
  "L5/6 NP"         = "#F9ED85",
  "L6b"             = "#FFCAE9",
  "Chandelier"      = "#A9E4DE",
  "Sst Chodl"       = "#83B8C6",
  "Lamp5"           = "#D0E444",    
  "Lamp5 Lhx6"      = "#D0EDEF",
  "Pax6"            = "#E1D3F8",
  "Vip"             = "#B0B5ED", 
  "Sncg"            = "#D6DFE8",
  "Pvalb"           = "#8FC0E3",    
  "Sst"             = "#7FB3D5",
  
  # Previously added colors
  "OPC"             = "#8AC926",    # Vibrant green
  "Oligodendrocyte" = "#FF7F51",    # Bright orange-red
  "Microglia-PVM"   = "#6A4C93",    # Deep purple
  "Astrocyte"       = "#1982C4",    # Rich blue
  "VLMC"            = "#FF595E",    # Vibrant red-pink
  "Endothelial"     = "#52B69A",    # Teal green
  
  # Add Tumor explicitly
  "Tumor"           = "#BEBEBE"     # Gray
)

plot3 <- DimPlot(obj, label = FALSE, group.by = "refined_anot", cols = my_colors_3) + NoAxes()
plot3



library(Seurat)
library(ggplot2)
library(viridis)

# Define your custom magma color palette explicitly
magma_colors <- c("#000000", magma(256))

# Your marker genes stored clearly in a vector
General_Markers <- c("EGFR", "GFAP", "IDH1", "PECAM1", "RBFOX3",
                     "SOX6", "SOX9", "SOX10", "TMEM119", "C1QA", 
                     "CD68", "TP53")

# Set your save path clearly (ensure folder exists)
save_path <- "MapMyCells/Marker_Genes/MarkerGenes_Magma/"

# Plot marker genese
for (gene in General_Markers) {
  

  p <- FeaturePlot(obj, 
                   features = gene,
                   order = TRUE, 
                   raster = FALSE,
                   min.cutoff = 0,
                   max.cutoff = 3.5) +
    scale_color_gradientn(colors = magma_colors,
                          limits = c(0, 3.5),
                          breaks = seq(0, 3.5, by = 0.5),
                          oob = scales::squish) +
    NoAxes()
  

  file_name <- paste0(save_path, gene, "_AllCells_UMAP_Magma.pdf")
  
  # Save plot 
  ggsave(filename = file_name,
         plot = p,
         width = 6.93,
         height = 5.3,
         units = "in",
         device = "pdf")
  
  cat("Plot saved successfully for:", gene, "\n")
}


#Neuron Sub pops
#Excit
# excitatory neuron markers
excitatory_neuron_markers <- c("SLC17A6", "CUX2", "RORB", "IL1RAPL2", 
                               "ZNF804B", "FEZF2", "NWD2", "TSHZ2", "SYT6")

# Explicitly define your new save path (ensure it exists)
save_path <- "MapMyCells/Marker_Genes/MarkerGenes_ExcitNeurons_Magma/"

# plotting excit markers
for (gene in excitatory_neuron_markers) {
  

  p <- FeaturePlot(obj, 
                   features = gene,
                   order = TRUE, 
                   raster = FALSE,
                   min.cutoff = 0,
                   max.cutoff = 3.5) +
    scale_color_gradientn(colors = magma_colors,
                          limits = c(0, 3.5),
                          breaks = seq(0, 3.5, by = 0.5),
                          oob = scales::squish) +
    NoAxes()
  

  file_name <- paste0(save_path, gene, "_AllCells_UMAP_Magma.pdf")
  

  ggsave(filename = file_name,
         plot = p,
         width = 6.93,
         height = 5.3,
         units = "in",
         device = "pdf")
  
  cat("Excitatory marker plot saved successfully for:", gene, "\n")
}



#inhibitory markers
inhib_neuron_markers <- c("GAD1", "MEIS2", "ADARB2", "VIP", "CCK",
                          "LAMP5", "NXPH1", "SST", "PVALB")

# Specify path to inhib plots
save_path <- "MapMyCells/Marker_Genes/MarkerGenes_InhibNeurons_Magma/"

# plotting inhibitory neuron markers
for (gene in inhib_neuron_markers) {
  
  # Create plot 
  p <- FeaturePlot(obj, 
                   features = gene,
                   order = TRUE, 
                   raster = FALSE,
                   min.cutoff = 0,
                   max.cutoff = 3.5) +
    scale_color_gradientn(colors = magma_colors,
                          limits = c(0, 3.5),
                          breaks = seq(0, 3.5, by = 0.5),
                          oob = scales::squish) +
    NoAxes()
  
  # File name
  file_name <- paste0(save_path, gene, "_AllCells_UMAP_Magma.pdf")
  
  # Save plot
  ggsave(filename = file_name,
         plot = p,
         width = 6.93,
         height = 5.3,
         units = "in",
         device = "pdf")
  
  cat("Inhibitory marker plot saved successfully for:", gene, "\n")
}
