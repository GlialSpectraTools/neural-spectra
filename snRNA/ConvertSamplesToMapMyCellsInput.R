library(Seurat)
library(rliger)
library(scCustomize)
library(qs)
library(anndata)
library(reticulate)
reticulate::py_install("anndata")

##############
# Read each sample and convert to MapMyCells input
##

# Define paths
input_folder <- "Seurat_Objects"
output_folder <- "MapMyCells/Input_sample_level"

# List all relevant .rds files (starting with "SF")
rds_files <- list.files(path = input_folder, pattern = "^SF.*\\.rds$", full.names = FALSE)

# Loop through each file and convert to h5ad (compataply with mapmycells)
for (file in rds_files) {
  # Read RDS file
  obj <- readRDS(file.path(input_folder, file))
  
  # Set default assay to RNA
  DefaultAssay(obj) <- "RNA"
  
  # Define output file name (replace .rds with _RNA.h5ad)
  output_file_name <- paste0(tools::file_path_sans_ext(file), "_RNA.h5ad")
  
  # Convert to .h5ad
  as.anndata(x = obj, file_path = output_folder, file_name = output_file_name)
  
  # Print progress
  print(paste("Converted:", file, "to", output_file_name))
}


#### Then use MTG and hierachal mapping in mapmycellsgui