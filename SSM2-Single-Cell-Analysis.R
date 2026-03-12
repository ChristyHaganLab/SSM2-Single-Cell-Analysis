# Single Cell RNA-seq Analysis Script for Progesterone and Placebo Treatments
# Date: October 2025
# Authors: Eilidh Chowanec
# NOTES: This code was adapted from Zhang et al. with the help of claude.ai and edited by Jeffery Vahrenkamp
# 
# This script processes 10X Genomics single cell RNA-seq data from 2 treatment groups
# (Progesterone and Placebo) across 2 batches:
# - Batch 1: 1 Placebo sample, 1 Progesterone sample
# - Batch 2: 1 Placebo sample, 2 Progesterone samples
# Total: 2 Placebo samples, 3 Progesterone samples (5 samples total)
#
# Note: This script includes batch correction using Harmony

# ============================================================================
### SETUP INSTRUCTIONS - FOLLOW THESE STEPS BEFORE RUNNING THE SCRIPT
# ============================================================================

### OVERVIEW:
### This guide will help you set up the proper directory structure and configure 
### RStudio for analyzing your 10X Genomics single cell RNA-seq data with 2 
### treatment groups (Progesterone and Placebo) across 2 batches.
###
### EXPERIMENTAL DESIGN:
### - Batch 1: 1 Placebo + 1 Progesterone sample
### - Batch 2: 1 Placebo + 2 Progesterone samples
### - Total: 2 Placebo samples + 3 Progesterone samples = 5 samples

### STEP 1: CREATE THE MAIN PROJECT DIRECTORY
###
### On Windows:
### 1. Navigate to your Desktop
### 2. Right-click on an empty area
### 3. Select "New" → "Folder"
### 4. Name the folder exactly: SSM2sc
###
### On Mac:
### 1. Navigate to your Desktop
### 2. Right-click on an empty area
### 3. Select "New Folder"
### 4. Name the folder exactly: SSM2sc
###
### On Linux:
### 1. Open terminal and navigate to Desktop:
###    cd ~/Desktop
###    mkdir SSM2sc

### STEP 2: CREATE BATCH AND SAMPLE SUBFOLDERS
###
### Inside the SSM2sc folder, you need to create a folder structure that
### reflects your experimental design with 2 batches and 5 samples total.
###
### RECOMMENDED FOLDER STRUCTURE:
### Desktop/
### └── SSM2sc/
###     ├── Batch1_Placebo/
###     ├── Batch1_Progesterone/
###     ├── Batch2_Placebo/
###     ├── Batch2_Progesterone1/
###     └── Batch2_Progesterone2/
###
### Method 1 - Using File Explorer/Finder:
### 1. Double-click to open the SSM2sc folder
### 2. Create 5 new folders with these exact names:
###    - Batch1_Placebo
###    - Batch1_Progesterone
###    - Batch2_Placebo
###    - Batch2_Progesterone1
###    - Batch2_Progesterone2
###
### Method 2 - Using Command Line:
### cd ~/Desktop/SSM2sc
### mkdir "Batch1_Placebo"
### mkdir "Batch1_Progesterone"
### mkdir "Batch2_Placebo"
### mkdir "Batch2_Progesterone1"
### mkdir "Batch2_Progesterone2"
###
### IMPORTANT NAMING NOTES:
### - Use consistent naming that includes both batch and treatment information
### - Avoid spaces in folder names (use underscores instead)
### - Be consistent with capitalization

### STEP 3: DOWNLOAD AND ORGANIZE YOUR 10X DATA FILES
###
### For each sample, you need to download the 10X Genomics output files 
### into their respective folders.
###
### Required Files for Each Sample:
### Each sample folder must contain these 3 files:
###
### 1. barcodes.tsv.gz (or barcodes.tsv)
###    - Contains cell barcode sequences
###    
### 2. features.tsv.gz (or genes.tsv for older 10X format)
###    - Contains gene information (gene IDs and symbols)
###    
### 3. matrix.mtx.gz (or matrix.mtx)
###    - Contains the count matrix (genes × cells)
###
### File Placement Instructions:
### Place the 10X output files in each of your 5 sample folders:
###
### 1. Batch1_Placebo: Place files in Desktop/SSM2sc/Batch1_Placebo/
### 2. Batch1_Progesterone: Place files in Desktop/SSM2sc/Batch1_Progesterone/
### 3. Batch2_Placebo: Place files in Desktop/SSM2sc/Batch2_Placebo/
### 4. Batch2_Progesterone1: Place files in Desktop/SSM2sc/Batch2_Progesterone1/
### 5. Batch2_Progesterone2: Place files in Desktop/SSM2sc/Batch2_Progesterone2/
### 6. Place additional like SSM2sc_filtered.RDS in the SSM2sc folder
###
### Final Directory Structure with Files:
### Desktop/
### └── SSM2sc/
###     ├── Batch1_Placebo/
###     │   ├── barcodes.tsv.gz
###     │   ├── features.tsv.gz
###     │   └── matrix.mtx.gz
###     ├── Batch1_Progesterone/
###     │   ├── barcodes.tsv.gz
###     │   ├── features.tsv.gz
###     │   └── matrix.mtx.gz
###     ├── Batch2_Placebo/
###     │   ├── barcodes.tsv.gz
###     │   ├── features.tsv.gz
###     │   └── matrix.mtx.gz
###     ├── Batch2_Progesterone1/
###     │   ├── barcodes.tsv.gz
###     │   ├── features.tsv.gz
###     │   └── matrix.mtx.gz
###     └── Batch2_Progesterone2/
###         ├── barcodes.tsv.gz
###         ├── features.tsv.gz
###         └── matrix.mtx.gz

### STEP 4: SET UP RSTUDIO WORKING DIRECTORY
###
### Method 1 - Using RStudio Interface (Recommended):
### 1. Open RStudio
### 2. Create a New Project:
###    - File → New Project
###    - Choose "Existing Directory"
###    - Browse to and select your SSM2sc folder
###    - Click "Create Project"
### 3. Verify Setup:
###    - In the RStudio console, type: getwd()
###    - It should return something like: "/path/to/folder/SSM2sc"
###
### Method 2 - Using R Commands:
### 1. Open RStudio
### 2. Set Working Directory:
###    For Windows users: setwd("C:/Users/YourName/Desktop/SSM2sc")
###    For Mac users: setwd("~/Desktop/SSM2sc")
###    For Linux users: setwd("~/Desktop/SSM2sc")
### 3. Verify: getwd() and list.dirs()
###
### Method 3 - Using RStudio's Files Panel:
### 1. Open RStudio
### 2. In bottom-right panel, click "Files" tab
### 3. Navigate to Desktop → SSM2sc
### 4. Click gear icon (More) → "Set As Working Directory"

### STEP 5: VERIFY YOUR SETUP BEFORE RUNNING ANALYSIS
###
### Run these commands in RStudio to verify everything is correct:
###
### getwd()  # Check working directory
### list.dirs(recursive = FALSE)  # Should show your 5 sample folders
###
### # Check if each sample folder contains the required files
### samples <- c("Batch1_Placebo", "Batch1_Progesterone", "Batch2_Placebo", 
###              "Batch2_Progesterone1", "Batch2_Progesterone2")
### for(folder in samples) {
###   cat("Files in", folder, ":\n")
###   print(list.files(folder))
###   cat("\n")
### }
###
### Expected Output for Each Folder:
### Files in Batch1_Placebo :
### [1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"   
###
### Files in Batch1_Progesterone :
### [1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"   
###
### (and so on for all 5 folders...)

### TROUBLESHOOTING:
###
### Common Issues:
### 1. "Directory does not exist" error:
###    - Double-check folder names match exactly what you created
###    - Ensure no extra spaces in folder names
###    - Check capitalization (R is case-sensitive)
###
### 2. "Cannot find files" error:
###    - Verify all 3 files (barcodes, features, matrix) are in each folder
###    - Check file extensions (.tsv.gz or .tsv, .mtx.gz or .mtx)
###    - Make sure files aren't in a subdirectory within the sample folder
###
### 3. Working directory issues:
###    - Run getwd() to check current directory
###    - Use setwd() to change if needed
###    - Make sure you're in the SSM2sc folder, not a subfolder
###
### 4. File permission errors:
###    - Ensure you have read/write permissions for the SSM2sc folder
###    - On Mac/Linux, you may need to adjust permissions using chmod

### READY TO RUN? 
### Once you've completed all setup steps above, you can run the analysis script below!

# ============================================================================
# END OF SETUP INSTRUCTIONS - ANALYSIS SCRIPT BEGINS HERE
# ============================================================================

# Set seed for reproducibility
set.seed(123)

helpMessage = function(){
	args <- commandArgs(FALSE)
	file_arg <- args[grepl("^--file=", args)]
	script_path <- basename(sub("^--file=", "", file_arg))
	message("
Usage:
Rscript ",script_path," [-h 2 3 4 5 6 7]


ARGUMENTS:
default			If no value is given, generates all figures images
-h/--help		Prints this help message
2-7			Generates that figures images. Multiple numbers can be entered.

EXAMPLE:
To generate images for figures 3 and 7

Rscript ",script_path," 3 7

")
q()

}




if ( identical(parent.frame(), .GlobalEnv) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("-h" %in% args || "--help" %in% args){
    helpMessage()
  }
}

loadLibrary <- function(package_name){
	package_name=as.character(substitute(package_name))
	setRepositories(graphics=F,ind=c(1,2,3,4,5,6,7))
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name, dependencies=T);
    }
     suppressPackageStartupMessages(library(package=package_name,character.only=T));
}

# Load required libraries
loadLibrary(R.utils)
loadLibrary(Seurat)
loadLibrary(Matrix)
loadLibrary(stringr)
loadLibrary(ggplot2)
loadLibrary(patchwork)
loadLibrary(plyr)
loadLibrary(dplyr)
loadLibrary(harmony)  # For batch correction

# Increase memory limit (Windows only - Mac/Linux users see below)
if(.Platform$OS.type == "windows") {
  # Check current memory limit
  current_limit <- memory.limit()
  print(paste("Current memory limit:", current_limit, "MB"))
  
  # Set new memory limit to 32GB (32000 MB) - adjust as needed based on your RAM
  memory.limit(size = 32000)
  new_limit <- memory.limit()
  print(paste("New memory limit:", new_limit, "MB"))
} else {
  print("You're on Mac/Linux - R will use system memory automatically.")
  print("If you still get memory errors, try these solutions:")
  print("1. Close other applications to free RAM")
  print("2. Restart R session: Cmd+Shift+F10 (RStudio)")
  print("3. Scale only variable features instead of all genes")
  print("4. Use gc() to force garbage collection")
  
  # Force garbage collection to free memory
  gc()
}
base_dir <- file.path(".")
 if (!dir.exists(base_dir)) {
      stop(paste("ERROR: Directory does not exist:", base_dir,
                 "\nPlease update the base_dir variable to match your system."))
    }

figures_dir <- file.path(base_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir)


getData <- function(fileName,objName){
	ftype=tools::file_ext(fileName)
	if (!file.exists(file.path(base_dir,fileName))){
	  print(paste0(fileName," is not present, Loading from online resources"))
      fileName=paste0("https://home.chpc.utah.edu/~u6004424/Chowanec_Data_2025/",fileName)
	  print(paste0("To speed up future runs download file from ",fileName)) 
      fileName=url(fileName)
	}else{
		print(paste0("Loaded ",fileName))
	}
	if (ftype=="xlsx"){
		assign(objName, read_excel(fileName), envir = .GlobalEnv)
	}
	if (ftype=="RDS" || ftype=="rds"){
		assign(objName, readRDS(fileName), envir = .GlobalEnv);
	}
	if (ftype=="csv"){
		assign(objName, read.csv(fileName), envir = .GlobalEnv)
	}
}

setup <- function(){

	# ============================================================================
	# SETUP AND FILTERING
	# ============================================================================

	# Set data directory
	# IMPORTANT: Update this path to match YOUR computer's directory structure
	# Examples:
	#   Windows: base_dir <- "C:/Users/YourUsername/Desktop/SSM2sc"
	#   Mac/Linux: base_dir <- file.path(path.expand("~"), "Desktop", "SSM2sc")
	# 
	# Or uncomment and modify one of these:
	# base_dir <- "C:/Users/YourUsername/Desktop/SSM2sc"  # Windows - UPDATE USERNAME
	# base_dir <- "/Users/YourUsername/Desktop/SSM2sc"     # Mac - UPDATE USERNAME
	# base_dir <- file.path(path.expand("~"), "Desktop", "SSM2sc")  # Mac/Linux (may work)

	# Define sample names and metadata
	# Batch 1: 1 Placebo, 1 Progesterone
	# Batch 2: 1 Placebo, 2 Progesterone
	sample_names <- c("Batch1_Placebo", "Batch1_Progesterone", 
					  "Batch2_Placebo", "Batch2_Progesterone1", "Batch2_Progesterone2")

	# Define treatment and batch for each sample
	sample_treatments <- c("placebo", "progesterone", 
						   "placebo", "progesterone", "progesterone")
	sample_batches <- c("batch1", "batch1", 
						"batch2", "batch2", "batch2")

	print("Loading 10X data for all samples...")
	print(paste("Total samples to load:", length(sample_names)))

	# Load and filter expression matrices
	expression_matrices <- list()
	for(i in 1:length(sample_names)) {
	  sample_name <- sample_names[i]
	  data_dir <- file.path(base_dir, sample_name)
	  
	  print(paste("Loading sample", i, "of", length(sample_names), ":", sample_name))
	  
	  # Read 10X data
	  expression_matrix <- Read10X(data.dir = data_dir)
	  
	  # Filter out genes expressed in very few cells (quality control)
	  keep_feature <- Matrix::rowSums(expression_matrix)
	  expression_matrix <- expression_matrix[keep_feature > 30,]
	  
	  print(paste("  - Loaded", ncol(expression_matrix), "cells"))
	  print(paste("  - Retained", nrow(expression_matrix), "genes after filtering"))
	  
	  expression_matrices[[sample_name]] <- expression_matrix
	}

	print("All samples loaded successfully!")

	# Merge matrices with unique sample tags
	merged_matrices <- list()
	for(i in 1:length(sample_names)) {
	  sample_name <- sample_names[i]
	  temp_matrix <- expression_matrices[[sample_name]]
	  
	  # Add unique suffix to cell barcodes to distinguish samples
	  colnames(temp_matrix) <- paste0(colnames(temp_matrix), "-", i)
	  merged_matrices[[i]] <- temp_matrix
	}

	print("Finding common genes across all samples...")

	# Find common genes across all samples and merge
	common_genes <- Reduce(intersect, lapply(merged_matrices, rownames))
	print(paste("Found", length(common_genes), "genes common to all samples"))

	merged_matrices_subset <- lapply(merged_matrices, function(x) x[common_genes, ])
	expression_matrix <- do.call(cbind, merged_matrices_subset)

	print(paste("Final merged matrix:", nrow(expression_matrix), "genes x", 
				ncol(expression_matrix), "cells"))

	# Create Seurat object
	print("Creating Seurat object...")
	SSM2sc <- CreateSeuratObject(counts = expression_matrix,
								 project = "SSM2sc",
								 min.cells = 10,
								 min.features = 200)

	print(paste("Seurat object created with", ncol(SSM2sc), "cells"))

	# Add QC metrics
	print("Calculating QC metrics...")
	SSM2sc[["percent.mt"]] <- PercentageFeatureSet(SSM2sc, pattern = "^mt-")
	SSM2sc[["percent.rp"]] <- PercentageFeatureSet(SSM2sc, pattern = "^Rpl|^Rps")

	# Add metadata for treatment and batch
	print("Adding sample metadata...")
	samp_tag <- as.numeric(str_sub(rownames(SSM2sc@meta.data), -1))

	SSM2sc[["treatment"]] <- plyr::mapvalues(samp_tag,
											 from = 1:5,
											 to = sample_treatments)

	SSM2sc[["batch"]] <- plyr::mapvalues(samp_tag,
										 from = 1:5,
										 to = sample_batches)

	SSM2sc[["sampleID"]] <- plyr::mapvalues(samp_tag,
											from = 1:5,
											to = sample_names)

	# Define treatment colors for visualization
	treatment_colors <- c("progesterone" = "#E04B2D", 
						  "placebo" = "#4B2DE0")

	batch_colors <- c("batch1" = "#2DE04B",
					  "batch2" = "#E0A42D")

	# Quality filtering
	print("Applying quality filters...")
	print(paste("Cells before filtering:", ncol(SSM2sc)))

	SSM2sc <- subset(SSM2sc, 
					 subset = nFeature_RNA > 500 & 
					   nCount_RNA > 1000 &
					   percent.mt < 25 & 
					   percent.rp > 1)

	print(paste("Cells after filtering:", ncol(SSM2sc)))

	# Normalize and identify variable features
	print("Normalizing data...")
	SSM2sc <- NormalizeData(SSM2sc)

	print("Finding variable features...")
	SSM2sc <- FindVariableFeatures(SSM2sc, 
								   selection.method = "vst", 
								   nfeatures = 3000)

	# Scale data
	print("Scaling data (this may take a while)...")
	all.genes <- rownames(SSM2sc)
	SSM2sc <- ScaleData(SSM2sc, features = all.genes)

	# Run PCA
	print("Running PCA...")
	SSM2sc <- RunPCA(SSM2sc, 
					 features = VariableFeatures(object = SSM2sc))

	# Batch correction with Harmony
	print("Running Harmony batch correction...")
	SSM2sc <- RunHarmony(SSM2sc, "batch")

	# Clustering and UMAP using Harmony-corrected embeddings
	print("Finding neighbors and clusters...")
	SSM2sc <- FindNeighbors(SSM2sc, 
							reduction = "harmony",
							dims = 1:30)
	SSM2sc <- FindClusters(SSM2sc, resolution = 2.0)

	print("Running UMAP...")
	SSM2sc <- RunUMAP(SSM2sc, 
					  reduction = "harmony",
					  dims = 1:30, 
					  n.neighbors = 10, 
					  min.dist = 0.3)

	Idents(SSM2sc) <- SSM2sc$seurat_clusters

	# Remove outlier clusters
	print("Cluster cell counts:")
	cluster_counts <- table(SSM2sc$seurat_clusters)
	print(cluster_counts)

	# Remove specific outlier clusters
	exclude <- c("34", "36", "37", "38", "39", "40")
	print(paste("Removing outlier clusters:", paste(exclude, collapse = ", ")))

	SSM2sc[["keep"]] <- plyr::mapvalues(!(SSM2sc$seurat_clusters %in% exclude),
										c(FALSE, TRUE), c("discard", "keep"))
	SSM2sc <- subset(SSM2sc, subset = keep == "keep")

	print(paste("After removing outlier clusters:", ncol(SSM2sc), "cells remain"))

	# Create output directory
	# Note: figures_dir is built from base_dir, so it will automatically use the correct path
	figures_dir <- file.path(base_dir, "Figures")
	if (!dir.exists(figures_dir)) dir.create(figures_dir)

	# Verify the directory was created successfully
	if (!dir.exists(figures_dir)) {
	  stop(paste("ERROR: Could not create directory:", figures_dir,
				 "\nCheck that you have write permissions for this location."))
	}

	# Save filtered object
	print("Saving filtered Seurat object...")
	saveRDS(SSM2sc, file.path(base_dir, "SSM2sc_filtered.RDS"))

	print("=== FILTERED OBJECT SAVED ===")
	print(paste("Filtered object saved as:", file.path(base_dir, "SSM2sc_filtered.RDS")))
	print(paste("Final cell count:", ncol(SSM2sc)))
	print(paste("Final gene count:", nrow(SSM2sc)))
	print(paste("Number of clusters:", length(unique(SSM2sc$seurat_clusters))))
	print("Setup and filtering complete!")

	# Find marker genes for all clusters
	print("Finding marker genes for all clusters (this may take a while)...")
	all_markers <- FindAllMarkers(SSM2sc,
								  only.pos = TRUE, # Only positive markers
								  min.pct = 0.25, # Gene must be expressed in at least 25% of cells
								  logfc.threshold = 0.25, # Minimum log fold change
								  test.use = "wilcox") # Statistical test

	# View the results
	print("Top 20 markers across all clusters:")
	head(all_markers, 20)

	# Check the structure
	str(all_markers)

	# Generate top 10 markers per cluster
	top10_markers <- all_markers %>%
	  group_by(cluster) %>%
	  top_n(n = 10, wt = avg_log2FC) %>%
	  arrange(cluster, desc(avg_log2FC))

	# Generate top 5 markers per cluster
	top5_markers <- all_markers %>%
	  group_by(cluster) %>%
	  top_n(n = 5, wt = avg_log2FC) %>%
	  arrange(cluster, desc(avg_log2FC))

	# Save all markers
	print("Saving marker gene results...")
	write.csv(all_markers, file.path(figures_dir, "all_cluster_markers.csv"), row.names = FALSE)

	# Save top 10 markers
	write.csv(top10_markers, file.path(figures_dir, "top10_cluster_markers.csv"), row.names = FALSE)

	# Save top 5 markers
	write.csv(top5_markers, file.path(figures_dir, "top5_cluster_markers.csv"), row.names = FALSE)

	print("=== MARKER GENES SAVED ===")
	print(paste("All markers saved as:", file.path(figures_dir, "all_cluster_markers.csv")))
	print(paste("Top 10 markers saved as:", file.path(figures_dir, "top10_cluster_markers.csv")))
	print(paste("Top 5 markers saved as:", file.path(figures_dir, "top5_cluster_markers.csv")))

	# ============================================================================
	# END OF SETUP AND FILTERING SECTION
	# ============================================================================
}


figure2 <- function(){
	print("Generating figure2")
		getData("SSM2sc_filtered.RDS","SSM2sc")

	# ============================================================================
	# FIGURE 2: CELL TYPE IDENTIFICATION AND VISUALIZATION
	# ============================================================================

	# Create Figure 2 subdirectory
	# Note: figure2_dir is built from base_dir, so it will automatically use the correct path
	figure2_dir <- file.path(base_dir, "Figures", "Figure_2")
	if (!dir.exists(figure2_dir)) dir.create(figure2_dir, recursive = TRUE)

	# Verify the directory was created successfully
	if (!dir.exists(figure2_dir)) {
	  stop(paste("ERROR: Could not create directory:", figure2_dir,
				 "\nCheck that you have write permissions for this location."))
	}

	print("=== FIGURE 2: RENAMING CLUSTERS WITH CELL TYPE LABELS ===")

	# Define cluster to cell type mapping (based on marker gene analysis)
	cluster_to_celltype <- c(
	  "0" = "CD4+ 2",
	  "1" = "Tumor 1",
	  "2" = "Macrophage 1",
	  "3" = "Tumor 3",
	  "4" = "Tumor 2",
	  "5" = "Tumor 4",
	  "6" = "Tumor 5",
	  "7" = "CD8+ 1",
	  "8" = "CD4+ 1",
	  "9" = "Tumor 6",
	  "10" = "Tumor 7",
	  "11" = "CD4+ Treg",
	  "12" = "Tumor 8",
	  "13" = "Tumor 9",
	  "14" = "Macrophage 2",
	  "15" = "Epithelial 1",
	  "16" = "Epithelial 2",
	  "17" = "B cell",
	  "18" = "Tumor 11",
	  "19" = "Tumor 10",
	  "20" = "Fibroblast 1",
	  "21" = "Macrophage 3",
	  "22" = "Macrophage 4",
	  "23" = "Neutrophil",
	  "24" = "NK",
	  "25" = "Macrophage 5",
	  "26" = "Basal",
	  "27" = "cDC",
	  "28" = "CD8+ 2",
	  "29" = "Epithelial 3",
	  "30" = "pDC",
	  "31" = "Epithelial 4",
	  "32" = "Fibroblast 2",
	  "33" = "Macrophage 6",
	  "35" = "Endothelial"
	)

	# Define hex colors for each cell type (based on Excel file)
	celltype_colors <- c(
	  "CD4+ 2" = "#FA0087",
	  "Tumor 1" = "#1C7F93",
	  "Macrophage 1" = "#16FF32",
	  "Tumor 3" = "#683B79",
	  "Tumor 2" = "#FC1CBF",
	  "Tumor 4" = "#F6222E",
	  "Tumor 5" = "#D85FF7",
	  "CD8+ 1" = "#2ED9FF",
	  "CD4+ 1" = "#F8A19F",
	  "Tumor 6" = "#C075A6",
	  "Tumor 7" = "#B10DA1",
	  "CD4+ Treg" = "#5A5156",
	  "Tumor 8" = "#B00068",
	  "Tumor 9" = "#782AB6",
	  "Macrophage 2" = "#3283FE",
	  "Epithelial 1" = "#1C8356",
	  "Epithelial 2" = "#FE00FA",
	  "B cell" = "#1CFFCE",
	  "Tumor 11" = "#822E1C",
	  "Tumor 10" = "#66B0FF",
	  "Fibroblast 1" = "#85660D",
	  "Macrophage 3" = "#7ED7D1",
	  "Macrophage 4" = "#325A9B",
	  "Neutrophil" = "#90AD1C",
	  "NK" = "#FEAF16",
	  "Macrophage 5" = "#DEA0FD",
	  "Basal" = "#1CBE4F",
	  "cDC" = "#B5EFB5",
	  "CD8+ 2" = "#E4E1E3",
	  "Epithelial 3" = "#F7E1A0",
	  "pDC" = "#BDCDFF",
	  "Epithelial 4" = "#AAF400",
	  "Fibroblast 2" = "#FBE426",
	  "Macrophage 6" = "#C4451C",
	  "Endothelial" = "#AA0DFE"
	)

	# Add cell type labels to metadata
	SSM2sc$cell_type_secondary <- plyr::mapvalues(
	  x = as.character(SSM2sc$seurat_clusters),
	  from = names(cluster_to_celltype),
	  to = as.character(cluster_to_celltype)
	)

	# Verify the mapping
	print("Cell type distribution:")
	print(table(SSM2sc$cell_type_secondary))

	# Create order vector for colors (matching the order in the Seurat object)
	unique_celltypes <- unique(SSM2sc$cell_type_secondary)
	order_vec <- match(unique_celltypes, names(celltype_colors))

	# Create UMAP by cell type
	print("Creating UMAP visualization...")
	p00 <- DimPlot(SSM2sc, 
				   group.by = "cell_type_secondary",
				   cols = celltype_colors[order_vec],
				   combine = TRUE, 
				   label = FALSE, 
				   pt.size = 0.5) & 
	  NoAxes() &
	  labs(title = "UMAP by Cell Type") &
	  theme(plot.title = element_text(hjust = 0, size = 30, face = "plain"), 
			legend.text = element_text(size = 28),
			axis.title.y.right = element_text(size = 20)) &
	  guides(color = guide_legend(override.aes = list(size = 8), ncol = 2))

	# Save the UMAP by cell type
	ggsave(file.path(figure2_dir, "UMAP_by_celltype.png"), 
		   p00, 
		   width = 16, 
		   height = 12, 
		   dpi = 300)

	print("UMAP by cell type saved!")

	treatment_colors <- c("progesterone" = "#E04B2D","placebo" = "#4B2DE0")
	batch_colors <- c("batch1" = "#2DE04B","batch2" = "#E0A42D")


	# Create UMAP by treatment
	print("Creating UMAP by treatment...")
	p01 <- DimPlot(SSM2sc, 
				   group.by = "treatment", 
				   cols = treatment_colors,
				   combine = TRUE, 
				   label = FALSE, 
				   pt.size = 0.4) &
	  NoAxes() &
	  ylim(-13, 13) &
	  labs(title = "Treatment") & 
	  theme(plot.title = element_text(hjust = 0, size = 30, face = "plain"), 
			legend.text = element_text(size = 28),
			axis.title.y.right = element_text(size = 20)) &
	  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1))

	# Save the treatment UMAP plot
	ggsave(file.path(figure2_dir, "UMAP_by_treatment.png"), 
		   p01, 
		   width = 16, 
		   height = 12, 
		   dpi = 300)

	print("UMAP by treatment saved!")

	# Create cell type category mapping for primary cell types
	print("Creating primary cell type category mapping...")

	# Define categories for each cell type 
	category_mapping <- c(
	  "Tumor 1" = "Tumor",
	  "Tumor 2" = "Tumor",
	  "Tumor 3" = "Tumor",
	  "Tumor 4" = "Tumor",
	  "Tumor 5" = "Tumor",
	  "Tumor 6" = "Tumor",
	  "Tumor 7" = "Tumor",
	  "Tumor 8" = "Tumor",
	  "Tumor 9" = "Tumor",
	  "Tumor 10" = "Tumor",
	  "Tumor 11" = "Tumor",
	  "CD4+ 1" = "Immune",
	  "CD4+ 2" = "Immune",
	  "CD4+ Treg" = "Immune",
	  "CD8+ 1" = "Immune",
	  "CD8+ 2" = "Immune",
	  "B cell" = "Immune",
	  "NK" = "Immune",
	  "Neutrophil" = "Immune",
	  "Macrophage 1" = "Immune",
	  "Macrophage 2" = "Immune",
	  "Macrophage 3" = "Immune",
	  "Macrophage 4" = "Immune",
	  "Macrophage 5" = "Immune",
	  "Macrophage 6" = "Immune",
	  "cDC" = "Immune",
	  "pDC" = "Immune",
	  "Fibroblast 1" = "Stromal",
	  "Fibroblast 2" = "Stromal",
	  "Endothelial" = "Stromal",
	  "Epithelial 1" = "Non-tumor epithelial",
	  "Epithelial 2" = "Non-tumor epithelial",
	  "Epithelial 3" = "Non-tumor epithelial",
	  "Epithelial 4" = "Non-tumor epithelial",
	  "Basal" = "Non-tumor epithelial"
	)

	# Add cell type category to metadata
	SSM2sc$cell_type_primary_final <- plyr::mapvalues(
	  x = as.character(SSM2sc$cell_type_secondary),
	  from = names(category_mapping),
	  to = as.character(category_mapping)
	)

	# Convert to factor with specified level order
	SSM2sc$cell_type_primary_final <- factor(
	  SSM2sc$cell_type_primary_final,
	  levels = c("Immune", "Tumor", "Stromal", "Non-tumor epithelial")
	)

	# Verify the category mapping
	print("Primary cell type distribution:")
	print(table(SSM2sc$cell_type_primary_final))

	# Define colors for primary categories (matching your specified colors)
	primary_colors_final <- c(
	  "Immune" = "#a6cee3",              # Light blue
	  "Tumor" = "#e31a1c",               # Red
	  "Stromal" = "#b2df8a",             # Light green
	  "Non-tumor epithelial" = "#fb9a99" # Light pink
	)

	# Create side-by-side UMAP split by treatment (p02)
	print("Creating side-by-side UMAP by primary cell type and treatment...")

	p_placebo <- DimPlot(SSM2sc, 
						 group.by = "cell_type_primary_final", 
						 cols = primary_colors_final,
						 cells = which(SSM2sc@meta.data$treatment == "placebo"),
						 combine = TRUE, 
						 label = FALSE, 
						 pt.size = 0.6) &
	  NoAxes() &
	  labs(title = "Placebo") &
	  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
			legend.text = element_text(size = 16),
			legend.title = element_text(size = 18)) &
	  guides(color = guide_legend(override.aes = list(size = 5), title = "Primary Cell Type"))

	p_progesterone <- DimPlot(SSM2sc, 
							  group.by = "cell_type_primary_final", 
							  cols = primary_colors_final,
							  cells = which(SSM2sc@meta.data$treatment == "progesterone"),
							  combine = TRUE, 
							  label = FALSE, 
							  pt.size = 0.6) &
	  NoAxes() &
	  labs(title = "Progesterone") &
	  theme(plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
			legend.text = element_text(size = 16),
			legend.title = element_text(size = 18)) &
	  guides(color = guide_legend(override.aes = list(size = 5), title = "Primary Cell Type"))

	# Combine side-by-side plots
	p02 <- p_placebo + p_progesterone + 
	  plot_layout(guides = "collect") +
	  plot_annotation(title = "Primary Cell Type",
					  theme = theme(plot.title = element_text(hjust = 0.5, size = 28, face = "bold")))

	# Save the combined plot
	ggsave(file.path(figure2_dir, "UMAP_primary_celltype_by_treatment.png"), 
		   p02, 
		   width = 16, 
		   height = 8, 
		   dpi = 300)

	print("Side-by-side UMAP (primary cell type by treatment) saved!")

	# Create stacked bar plot of primary cell type proportions by treatment
	print("Creating stacked bar plot of cell type proportions...")

	# Calculate proportions by treatment
	cluster_treatment_counts <- SSM2sc@meta.data %>%
	  group_by(treatment, cell_type_primary_final) %>%
	  summarise(count = n(), .groups = 'drop') %>%
	  group_by(treatment) %>%
	  mutate(proportion = count / sum(count))

	# Create stacked bar plot
	p03 <- ggplot(cluster_treatment_counts, aes(x = treatment, y = proportion, fill = cell_type_primary_final)) +
	  geom_bar(stat = "identity", position = "stack", width = 0.6) +
	  scale_fill_manual(values = primary_colors_final) +
	  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00), 
						 labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
	  labs(x = "", y = "Proportion", fill = "Primary Cell Type") +
	  theme_minimal() +
	  theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "white", color = NA),
		plot.background = element_rect(fill = "white", color = NA),
		axis.text.x = element_text(size = 14, color = "black"),
		axis.text.y = element_text(size = 12, color = "black"),
		axis.ticks = element_line(color = "black"),
		legend.position = "right",
		legend.text = element_text(size = 11),
		legend.margin = margin(l = 20),
		plot.margin = margin(t = 10, r = 20, b = 10, l = 50)
	  ) +
	  coord_flip()

	# Display and save plot
	print(p03)
	ggsave(file.path(figure2_dir, "Stacked_barplot_primary_celltype_by_treatment.png"), 
		   p03, 
		   width = 10, 
		   height = 4, 
		   dpi = 300)

	# Print proportion table
	print("Cell type proportions by treatment:")
	print(cluster_treatment_counts)

	# Save proportion table as CSV
	write.csv(cluster_treatment_counts, 
			  file.path(figure2_dir, "primary_celltype_proportions_by_treatment.csv"), 
			  row.names = FALSE)

	print("Stacked bar plot saved!")

	# ============================================================================
	# STATISTICAL ANALYSIS: Fisher's Exact Test with FDR Correction
	# ============================================================================

	print("Performing statistical analysis...")

	# Get treatment totals
	treatment_totals <- table(SSM2sc@meta.data$treatment)
	treatments <- names(treatment_totals)

	cat("\nTotal cells per treatment:\n")
	print(treatment_totals)

	# Get all primary cell types
	primary_cell_types <- levels(SSM2sc$cell_type_primary_final)

	# Initialize results dataframe
	results <- data.frame(
	  cell_type = character(),
	  placebo_count = numeric(),
	  progesterone_count = numeric(),
	  placebo_prop = numeric(),
	  progesterone_prop = numeric(),
	  fold_change = numeric(),
	  p_value = numeric(),
	  stringsAsFactors = FALSE
	)

	# Perform Fisher's exact test for each cell type
	for(cell_type in primary_cell_types) {
	  
	  # Get counts for this cell type
	  placebo_this_type <- sum(SSM2sc@meta.data$treatment == "placebo" & 
								 SSM2sc@meta.data$cell_type_primary_final == cell_type)
	  progesterone_this_type <- sum(SSM2sc@meta.data$treatment == "progesterone" & 
									  SSM2sc@meta.data$cell_type_primary_final == cell_type)
	  
	  # Get counts for all other cell types combined
	  placebo_other_types <- treatment_totals["placebo"] - placebo_this_type
	  progesterone_other_types <- treatment_totals["progesterone"] - progesterone_this_type
	  
	  # Create contingency table
	  contingency_table <- matrix(c(placebo_this_type, placebo_other_types,
									progesterone_this_type, progesterone_other_types),
								  nrow = 2, byrow = TRUE)
	  
	  # Perform Fisher's exact test
	  fisher_test <- fisher.test(contingency_table)
	  
	  # Calculate proportions
	  placebo_prop <- placebo_this_type / treatment_totals["placebo"]
	  progesterone_prop <- progesterone_this_type / treatment_totals["progesterone"]
	  
	  # Calculate fold change (progesterone/placebo)
	  fold_change <- ifelse(placebo_prop > 0, progesterone_prop / placebo_prop, NA)
	  
	  # Add to results
	  results <- rbind(results, data.frame(
		cell_type = cell_type,
		placebo_count = placebo_this_type,
		progesterone_count = progesterone_this_type,
		placebo_prop = placebo_prop,
		progesterone_prop = progesterone_prop,
		fold_change = fold_change,
		p_value = fisher_test$p.value
	  ))
	}

	# Apply FDR correction
	results$p_adj <- p.adjust(results$p_value, method = "fdr")

	# Add significance indicators
	results$significance <- ifelse(results$p_adj < 0.001, "***",
								   ifelse(results$p_adj < 0.01, "**",
										  ifelse(results$p_adj < 0.05, "*", "ns")))

	# Sort by adjusted p-value
	results <- results[order(results$p_adj), ]

	# Display results
	cat("\n=== PRIMARY CELL TYPE STATISTICAL RESULTS ===\n")
	cat("Fisher's Exact Test with FDR correction\n")
	cat("Comparing progesterone vs placebo for each primary cell type\n\n")

	print(results)

	# Save results
	write.csv(results, 
			  file.path(figure2_dir, "primary_celltype_statistical_comparison.csv"), 
			  row.names = FALSE)

	# Identify significant differences
	significant_results <- results[results$p_adj < 0.05, ]

	if(nrow(significant_results) > 0) {
	  cat("\n=== SIGNIFICANT DIFFERENCES (FDR p < 0.05) ===\n")
	  for(i in 1:nrow(significant_results)) {
		row <- significant_results[i, ]
		direction <- ifelse(row$fold_change > 1, "INCREASED", "DECREASED")
		cat(sprintf("%s: %s in progesterone vs placebo (FC=%.3f, FDR p=%.2e)\n", 
					row$cell_type, direction, row$fold_change, row$p_adj))
	  }
	} else {
	  cat("\n=== NO SIGNIFICANT DIFFERENCES (FDR p < 0.05) ===\n")
	}

	print("Statistical analysis complete!")

	# Save updated Seurat object with cell type labels and categories
	saveRDS(SSM2sc, file.path(base_dir, "SSM2sc_with_celltypes.RDS"))

	print("=== FIGURE 2 COMPLETE ===")
	print(paste("Updated object saved as:", file.path(base_dir, "SSM2sc_with_celltypes.RDS")))
	print(paste("UMAP by cell type saved to:", file.path(figure2_dir, "UMAP_by_celltype.png")))
	print(paste("UMAP by treatment saved to:", file.path(figure2_dir, "UMAP_by_treatment.png")))
	print(paste("UMAP primary cell type by treatment saved to:", file.path(figure2_dir, "UMAP_primary_celltype_by_treatment.png")))
	print(paste("Stacked bar plot saved to:", file.path(figure2_dir, "Stacked_barplot_primary_celltype_by_treatment.png")))
	print(paste("Statistical results saved to:", file.path(figure2_dir, "primary_celltype_statistical_comparison.csv")))

	# ============================================================================
	# SUPPLEMENTAL FIGURES EPCAM PTPRC EXPRESSION
	# ============================================================================

	# Create Supplemental Figures subdirectory
	supplemental_dir <- file.path(base_dir, "Figures", "Supplemental_figures")
	if (!dir.exists(supplemental_dir)) dir.create(supplemental_dir, recursive = TRUE)

	# Verify the directory was created successfully
	if (!dir.exists(supplemental_dir)) {
	  stop(paste("ERROR: Could not create directory:", supplemental_dir,
				 "\nCheck that you have write permissions for this location."))
	}

	print("=== SUPPLEMENTAL FIGURES ===")

	# Create feature plots for Epcam and Ptprc expression
	print("Creating UMAP for Epcam and Ptprc expression...")

	# Check if genes exist in the dataset
	genes_to_plot <- c("Epcam", "Ptprc")
	genes_present <- genes_to_plot[genes_to_plot %in% rownames(SSM2sc)]

	if(length(genes_present) == 0) {
	  print("Warning: Neither Epcam nor Ptprc found in dataset")
	} else {
	  print(paste("Found genes:", paste(genes_present, collapse = ", ")))
	  
	  # Create feature plot for both genes with blue-to-red gradient
	  # order = TRUE plots high-expressing cells on top
	  p_epcam_ptprc <- FeaturePlot(SSM2sc, 
								   features = genes_present,
								   pt.size = 0.3,
								   ncol = 2,
								   order = TRUE,
								   cols = c("blue", "lightblue", "white", "orange", "red")) &
		NoAxes() &
		theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
	  
	  # Display the plot
	  print(p_epcam_ptprc)
	  
	  # Save the plot
	  ggsave(file.path(supplemental_dir, "UMAP_Epcam_Ptprc_expression.png"), 
			 p_epcam_ptprc, 
			 width = 14, 
			 height = 6, 
			 dpi = 300)
	  
	  print("UMAP for Epcam and Ptprc saved to Supplemental_figures!")
	}


	print("=== SUPPLEMENTAL FIGURES COMPLETE ===")
	print(paste("Epcam and Ptprc UMAP saved to:", file.path(supplemental_dir, "UMAP_Epcam_Ptprc_expression.png")))

}



figure3 <- function(){
	print("Generating figure3")
		getData("SSM2sc_with_celltypes.RDS","SSM2sc")


	# ============================================================================
	# FIGURE 3: VOLCANO PLOT - TUMOR CELLS DIFFERENTIAL EXPRESSION
	# ============================================================================

	# Create Figure 3 subdirectory
	figure3_dir <- file.path(base_dir, "Figures", "Figure_3")
	if (!dir.exists(figure3_dir)) dir.create(figure3_dir, recursive = TRUE)

	# Verify the directory was created successfully
	if (!dir.exists(figure3_dir)) {
	  stop(paste("ERROR: Could not create directory:", figure3_dir,
				 "\nCheck that you have write permissions for this location."))
	}

	print("=== FIGURE 3: TUMOR DIFFERENTIAL EXPRESSION ANALYSIS ===")

	# Subset to tumor cells only (all Tumor 1 through Tumor 11 clusters)
	print("Subsetting to tumor cells...")
	SSM2sc_tumor <- subset(SSM2sc, subset = cell_type_primary_final == "Tumor")

	# Check how many cells we have per treatment in tumor cells
	print("Tumor cells per treatment:")
	print(table(SSM2sc_tumor$treatment))

	# Filter mitochondrial and ribosomal genes
	print("Filtering mitochondrial and ribosomal genes...")
	rm.ind <- grep("^mt-|^Rp[ls]", rownames(SSM2sc_tumor))
	if(length(rm.ind) > 0) {
	  keep <- rownames(SSM2sc_tumor)[-rm.ind]
	  SSM2sc_tumor <- SSM2sc_tumor[keep,]
	}

	# Filter low-expression genes (Seurat v5 compatible)
	print("Filtering low-expression genes...")
	counts_data <- GetAssayData(SSM2sc_tumor, assay = "RNA", layer = "counts")
	keep.ind <- rowSums(counts_data)
	keep.ind <- names(keep.ind[keep.ind >= 20])
	SSM2sc_tumor <- SSM2sc_tumor[keep.ind,]

	# Save background genes for enrichment analysis
	write.csv(keep.ind, file.path(figure3_dir, "tumor_background.csv"), row.names = FALSE)
	print(paste("Background genes saved:", length(keep.ind), "genes"))

	# Set identities to treatment
	Idents(SSM2sc_tumor) <- "treatment"

	# Perform differential expression analysis
	print("Performing differential expression analysis (this may take several minutes)...")
	print("Comparing: Progesterone vs Placebo")

	cluster0.markers0 <- FindMarkers(SSM2sc_tumor, 
									 ident.1 = "progesterone",
									 ident.2 = "placebo",
									 group.by = "treatment",
									 only.pos = FALSE,
									 test.use = "MAST",
									 latent.vars = c("nCount_RNA"))

	# Filter for significant results
	cluster0.markers <- cluster0.markers0[cluster0.markers0$p_val_adj < 0.05 & 
											abs(cluster0.markers0$avg_log2FC) > 0.25,]

	# Add gene names as a column
	cluster0.markers <- tibble::rownames_to_column(cluster0.markers, var = "Gene")

	# Save results
	saveTag <- file.path(figure3_dir, "tumor_progesterone_vs_placebo.csv")
	write.csv(cluster0.markers, saveTag, row.names = FALSE)

	# Print summary
	cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
	cat("Total tumor cells analyzed:", ncol(SSM2sc_tumor), "\n")
	cat("Found", nrow(cluster0.markers), "significantly differentially expressed genes\n")
	cat("Upregulated in progesterone:", sum(cluster0.markers$avg_log2FC > 0), "\n")
	cat("Downregulated in progesterone:", sum(cluster0.markers$avg_log2FC < 0), "\n")
	cat("Results saved to:", saveTag, "\n\n")

	# ============================================================================
	# CREATE VOLCANO PLOT
	# ============================================================================

	print("Creating volcano plot...")

	# Load EnhancedVolcano library
	if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
	  if (!requireNamespace("BiocManager", quietly = TRUE)) {
		install.packages("BiocManager")
	  }
	  BiocManager::install("EnhancedVolcano")
	}
	loadLibrary(EnhancedVolcano)

	# Create color scheme based on fold change direction
	keyvals <- ifelse(
	  cluster0.markers$avg_log2FC < 0, 'royalblue',
	  ifelse(cluster0.markers$avg_log2FC > 0, 'red3', 'black'))
	keyvals[is.na(keyvals)] <- 'black'
	names(keyvals)[keyvals == 'red3'] <- 'high'
	names(keyvals)[keyvals == 'black'] <- 'mid'
	names(keyvals)[keyvals == 'royalblue'] <- 'low'

	# Select downregulated genes (lower in progesterone vs placebo)
	sel_labs <- (cluster0.markers$Gene)[which(names(keyvals) %in% c('low'))]

	# Filter for specific immune-related gene families
	sel_labs <- grep("Ifi|Isg|Irf|Cxc|Zbp|Gbp|Stat|Irgm|Irg|Ly6|Ddx|Pd",
					 sel_labs, value = TRUE)

	# Remove Stat1 from the labels
	sel_labs <- sel_labs[sel_labs != "Stat1"]

	# Print how many genes will be labeled
	cat("Number of immune-related genes to label:", length(sel_labs), "\n")
	cat("Genes to be labeled:\n")
	print(sel_labs)

	# Create the volcano plot
	p1 <- EnhancedVolcano(cluster0.markers,
						  lab = cluster0.markers$Gene,
						  selectLab = sel_labs,
						  x = "avg_log2FC",
						  y = "p_val_adj",
						  xlim = c(-2, 2),
						  ylim = c(0, 260),
						  caption = NULL,
						  axisLabSize = 30,
						  xlab = bquote(~Log[2]~ "fold change"),
						  ylab = bquote(~-Log[10]~"adjusted"~italic(p)),
						  title = NULL,
						  subtitle = "",
						  pCutoff = 0.05,
						  pointSize = 2.0,
						  labSize = 9,
						  labCol = 'black',
						  labFace = 'plain',
						  drawConnectors = TRUE,
						  arrowheads = FALSE,
						  boxedLabels = FALSE,
						  colCustom = keyvals,
						  FCcutoff = 0.25,
						  legendPosition = "none")

	# Display the plot
	print(p1)

	### VOLCANO PLOT - UP GENES

# Create color scheme based on fold change direction
keyvals <- ifelse(
  cluster0.markers$avg_log2FC < 0, 'royalblue',
  ifelse(cluster0.markers$avg_log2FC > 0, 'red3', 'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red3'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

# New gene list to label
UP_labs <- c("Cdc14a", "Smc3", "Cdc14b", "Pten", "Cep290", "Mapre1",
             "Lzts2", "Ift46", "Mns1", "Ube2b", "Aaas", "Mark2", "Map4", "Rae1",
             "C2cd3", "Stmn1", "Bbs4", "Chmp4b", "Cep70", "Dnm2", "Hspa1b",
             "Aurka", "Chmp3", "Cenph", "Ddb1", "Clasp2", "Kif2a", "Wnt4")


# Remove the 4 genes from selectLab since we'll add them manually
UP_labs_main <- UP_labs[!UP_labs %in% c("Wnt4", "Cdc14a")]

v1 <- EnhancedVolcano(cluster0.markers,
                      lab = cluster0.markers$Gene,
                      selectLab = UP_labs_main,
                      x = "avg_log2FC",
                      y = "p_val_adj",
                      xlim = c(-2, 2),
                      ylim = c(0, 260),
                      caption = NULL,
                      axisLabSize = 30,
                      xlab = bquote(~Log[2]~ "fold change"),
                      ylab = bquote(~-Log[10]~"adjusted"~italic(p)),
                      title = NULL,
                      subtitle = "",
                      pCutoff = 0.05,
                      pointSize = 2.0,
                      labSize = 9,
                      labCol = 'black',
                      labFace = 'plain',
                      drawConnectors = TRUE,
                      arrowheads = FALSE,
                      boxedLabels = FALSE,
                      colCustom = keyvals,
                      FCcutoff = 0.25,
                      legendPosition = "none",
                      max.overlaps = Inf) +
  annotate("text", x = 1.05, y = 255, label = "Wnt4", size = 9) +
  annotate("text", x = 1.64, y = 250, label = "Cdc14a", size = 9)

   print(v1)

	# Save the plot
	ggsave(file.path(figure3_dir, "tumor_progesterone_vs_placebo_volcano.pdf"),
		   p1, width = 9, height = 12)

	ggsave(file.path(figure3_dir, "tumor_progesterone_vs_placebo_volcano.png"),
		   p1, width = 9, height = 12, dpi = 300)

	ggsave(file.path(figure3_dir, "Volcano_up_regulegulated_genes.pdf"),
		   v1,width=9,height=12)

	ggsave(file.path(figure3_dir, "Volcano_up_regulegulated_genes.png"),
		   v1, width = 9, height = 12, dpi = 300)

	print("=== FIGURE 3 COMPLETE ===")
	print(paste("Differential expression results saved to:", saveTag))
	print(paste("Volcano plot saved to:", file.path(figure3_dir, "tumor_progesterone_vs_placebo_volcano.png")))

	# ============================================================================
	# SUPPLEMENTAL VOLCANO PLOT ANALYSIS - SETUP INSTRUCTIONS
	# ============================================================================
	#
	# Before running this code, you need to set up your file structure and 
	# download the required data files.
	#
	# STEP 1: CREATE DIRECTORY STRUCTURE
	# -----------------------------------
	# Create the following folder structure on your Desktop:
	#
	# Desktop/
	# └── SSM2sc/
	#     ├── Figures/
	#     │   ├── Figure_3/
	#     │   │   └── tumor_progesterone_vs_placebo.csv
	#     │   └── Supplemental_figures/  (will be created automatically)
	#     └── overlapping_genes_estrogen_progesterone.csv
	#
	# STEP 2: DOWNLOAD AND SAVE REQUIRED FILES
	# -----------------------------------------
	# File 1: Differential Expression Data
	#   - File: tumor_progesterone_vs_placebo.csv
	#   - Download from: [Your data source/repository link]
	#   - Save to: Desktop/SSM2sc/Figures/Figure_3/tumor_progesterone_vs_placebo.csv
	#
	# File 2: Overlapping Estrogen Response Genes
	#   - File: overlapping_genes_estrogen_progesterone.csv
	#   - Download from: [Your data source/repository link]
	#   - Save to: Desktop/SSM2sc/overlapping_genes_estrogen_progesterone.csv
	#   - IMPORTANT: This file should be in the ROOT SSM2sc folder, NOT in a subfolder
	#
	# STEP 3: VERIFY FILE PLACEMENT
	# ------------------------------
	# Windows users verify:
	#   C:/Users/[YourUsername]/Desktop/SSM2sc/overlapping_genes_estrogen_progesterone.csv
	#   C:/Users/[YourUsername]/Desktop/SSM2sc/Figures/Figure_3/tumor_progesterone_vs_placebo.csv
	#
	# Mac/Linux users verify:
	#   ~/Desktop/SSM2sc/overlapping_genes_estrogen_progesterone.csv
	#   ~/Desktop/SSM2sc/Figures/Figure_3/tumor_progesterone_vs_placebo.csv
	#
	# STEP 4: RUN THE CODE
	# --------------------
	# The code will:
	#   1. Automatically install required packages (EnhancedVolcano) if needed
	#   2. Read both CSV files from the correct locations
	#   3. Generate the volcano plot
	#   4. Automatically create the Supplemental_figures folder if it doesn't exist
	#   5. Save output files to: Desktop/SSM2sc/Figures/Supplemental_figures/
	#
	# EXPECTED OUTPUT
	# ---------------
	# Two files will be created:
	#   - volcano_progesterone_estrogen_overlap.pdf
	#   - volcano_progesterone_estrogen_overlap.png
	# Both saved in: Desktop/SSM2sc/Figures/Supplemental_figures/
	#
	# ============================================================================
	# BEGIN CODE FOR SUPPLEMENTAL VOLCANO 
	# ============================================================================

	print("Creating volcano plot with estrogen overlap genes...")

	# Load EnhancedVolcano library
	if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
	  if (!requireNamespace("BiocManager", quietly = TRUE)) {
		install.packages("BiocManager")
	  }
	  BiocManager::install("EnhancedVolcano")
	}
	loadLibrary(EnhancedVolcano)


	# Read the differential expression data
	cluster0.markers <- read.csv(file.path(base_dir, "Figures", "Figure_3", "tumor_progesterone_vs_placebo.csv"))

	# Read the overlapping estrogen response genes
	getData("overlapping_genes_estrogen_progesterone.csv","overlap_genes")
	genes_to_label <- overlap_genes$Gene

	cat("Total overlapping estrogen response genes:", length(genes_to_label), "\n")
	cat("Genes to be labeled:\n")
	print(genes_to_label)

	# Create color scheme based on fold change direction
	keyvals <- ifelse(
	  cluster0.markers$avg_log2FC < 0, 'royalblue',
	  ifelse(cluster0.markers$avg_log2FC > 0, 'red3', 'black'))
	keyvals[is.na(keyvals)] <- 'black'
	names(keyvals)[keyvals == 'red3'] <- 'high'
	names(keyvals)[keyvals == 'black'] <- 'mid'
	names(keyvals)[keyvals == 'royalblue'] <- 'low'

	# Use the overlapping estrogen response genes for labeling
	sel_labs <- genes_to_label

	# Count how many overlap genes are up vs down
	overlap_up <- sum(cluster0.markers$Gene %in% genes_to_label & cluster0.markers$avg_log2FC > 0)
	overlap_down <- sum(cluster0.markers$Gene %in% genes_to_label & cluster0.markers$avg_log2FC < 0)

	cat("\nBreakdown of overlap genes:\n")
	cat("  Upregulated in progesterone:", overlap_up, "\n")
	cat("  Downregulated in progesterone:", overlap_down, "\n")

	# Create the volcano plot
	p1 <- EnhancedVolcano(cluster0.markers,
						  lab = cluster0.markers$Gene,
						  selectLab = sel_labs,
						  x = "avg_log2FC",
						  y = "p_val_adj",
						  xlim = c(-2, 2),
						  ylim = c(0, 260),
						  caption = NULL,
						  axisLabSize = 30,
						  xlab = bquote(~Log[2]~ "fold change"),
						  ylab = bquote(~-Log[10]~"adjusted"~italic(p)),
						  title = "Progesterone vs Placebo - Estrogen Response Overlap",
						  subtitle = paste(length(genes_to_label), "overlapping genes labeled"),
						  pCutoff = 0.05,
						  pointSize = 2.0,
						  labSize = 9,
						  labCol = 'black',
						  labFace = 'plain',
						  drawConnectors = TRUE,
						  arrowheads = FALSE,
						  boxedLabels = FALSE,
						  colCustom = keyvals,
						  FCcutoff = 0.25,
						  legendPosition = "none")

	# Display the plot
	print(p1)

	# Set output directory and create if it doesn't exist
	output_dir <- file.path(base_dir, "Figures", "Supplemental_figures")
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

	# Save the plot
	ggsave(file.path(output_dir, "volcano_progesterone_estrogen_overlap.pdf"),
		   p1, width = 9, height = 12)
	ggsave(file.path(output_dir, "volcano_progesterone_estrogen_overlap.png"),
		   p1, width = 9, height = 12, dpi = 300)

	cat("\n=== VOLCANO PLOT COMPLETE ===\n")
	cat("Plot saved to:", output_dir, "\n")
	cat("  - volcano_progesterone_estrogen_overlap.pdf\n")
	cat("  - volcano_progesterone_estrogen_overlap.png\n")

	# ============================================================================
	# GO ENRICHMENT ANALYSIS - DOWNREGULATED PATHWAYS
	# ============================================================================

	print("=== PERFORMING GO ENRICHMENT ANALYSIS ===")

	# Install required packages if needed (do this first, before loading)
	print("Checking and installing required packages...")

	if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
	  print("Installing clusterProfiler...")
	  BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
	}

	if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
	  print("Installing org.Mm.eg.db...")
	  BiocManager::install("org.Mm.eg.db", update = FALSE, ask = FALSE)
	}

	if (!requireNamespace("enrichplot", quietly = TRUE)) {
	  print("Installing enrichplot...")
	  BiocManager::install("enrichplot", update = FALSE, ask = FALSE)
	}

	if (!requireNamespace("stringr", quietly = TRUE)) {
	  print("Installing stringr...")
	  install.packages("stringr", repos = "https://cran.rstudio.com/", quiet = TRUE)
	}

	if (!requireNamespace("forcats", quietly = TRUE)) {
	  print("Installing forcats...")
	  install.packages("forcats", repos = "https://cran.rstudio.com/", quiet = TRUE)
	}

	# Now load all libraries after installations are complete
	print("Loading libraries...")
	loadLibrary(clusterProfiler)
	loadLibrary(org.Mm.eg.db)
	loadLibrary(enrichplot)
	loadLibrary(stringr)
	loadLibrary(forcats)

	# Extract downregulated and upregulated genes
	negs <- cluster0.markers[cluster0.markers$avg_log2FC < 0,]$Gene
	poses <- cluster0.markers[cluster0.markers$avg_log2FC > 0,]$Gene

	cat("Downregulated genes:", length(negs), "\n")
	cat("Upregulated genes:", length(poses), "\n")
	cat("Background genes:", length(keep.ind), "\n")

	# GO enrichment for downregulated genes
	print("Running GO enrichment for downregulated genes...")
	ggo.neg <- enrichGO(gene = negs,
						universe = keep.ind,
						OrgDb = org.Mm.eg.db,
						keyType = "SYMBOL", 
						ont = "BP",
						pAdjustMethod = "fdr",
						pvalueCutoff = 0.1,
						qvalueCutoff = 0.05)

	# GO enrichment for upregulated genes
	print("Running GO enrichment for upregulated genes...")
	ggo.pos <- enrichGO(gene = poses,
						universe = keep.ind,
						OrgDb = org.Mm.eg.db,
						keyType = "SYMBOL", 
						ont = "BP",
						pAdjustMethod = "fdr",
						pvalueCutoff = 0.1,
						qvalueCutoff = 0.05)

	# Save enrichment results
	write.csv(as.data.frame(ggo.neg), 
			  file.path(figure3_dir, "GO_enrichment_downregulated.csv"), 
			  row.names = FALSE)

	write.csv(as.data.frame(ggo.pos), 
			  file.path(figure3_dir, "GO_enrichment_upregulated.csv"), 
			  row.names = FALSE)

	# Print summary
	cat("\nDownregulated pathways found:", nrow(as.data.frame(ggo.neg)), "\n")
	cat("Upregulated pathways found:", nrow(as.data.frame(ggo.pos)), "\n")

	# Show top enriched terms
	if(nrow(as.data.frame(ggo.neg)) > 0) {
	  cat("\nTop 5 downregulated pathways:\n")
	  print(head(as.data.frame(ggo.neg)[,c("Description", "pvalue", "qvalue")], 5))
	}

	if(nrow(as.data.frame(ggo.pos)) > 0) {
	  cat("\nTop 5 upregulated pathways:\n")
	  print(head(as.data.frame(ggo.pos)[,c("Description", "pvalue", "qvalue")], 5))
	}

	# ============================================================================
	# CREATE CUSTOM GROUPED DOTPLOT FOR DOWNREGULATED PATHWAYS
	# ============================================================================

	print("\n=== CREATING CUSTOM GROUPED DOTPLOT ===")

	# Note: You should have a file called "GO downregulated pathways groupings .xlsx" in your Documents folder
	# This file should have 5 columns with headers as group names, and pathways listed under each column
	# Example structure:
	# | Immune response and interferon signaling | Mitochondrial respiration and ATP synthesis | ... |
	# | response to interferon-beta              | NADH dehydrogenase complex assembly         | ... |
	# | antigen processing...                    | mitochondrial respiratory chain...          | ... |

	loadLibrary(readxl)
	loadLibrary(tidyr)

	#PROBLEM
	#groupings_file <- "C:/Users/echowanec/Documents/GO downregulated pathways groupings .xlsx"
	groupings_file =paste0(base_dir,"/GO_downregulated_pathways_groupings.xlsx")
	if (!file.exists(groupings_file)){
		download.file("https://home.chpc.utah.edu/~u6004424/Chowanec_Data_2025/GO_downregulated_pathways_groupings.xlsx",groupings_file)
	}
	# Check if groupings file exists
	if(file.exists(groupings_file)) {
	  print("Found custom groupings file, using manual groupings...")
	  
	  # Read custom groupings (wide format with each group as a column)
	  groupings_wide <- read_excel(groupings_file)
	  
	  # Reshape from wide to long format
	  groupings <- groupings_wide %>%
		pivot_longer(cols = everything(),
					 names_to = "Group",
					 values_to = "Pathway") %>%
		filter(!is.na(Pathway))
	  
	  # Get GO data
	  go_data <- as.data.frame(ggo.neg)
	  go_data$score <- -log10(go_data$qvalue)
	  
	  # Merge with custom groupings
	  go_data_grouped <- merge(go_data, groupings, 
							   by.x = "Description", 
							   by.y = "Pathway",
							   all.x = FALSE)
	  
	  print(paste("Matched", nrow(go_data_grouped), "pathways to custom groups"))
	  
	} else {
	  print("Custom groupings file not found, using all downregulated pathways...")
	  print("To use custom groups, create an Excel file with pathway groupings")
	  
	  # Use all pathways without custom grouping
	  go_data_grouped <- as.data.frame(ggo.neg)
	  go_data_grouped$score <- -log10(go_data_grouped$qvalue)
	  go_data_grouped$Group <- "Downregulated"
	  go_data_grouped <- go_data_grouped %>% head(50)  # Top 50 pathways
	}

	# Plot order (bottom to top in the plot)
	group_order <- c("Sterol and lipid biosynthesis",
					 "Nucleotide and purine biosynthesis", 
					 "Mitochondrial respiration and ATP synthesis",
					 "Immune response and interferon signaling",
					 "Energy and nucleotide metabolism")

	# Only use groups that exist in the data
	group_order <- group_order[group_order %in% unique(go_data_grouped$Group)]

	go_data_grouped$Group <- factor(go_data_grouped$Group, levels = group_order)

	go_data_grouped <- go_data_grouped %>%
	  arrange(Group, desc(score))

	# Truncate pathway names to 70 characters
	go_data_grouped$Description <- str_sub(go_data_grouped$Description, 1, 70)

	# Handle duplicate names by adding asterisks
	while(any(duplicated(go_data_grouped$Description))) {
	  go_data_grouped$Description[duplicated(go_data_grouped$Description)] <-
		paste0(go_data_grouped$Description[duplicated(go_data_grouped$Description)], "*")
	}

	# Define specific colors for each group
	group_colors <- c("Sterol and lipid biosynthesis" = "#bebada",           # purple
					  "Nucleotide and purine biosynthesis" = "#fb8072",      # coral/red
					  "Mitochondrial respiration and ATP synthesis" = "#80b1d3",  # blue
					  "Immune response and interferon signaling" = "#fdb462",     # orange
					  "Energy and nucleotide metabolism" = "#b3de69",             # green
					  "Downregulated" = "#8dd3c7")  # Default color if no grouping

	# REVERSED legend order (top to bottom in legend matches plot top to bottom)
	legend_order <- rev(group_order)

	# Create the dotplot
	p_go_dotplot <- ggplot(go_data_grouped, 
						   aes(x = reorder(Description, as.numeric(Group)*1000 + score), 
							   y = score, 
							   color = Group)) +
	  geom_point(stat = "identity", size = 6) + 
	  ylim(0, max(go_data_grouped$score) + 2) +
	  scale_color_manual("",
						 values = group_colors,
						 breaks = legend_order) +
	  labs(title = "Downregulated Pathways: Progesterone vs Placebo",
		   y = "-log10(enrichment q-value)", 
		   x = "") +
	  coord_flip() + 
	  theme_bw() +
	  theme(axis.title.x = element_text(size = 20, hjust = 0.5),
			axis.text.y = element_text(size = 18, colour = "black"),
			axis.text.x = element_text(size = 20, colour = "black"),
			title = element_text(size = 16, hjust = 1),
			legend.text = element_text(size = 14),
			legend.direction = "vertical",
			legend.position = "right")

	# Display the plot
	print(p_go_dotplot)

	# Save the plot
	ggsave(file.path(figure3_dir, "downregulated_pathways_dotplot.pdf"),
		   p_go_dotplot, width = 16, height = 15)
	ggsave(file.path(figure3_dir, "downregulated_pathways_dotplot.png"),
		   p_go_dotplot, width = 16, height = 15, dpi = 300)

	# Save the processed data
	write.csv(go_data_grouped, 
			  file.path(figure3_dir, "processed_downregulated_pathways.csv"), 
			  row.names = FALSE)

	# Print summary
	cat("\n=== GO DOTPLOT SUMMARY ===\n")
	cat("Total pathways plotted:", nrow(go_data_grouped), "\n")

	if(length(group_order) > 1) {
	  cat("\nPathways per group:\n")
	  for(group in group_order) {
		n_pathways <- sum(go_data_grouped$Group == group, na.rm = TRUE)
		cat(group, ":", n_pathways, "pathways\n")
	  }
	}

	cat("\nGO enrichment analysis and dotplot complete!\n")
	cat("Files saved to:", figure3_dir, "\n")

	# ============================================================================
	# ATP/OXIDATIVE PHOSPHORYLATION HEATMAP
	# ============================================================================

	print("=== CREATING ATP/OXIDATIVE PHOSPHORYLATION HEATMAP (Fixed Labels) ===")

	loadLibrary(ComplexHeatmap)
	loadLibrary(circlize)
	loadLibrary(stringr)
	loadLibrary(grid)
	loadLibrary(dplyr)

	# Use existing data from previous run
	pathway_data <- as.data.frame(ggo.neg)

	# Define ATP/oxidative phosphorylation patterns
	atp_patterns <- c("ATP", "oxidative phosphorylation", "respiratory chain", 
					  "electron transport", "NADH", "cytochrome", "respiratory electron transport",
					  "mitochondrial ATP synthesis", "proton-transporting ATP synthase",
					  "energy derivation by oxidation", "cellular respiration",
					  "proton motive force")

	# Filter for ATP-related pathways
	atp_mask <- grepl(paste(atp_patterns, collapse = "|"), 
					  pathway_data$Description, ignore.case = TRUE)

	atp_pathways <- pathway_data[atp_mask, ]

	cat("\nFound", nrow(atp_pathways), "ATP/oxidative phosphorylation pathways\n")

	# Create pathway-gene matrix
	genelists <- atp_pathways$geneID
	genelists0 <- lapply(genelists, function(X) {
	  unlist(str_split(X, pattern = "/"))
	})
	allgenes <- Reduce(union, genelists0)
	genelists1 <- lapply(genelists0, function(X) {
	  as.numeric(allgenes %in% X)
	})

	genelists1 <- as.matrix(Reduce(rbind, genelists1))
	rownames(genelists1) <- atp_pathways$Description
	colnames(genelists1) <- allgenes

	# Get DE data
	de_atp <- cluster0.markers[cluster0.markers$Gene %in% colnames(genelists1), ]
	rownames(de_atp) <- de_atp$Gene
	de_atp$p_value <- -log10(de_atp$p_val_adj)

	annotation_data <- de_atp[colnames(genelists1), ]
	annotation_data <- as.data.frame(annotation_data[, c("p_value", "avg_log2FC")])

	# Select top genes if too many
	if(ncol(genelists1) > 50) {
	  gene_significance <- annotation_data$p_value
	  names(gene_significance) <- rownames(annotation_data)
	  top_genes <- names(sort(gene_significance, decreasing = TRUE))[1:50]
	  
	  genelists1 <- genelists1[, top_genes, drop = FALSE]
	  annotation_data <- annotation_data[top_genes, , drop = FALSE]
	}

	print(paste("Matrix dimensions:", nrow(genelists1), "×", ncol(genelists1)))

	# Color functions
	col_fun_pval <- colorRamp2(c(0, 4, 10, 100), 
							   c("#f7f7f7", "#d6604d", "#ca0020", "#67001f"))
	col_fun_fc <- colorRamp2(c(-1, -0.5, -0.25, 0), 
							 c("dodgerblue4", "dodgerblue", "lightskyblue1", "grey60"))

	# Create heatmap WITHOUT padding parameter in constructor
	ht_atp <- Heatmap(
	  genelists1,
	  width = unit(10, "inch"),
	  height = unit(4, "inch"),
	  row_title = "",
	  column_title = "ATP and Oxidative Phosphorylation",
	  clustering_distance_columns = "euclidean",
	  clustering_method_columns = "ward.D2",
	  col = colorRamp2(c(0, 1), c("#e0e0e0", "grey20")),
	  show_row_names = TRUE,
	  row_names_side = "left",
	  show_column_names = TRUE,
	  cluster_rows = FALSE,
	  show_row_dend = FALSE,
	  top_annotation = HeatmapAnnotation(
		"-log(p-value)" = annotation_data$p_value,
		"LogFC" = annotation_data$avg_log2FC,
		col = list(
		  "-log(p-value)" = col_fun_pval,
		  "LogFC" = col_fun_fc
		),
		show_legend = TRUE,
		show_annotation_name = TRUE,
		annotation_name_gp = gpar(fontsize = 10)
	  ),
	  show_column_dend = FALSE,
	  show_heatmap_legend = FALSE,
	  border_gp = gpar(col = "black", lty = 2),
	  rect_gp = gpar(col = "white", lwd = 1),
	  column_title_gp = gpar(fill = "white", col = "gray35", border = "white", fontsize = 14),
	  row_names_gp = gpar(fontsize = 11),
	  column_names_gp = gpar(fontsize = 9),
	  column_title_side = "top",
	  column_names_rot = 45,
	  heatmap_legend_param = list(direction = "horizontal")
	)

	# Save heatmap with wider dimensions and extra left padding in draw()
	pdf(file.path(figure3_dir, "ATP_oxidative_phosphorylation_heatmap.pdf"), 
		width = 20, height = 8)  # Even wider for pathway names

	# Draw with extra left padding
	draw(ht_atp, padding = unit(c(2, 2, 2, 80), "mm"))  # 80mm left padding for full pathway names

	dev.off()

	print(paste("\nHeatmap saved to:", 
				file.path(figure3_dir, "ATP_oxidative_phosphorylation_heatmap.pdf")))
	print("Pathway labels should now be fully visible!")

	cat("\n=== SUMMARY ===\n")
	cat("Pathways shown:", nrow(genelists1), "\n")
	cat("Genes shown:", ncol(genelists1), "\n")
	cat("Analysis complete!\n")

	# ============================================================================
	# CREATE CUSTOM GROUPED DOTPLOT FOR UPREGULATED PATHWAYS
	# ============================================================================

	print("\n=== CREATING CUSTOM GROUPED DOTPLOT FOR UPREGULATED PATHWAYS ===")

	# Try to find the upregulated groupings file in multiple locations
	# Note: trying both with and without space before .xlsx
	possible_paths_up <- c(
	  file.path(Sys.getenv("USERPROFILE"), "Documents", "GO upregulated pathways groupings.xlsx"),
	  file.path(Sys.getenv("USERPROFILE"), "Documents", "GO upregulated pathways groupings .xlsx"),
	  file.path(Sys.getenv("HOME"), "Documents", "GO upregulated pathways groupings.xlsx"),
	  file.path(figure3_dir, "GO upregulated pathways groupings.xlsx"),
	  file.path(getwd(), "GO upregulated pathways groupings.xlsx")
	)

	groupings_file_up <- NULL
	for(path in possible_paths_up) {
	  if(file.exists(path)) {
		groupings_file_up <- path
		print(paste("Found upregulated groupings file at:", groupings_file_up))
		break
	  }
	}

	# Check if groupings file exists
	if(!is.null(groupings_file_up) && file.exists(groupings_file_up)) {
	  print("Using custom pathway groupings from Excel file")
	  
	  # Read custom groupings (wide format with each group as a column)
	  groupings_wide_up <- read_excel(groupings_file_up)
	  
	  # Reshape from wide to long format
	  groupings_up <- groupings_wide_up %>%
		pivot_longer(cols = everything(),
					 names_to = "Group",
					 values_to = "Pathway") %>%
		filter(!is.na(Pathway))
	  
	  print("Sample of groupings:")
	  print(head(groupings_up, 10))
	  
	  # Get GO data for upregulated pathways
	  go_data_up <- as.data.frame(ggo.pos)
	  go_data_up$score <- -log10(go_data_up$qvalue)
	  
	  # Merge with custom groupings
	  go_data_grouped_up <- merge(go_data_up, groupings_up, 
								  by.x = "Description", 
								  by.y = "Pathway",
								  all.x = FALSE)
	  
	  print(paste("Matched", nrow(go_data_grouped_up), "pathways to custom groups"))
	  
	  if(nrow(go_data_grouped_up) == 0) {
		print("ERROR: No pathways matched between GO results and Excel file")
		print("Check that pathway names in Excel match GO term descriptions exactly")
	  } else {
		
		# Get unique group names
		group_names_up <- unique(groupings_up$Group)
		print(paste("Found", length(group_names_up), "groups:"))
		print(group_names_up)
		
		# Set group order (you can customize this)
		group_order_up <- group_names_up  # Use order from Excel columns
		
		go_data_grouped_up$Group <- factor(go_data_grouped_up$Group, levels = group_order_up)
		
		go_data_grouped_up <- go_data_grouped_up %>%
		  arrange(Group, desc(score))
		
		# Truncate pathway names to 70 characters
		go_data_grouped_up$Description <- str_sub(go_data_grouped_up$Description, 1, 70)
		
		# Handle duplicate names by adding asterisks
		while(any(duplicated(go_data_grouped_up$Description))) {
		  go_data_grouped_up$Description[duplicated(go_data_grouped_up$Description)] <-
			paste0(go_data_grouped_up$Description[duplicated(go_data_grouped_up$Description)], "*")
		}
		
		# Define colors for groups (Bright & Distinct palette)
		group_colors_up <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A")
		names(group_colors_up) <- group_order_up[1:min(length(group_order_up), 4)]
		
		# REVERSED legend order (top to bottom in legend matches plot top to bottom)
		legend_order_up <- rev(group_order_up)
		
		# Create the dotplot
		p_go_dotplot_up <- ggplot(go_data_grouped_up, 
								  aes(x = reorder(Description, as.numeric(Group)*1000 + score), 
									  y = score, 
									  color = Group)) +
		  geom_point(stat = "identity", size = 6) + 
		  ylim(0, max(go_data_grouped_up$score) + 2) +
		  scale_color_manual("",
							 values = group_colors_up,
							 breaks = legend_order_up) +
		  labs(title = "Upregulated Pathways: Progesterone vs Placebo",
			   y = "-log10(enrichment q-value)", 
			   x = "") +
		  coord_flip() + 
		  theme_bw() +
		  theme(axis.title.x = element_text(size = 20, hjust = 0.5),
				axis.text.y = element_text(size = 18, colour = "black"),
				axis.text.x = element_text(size = 20, colour = "black"),
				title = element_text(size = 16, hjust = 1),
				legend.text = element_text(size = 14),
				legend.direction = "vertical",
				legend.position = "right")
		
		# Display the plot
		print(p_go_dotplot_up)
		
		# Save the plot
		ggsave(file.path(figure3_dir, "upregulated_pathways_dotplot_custom_groups.pdf"),
			   p_go_dotplot_up, width = 16, height = 15)
		ggsave(file.path(figure3_dir, "upregulated_pathways_dotplot_custom_groups.png"),
			   p_go_dotplot_up, width = 16, height = 15, dpi = 300)
		
		# Save the processed data
		write.csv(go_data_grouped_up, 
				  file.path(figure3_dir, "processed_upregulated_pathways_custom_groups.csv"), 
				  row.names = FALSE)
		
		# Print summary
		cat("\n=== UPREGULATED PATHWAYS SUMMARY ===\n")
		cat("Total pathways plotted:", nrow(go_data_grouped_up), "\n\n")
		
		cat("Pathways per group:\n")
		for(group in group_order_up) {
		  n_pathways <- sum(go_data_grouped_up$Group == group, na.rm = TRUE)
		  if(n_pathways > 0) {
			cat(group, ":", n_pathways, "pathways\n")
		  }
		}
		
		cat("\nPlots saved to:", figure3_dir, "\n")
		cat("Upregulated pathways dotplot complete!\n")
	  }
	  
	} else {
	  print("Custom upregulated groupings file not found in the following locations:")
	  for(path in possible_paths_up) {
		print(paste("  -", path))
	  }
	  print("\nSkipping upregulated pathways dotplot")
	  print("To create this plot, add a file named 'GO upregulated pathways groupings.xlsx'")
	  print("with columns as group names and pathways listed below each group header")
	}


	# ============================================================================
	# SUPPLEMENTAL FIGURE: GENE EXPRESSION DOT PLOT
	# ============================================================================

	print("=== CREATING GENE EXPRESSION DOT PLOT ===")

	# Define genes of interest
	genes_to_plot <- c("Wnt4", "Myc", "Pgr", "Esr1", "Egfl7", "Col1a1", "Fcmr", 
					   "Klrb1a", "Cd3e", "Cd4", "Foxp3", "Cd8a", "Siglech", "Clec9a", 
					   "Adgre1", "Cd68", "Itgam", "Fcgr3", "Mki67")

	# Define custom order for y-axis (bottom to top) - ALL cell types from reloaded data
	cell_type_order <- c("Tumor 1", "Tumor 2", "Tumor 3", "Tumor 4", "Tumor 5", 
						 "Tumor 6", "Tumor 7", "Tumor 8", "Tumor 9", "Tumor 10", "Tumor 11",
						 "Epithelial 1", "Epithelial 2", "Epithelial 3", "Epithelial 4",
						 "Endothelial", 
						 "Fibroblast 1", "Fibroblast 2",
						 "Basal",
						 "B cell", "NK", 
						 "CD4+ 1", "CD4+ 2", "CD4+ Treg", 
						 "CD8+ 1", "CD8+ 2",
						 "Neutrophil",
						 "pDC", "cDC",
						 "Macrophage 1", "Macrophage 2", "Macrophage 3", 
						 "Macrophage 4", "Macrophage 5", "Macrophage 6", "Unknown")

	# Reverse the order so it displays bottom to top correctly
	cell_type_order <- rev(cell_type_order)

	# Set factor levels to control order BEFORE setting Idents
	SSM2sc$cell_type_secondary <- factor(SSM2sc$cell_type_secondary, 
										 levels = cell_type_order)

	SSM2sc$cell_type_secondary[is.na(SSM2sc$cell_type_secondary)] <- "Unknown"

	# Set identity to secondary cell type labels
	Idents(SSM2sc) <- "cell_type_secondary"

	# Create dot plot
	p_dotplot <- DotPlot(SSM2sc, 
						 features = genes_to_plot,
						 cols = c("lightgrey", "blue"),
						 dot.scale = 8) +
	  RotatedAxis() +
	  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
			axis.text.y = element_text(size = 12),
			axis.title = element_text(size = 14),
			legend.text = element_text(size = 10),
			legend.title = element_text(size = 12)) +
	  labs(title = "Gene Expression Across Secondary Cell Type Labels",
		   x = "Features",
		   y = "Secondary Cell Type Labels")

	# Display the plot
	print(p_dotplot)

	supplemental_dir <- file.path(base_dir, "Figures", "Supplemental_figures")
    if (!dir.exists(supplemental_dir)) dir.create(supplemental_dir, recursive = TRUE)
    if (!dir.exists(supplemental_dir)) {
        stop(paste("ERROR: Could not create directory:", supplemental_dir,"\nCheck that you have write permissions for this location."))
    }

	# Save the dot plot
	ggsave(file.path(supplemental_dir, "gene_expression_dotplot.pdf"),
		   p_dotplot, width = 12, height = 10)

	ggsave(file.path(supplemental_dir, "gene_expression_dotplot.png"),
		   p_dotplot, width = 12, height = 10, dpi = 300)

	print("=== SUPPLEMENTAL DOT PLOT COMPLETE ===")
	print(paste("Dot plot saved to:", file.path(supplemental_dir, "gene_expression_dotplot.png")))

}

Figure4 <- function() {

	# === FIGURE 4: TF ACTIVITY ANALYSIS: IMMUNE vs TUMOR (PLACEBO vs PROGESTERONE) ===

	library(decoupleR)
	library(dorothea)
	library(Seurat)
	library(dplyr)
	library(tidyr)
	library(tibble)

	Figure4_dir <- file.path(base_dir, "Figures", "Figure_4")
	if (!dir.exists(Figure4_dir)) dir.create(Figure4_dir, recursive = TRUE)
	
	# --- 1. Get DoRothEA TF-target network  ---
	net <- dorothea_mm %>%
	  filter(confidence %in% c("A", "B", "C")) %>%
	  rename(source = tf, mor = mor, target = target)
	
	# --- 2. Join layers on FULL object first, then subset ---
	SSM2sc[["RNA"]] <- JoinLayers(SSM2sc[["RNA"]])
	
	# Use cell_type_primary_final (not cell_type_primary)
	immune <- subset(SSM2sc, subset = cell_type_primary_final == "Immune")
	tumor <- subset(SSM2sc, subset = cell_type_primary_final == "Tumor")
	
	immune[["RNA"]] <- JoinLayers(immune[["RNA"]])
	tumor[["RNA"]] <- JoinLayers(tumor[["RNA"]])
	
	# --- 3. Run TF activity inference ---
	immune_acts <- run_wmean(
	  mat = as.matrix(immune[["RNA"]]$data),  
	  net = net,
	  .source = "source",
	  .target = "target",
	  .mor = "mor",
	  times = 100
	)
	
	tumor_acts <- run_wmean(
	  mat = as.matrix(tumor[["RNA"]]$data),
	  net = net,
	  .source = "source",
	  .target = "target",
	  .mor = "mor",
	  times = 100
	)
	
	# --- 4. Add TF activity scores as new assay ---
	immune_tf <- immune_acts %>% filter(statistic == "norm_wmean")
	tumor_tf <- tumor_acts %>% filter(statistic == "norm_wmean")
	
	immune_tf_mat <- immune_tf %>%
	  pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
	  column_to_rownames("condition") %>%
	  as.matrix() %>% t()
	
	tumor_tf_mat <- tumor_tf %>%
	  pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
	  column_to_rownames("condition") %>%
	  as.matrix() %>% t()
	
	immune[["tfact"]] <- CreateAssayObject(data = immune_tf_mat)
	tumor[["tfact"]] <- CreateAssayObject(data = tumor_tf_mat)
	
	# --- 5. Compare TF activity: Wilcoxon test (NOT FindMarkers) ---
	# FindMarkers fails on negative TF activity scores (log2FC produces NaN)
	# Wilcoxon test works directly on raw scores
	
	run_wilcox <- function(seurat_obj) {
	  tfs <- rownames(seurat_obj[["tfact"]])
	  scores <- FetchData(seurat_obj, vars = tfs, assay = "tfact")
	  scores$treatment <- seurat_obj$treatment
	  
	  results <- lapply(tfs, function(tf) {
	    w <- wilcox.test(scores[[tf]] ~ scores$treatment)
	    data.frame(
	      TF = tf,
	      mean_progesterone = mean(scores[[tf]][scores$treatment == "progesterone"]),
	      mean_placebo = mean(scores[[tf]][scores$treatment == "placebo"]),
	      diff = mean(scores[[tf]][scores$treatment == "progesterone"]) - 
	        mean(scores[[tf]][scores$treatment == "placebo"]),
	      p_val = w$p.value
	    )
	  }) %>% bind_rows()
	  
	  results$p_val_adj <- p.adjust(results$p_val, method = "BH")
	  results %>% arrange(p_val_adj)
	}
	
	immune_tf_results <- run_wilcox(immune)
	tumor_tf_results <- run_wilcox(tumor)
	
	# --- 6. Add DoRothEA confidence levels ---
	tf_confidence <- net %>%
	  group_by(source) %>%
	  summarise(best_confidence = min(confidence))
	
	immune_tf_results <- merge(immune_tf_results, tf_confidence, by.x = "TF", by.y = "source", all.x = TRUE)
	tumor_tf_results <- merge(tumor_tf_results, tf_confidence, by.x = "TF", by.y = "source", all.x = TRUE)
	
	# --- 7. Save ---
	write.csv(immune_tf_results, "immune_diff_TFs_wilcox.csv", row.names = FALSE)
	
	# Filter for significant TFs
	sig_immune <- immune_tf_results %>% filter(p_val_adj < 0.01) %>% pull(TF)
	sig_tumor <- tumor_tf_results %>% filter(p_val_adj < 0.01) %>% pull(TF)
	
	
	## IMMUNE AND TUMOR ONLY ## 
	
	immune_only2 <- setdiff(sig_immune, sig_tumor)
	tumor_only2 <- setdiff(sig_tumor, sig_immune)
	
	## DOTPLOT FOR TUMOR TFs - UNIQUE
	library(ggplot2)
	library(dplyr)
	
	# Calculate percent of cells expressing each TF
	calc_pct_expressing <- function(seurat_obj, tfs) {
	  scores <- FetchData(seurat_obj, vars = tfs, assay = "tfact")
	  sapply(tfs, function(tf) {
	    sum(scores[[tf]] > 0) / nrow(scores) * 100
	  })
	}
	
	# Get tumor-only TF results
	dot_df_tumor <- tumor_tf_results %>%
	  filter(TF %in% tumor_only2) %>%
	  filter(!TF %in% c("Mbd2", "Yy1")) %>%
	  dplyr::select(TF, diff, p_val_adj) %>%
	  mutate(
	    direction = ifelse(diff > 0, "Up", "Down"),
	    neg_log10_p = -log10(p_val_adj),
	    abs_diff = abs(diff)
	  )
	
	# Add percent expressing
	tumor_only2_plot <- tumor_only2[!tumor_only2 %in% c("Mbd2", "Yy1")]
	pct <- calc_pct_expressing(tumor, tumor_only2_plot)
	dot_df_tumor$pct_expressing <- pct[dot_df_tumor$TF]
	
	# Order TFs by p-value
	dot_df_tumor$TF <- factor(dot_df_tumor$TF, levels = dot_df_tumor$TF[order(dot_df_tumor$neg_log10_p)])
	dot_df_tumor$neg_log10_p[is.infinite(dot_df_tumor$neg_log10_p)] <- 350
	
	dotplot_tumor <- ggplot(dot_df_tumor, aes(x = neg_log10_p, y = TF)) +
	  geom_point(aes(size = pct_expressing, color = diff)) +
	  scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0,
	                        name = "Effect Size") +
	  scale_size_continuous(range = c(3, 10), name = "% Cells\nExpressing") +
	  theme_minimal() +
	  theme(axis.text.y = element_text(size = 16, face = "bold"),
	        axis.text.x = element_text(size = 16), 
	        axis.title.x = element_text(size = 14),
	        legend.title = element_text(size = 14),
	        legend.text = element_text(size = 14))+
	  xlab(expression(-Log[10]~"adjusted"~italic(p))) + ylab("") +
	  xlim(0, 375)
	
	print(dotplot_tumor)
	
	
	## DOTPLOT FOR IMMUNE TFs - UNIQUE
	dot_df_immune <- immune_tf_results %>%
	  filter(TF %in% immune_only2) %>%
	  filter(TF != "Mafg") %>% 
	  dplyr::select(TF, diff, p_val_adj) %>%
	  mutate(
	    direction = ifelse(diff > 0, "Up", "Down"),
	    neg_log10_p = -log10(p_val_adj),
	    abs_diff = abs(diff)
	    )
	
	# Add percent expressing
	immune_only2_plot <- immune_only2[immune_only2 != "Mafg"]
	pct_immune <- calc_pct_expressing(immune, immune_only2)
	dot_df_immune$pct_expressing <- pct_immune[dot_df_immune$TF]
	
	# Order TFs by p-value
	dot_df_immune$TF <- factor(dot_df_immune$TF, levels = dot_df_immune$TF[order(dot_df_immune$neg_log10_p)])
	dot_df_immune$neg_log10_p[is.infinite(dot_df_immune$neg_log10_p)] <- 350
	 
	  # After capping at 350, add small offsets
	  inf_idx <- which(dot_df_immune$neg_log10_p == 350)
	dot_df_immune$neg_log10_p[inf_idx] <- 350 + seq(0, length(inf_idx) - 1) * 5
	 
	  dotplot_immune <- ggplot(dot_df_immune, aes(x = neg_log10_p, y = TF)) +
	  geom_point(aes(size = pct_expressing, color = diff)) +
	  scale_color_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0,
	                              +                           name = "Effect Size") +
	  scale_size_continuous(range = c(3, 10), name = "% Cells\nExpressing") +
	  theme_minimal() +
	  theme(axis.text.y = element_text(size = 16, face = "bold"),
	              axis.text.x = element_text(size = 16), 
	              axis.title.x = element_text(size = 14),
	              legend.title = element_text(size = 14),
	              legend.text = element_text(size = 14))+
	  xlab(expression(-Log[10]~"adjusted"~italic(p))) + ylab("") +
	  xlim(0, 375)
	
	print(dotplot_immune)
	ggsave(file.path(Figure4_dir, "dotplot_immune_TF_largerfont.png"), dotplot_immune, width = 10, height = 8, dpi = 300)
  }
									 




Figure5 <- function(){
	print("Generating Figure5")
		getData("SSM2sc_with_celltypes.RDS","SSM2sc")


	# ============================================================================
	# Figure 5: T CELL PERCENTAGE OF IMMUNE CELLS BY TREATMENT GROUP
	# ============================================================================

	# Create Figure 5 subdirectory
	Figure5_dir <- file.path(base_dir, "Figures", "Figure_5")
	if (!dir.exists(Figure5_dir)) dir.create(Figure5_dir, recursive = TRUE)

	# Verify the directory was created successfully
	if (!dir.exists(Figure5_dir)) {
	  stop(paste("ERROR: Could not create directory:", Figure5_dir,
				 "\nCheck that you have write permissions for this location."))
	}

	supplemental_dir <- file.path(base_dir, "Figures", "Supplemental_figures")
	if (!dir.exists(supplemental_dir)) dir.create(supplemental_dir, recursive = TRUE)
	if (!dir.exists(supplemental_dir)) {
		stop(paste("ERROR: Could not create directory:", supplemental_dir,"\nCheck that you have write permissions for this location."))
	}



	print("=== Figure 5: T CELL PERCENTAGE OF IMMUNE CELLS BY TREATMENT ===")

	# Rename T cell labels in metadata for clearer visualization
	print("Renaming T cell subtypes in metadata...")
	SSM2sc$cell_type_secondary <- as.character(SSM2sc$cell_type_secondary)
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "CD4+ 1"] <- "CD4+ Naive"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "CD4+ 2"] <- "CD4+ CM"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "CD8+ 1"] <- "CD8+ CM"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "CD8+ 2"] <- "CD8+ TRM"

	# Save the updated Seurat object with renamed T cells
	saveRDS(SSM2sc, file = file.path(base_dir, "SSM2sc_with_celltypes.RDS"))
	print("Updated Seurat object saved with renamed T cell labels")

	# Subset to immune cells only
	SSM2sc_immune <- subset(SSM2sc, subset = cell_type_primary_final == "Immune")
	print(paste("Total immune cells:", ncol(SSM2sc_immune)))

	# Define T cell types (now using renamed labels)
	tcell_types <- c("CD4+ Naive", "CD4+ CM", "CD4+ Treg", "CD8+ CM", "CD8+ TRM")

	# Subset to T cells only
	SSM2sc_tcells <- subset(SSM2sc, subset = cell_type_secondary %in% tcell_types)
	print(paste("Total T cells:", ncol(SSM2sc_tcells)))
	print("T cells per treatment:")
	print(table(SSM2sc_tcells$treatment))

	# Calculate total immune cells per treatment
	immune_counts <- table(SSM2sc_immune$treatment)

	# Count T cells by treatment
	tcell_counts <- table(SSM2sc_tcells$treatment)

	# Create data frame with counts and percentages
	count_df <- data.frame(
	  Treatment = names(tcell_counts),
	  Count = as.numeric(tcell_counts),
	  Total_Immune = as.numeric(immune_counts[names(tcell_counts)]),
	  stringsAsFactors = FALSE
	)

	# Calculate percentage
	count_df$Percentage <- (count_df$Count / count_df$Total_Immune) * 100

	# Capitalize treatment names for plot
	count_df$Treatment <- factor(count_df$Treatment, 
								 levels = c("placebo", "progesterone"),
								 labels = c("Placebo", "Progesterone"))

	# Print summary
	print("T cell percentages:")
	print(count_df)

	# Define treatment colors
	treatment_colors <- c("Placebo" = "#4292C6", "Progesterone" = "#EF3B2C")

	# Create bar plot of T cell percentages
	p_tcell_percent <- ggplot(count_df, aes(x = Treatment, y = Percentage, fill = Treatment)) +
	  geom_bar(stat = "identity", width = 0.7) +
	  geom_text(aes(label = Count), vjust = -0.5, size = 8) +
	  scale_fill_manual(values = treatment_colors) +
	  labs(title = "T Cells as % of Total Immune Cells",
		   x = "Treatment",
		   y = "Percentage of Immune Cells (%)") +
	  theme_classic() +
	  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
			axis.title.x = element_text(size = 16),
			axis.title.y = element_text(size = 16),
			axis.text.x = element_text(size = 14),
			axis.text.y = element_text(size = 14),
			legend.title = element_text(size = 14),
			legend.text = element_text(size = 12),
			legend.position = "right") +
	  ylim(0, max(count_df$Percentage) * 1.15)

	print(p_tcell_percent)

	# Save plot
	ggsave(file.path(Figure5_dir, "tcell_percentage_by_treatment.pdf"),
		   p_tcell_percent, width = 8, height = 6)
	ggsave(file.path(Figure5_dir, "tcell_percentage_by_treatment.png"),
		   p_tcell_percent, width = 8, height = 6, dpi = 300)

	# Save the data table
	write.csv(count_df, 
			  file.path(Figure5_dir, "tcell_percentage_data.csv"), 
			  row.names = FALSE)

	print("=== Figure 5 COMPLETE ===")
	print(paste("T cell percentage bar plot saved to:", 
				file.path(Figure5_dir, "tcell_percentage_by_treatment.png")))

	# ============================================================================
	# SUPPLEMENTAL FIGURE: T CELL AND NK CELL MARKER GENES
	# ============================================================================

	print("=== T CELL AND NK CELL MARKER GENE ANALYSIS ===")

	# Define T cell and NK cell types from cell_type_secondary
	tcell_nk_types <- c("CD4+ 1", "CD4+ 2", "CD4+ Treg", "CD8+ 1", "CD8+ 2", "NK")

	print("Filtering for T cells and NK cells...")
	print(paste("Cell types included:", paste(tcell_nk_types, collapse = ", ")))

	# Subset to T cells and NK cells
	SSM2sc_tcells_nk <- subset(SSM2sc, subset = cell_type_secondary %in% tcell_nk_types)

	print(paste("Filtered from", ncol(SSM2sc), "to", ncol(SSM2sc_tcells_nk), "T cells and NK cells"))
	print("Cell counts by type:")
	print(table(SSM2sc_tcells_nk$cell_type_secondary))

	# Create renamed labels for plotting
	SSM2sc_tcells_nk$cell_type_renamed <- SSM2sc_tcells_nk$cell_type_secondary
	levels(SSM2sc_tcells_nk$cell_type_renamed) <- c(
	  "CD4+ 1" = "CD4+ Naive",
	  "CD4+ 2" = "CD4+ CM",
	  "CD4+ Treg" = "CD4+ Treg",
	  "CD8+ 1" = "CD8+ CM",
	  "CD8+ 2" = "CD8+ TRM",
	  "NK" = "NK"
	)

	# Map the values
	SSM2sc_tcells_nk$cell_type_renamed <- factor(
	  SSM2sc_tcells_nk$cell_type_secondary,
	  levels = c("CD4+ 1", "CD4+ 2", "CD4+ Treg", "CD8+ 1", "CD8+ 2", "NK"),
	  labels = c("CD4+ Naive", "CD4+ CM", "CD4+ Treg", "CD8+ CM", "CD8+ TRM", "NK")
	)

	# Filter mitochondrial and ribosomal genes
	print("Filtering genes...")
	rm.ind <- grep("^mt-|^Rp[ls]", rownames(SSM2sc_tcells_nk))
	if(length(rm.ind) > 0) {
	  keep <- rownames(SSM2sc_tcells_nk)[-rm.ind]
	  SSM2sc_tcells_nk <- SSM2sc_tcells_nk[keep,]
	  print(paste("Removed", length(rm.ind), "mitochondrial/ribosomal genes"))
	}

	# Filter lowly expressed genes (present in <20 cells)
	counts_data <- GetAssayData(SSM2sc_tcells_nk, assay = "RNA", layer = "counts")
	keep.ind <- rowSums(counts_data > 0)
	keep.ind <- names(keep.ind[keep.ind >= 20])
	SSM2sc_tcells_nk <- SSM2sc_tcells_nk[keep.ind,]
	print(paste("Kept", length(keep.ind), "genes expressed in >=20 cells"))

	# Set up cell type order (using original names for FindMarkers)
	order_vec <- c("CD4+ 1", "CD4+ 2", "CD4+ Treg", "CD8+ 1", "CD8+ 2", "NK")

	# Set factor levels for proper ordering
	SSM2sc_tcells_nk$cell_type_secondary <- factor(SSM2sc_tcells_nk$cell_type_secondary, 
												   levels = order_vec)
	Idents(SSM2sc_tcells_nk) <- "cell_type_secondary"

	# Find marker genes for each cell type
	print("Finding marker genes for each cell type...")
	main_list <- data.frame()

	for (cell_type in order_vec) {
	  print(paste("Processing:", cell_type))
	  
	  # Check if this cell type exists in the data
	  if(sum(SSM2sc_tcells_nk$cell_type_secondary == cell_type) > 10) {
		
		tryCatch({
		  # Find markers for this cell type vs all others
		  cluster_markers <- FindMarkers(SSM2sc_tcells_nk, 
										 ident.1 = cell_type,
										 only.pos = TRUE,
										 min.pct = 0.25,
										 logfc.threshold = 0.25,
										 test.use = "wilcox",
										 verbose = FALSE)
		  
		  # Filter for positive fold changes
		  cluster_markers <- cluster_markers[cluster_markers$avg_log2FC > 0, ]
		  
		  if(nrow(cluster_markers) > 0) {
			# Apply FDR correction
			cluster_markers$p_val_fdr <- p.adjust(cluster_markers$p_val, method = "fdr")
			
			# Order by FDR-corrected p-values, then by fold change
			cluster_markers <- cluster_markers[order(cluster_markers$p_val_fdr, 
													 -cluster_markers$avg_log2FC), ]
			
			# Filter for FDR significance (FDR p < 0.05)
			significant_markers <- cluster_markers[cluster_markers$p_val_fdr < 0.05, ]
			
			if(nrow(significant_markers) > 0) {
			  # Take top 10 FDR-significant markers
			  n_markers <- min(10, nrow(significant_markers))
			  cluster_markers <- significant_markers[1:n_markers, ]
			  
			  cluster_markers$Gene <- rownames(cluster_markers)
			  cluster_markers$cell <- cell_type
			  
			  main_list <- rbind(main_list, cluster_markers)
			  print(paste("  Found", n_markers, "FDR-significant markers"))
			} else {
			  print(paste("  No FDR-significant markers found for", cell_type))
			}
		  }
		}, error = function(e) {
		  print(paste("  Error processing", cell_type, ":", e$message))
		})
	  } else {
		print(paste("  Skipping", cell_type, "- too few cells"))
	  }
	}

	# Create outputs if markers were found
	if(nrow(main_list) > 0) {
	  print(paste("Found", nrow(main_list), "total FDR-significant marker genes"))
	  
	  # Create dot plot with renamed labels
	  print("Creating marker gene dot plot...")
	  
	  # Set identity to renamed cell types for plotting
	  Idents(SSM2sc_tcells_nk) <- "cell_type_renamed"
	  
	  p_tcell_markers <- DotPlot(SSM2sc_tcells_nk, 
								 features = unique(main_list$Gene),
								 cols = c(low = "#f7f7f7", high = "#6e016b"),
								 dot.min = 0.1) +
		scale_fill_gradient2(low = "#2166ac",
							 mid = "#f7f7f7",
							 high = "#b2182b") +
		coord_flip() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
			  axis.text.y = element_text(size = 10),
			  axis.title = element_blank(),
			  plot.title = element_text(hjust = 0.5)) +
		ggtitle("T Cell and NK Cell Marker Genes")
	  
	  # Display plot
	  print(p_tcell_markers)
	  
	  # Save outputs
	  ggsave(file.path(supplemental_dir, "tcell_nk_marker_dotplot.pdf"), 
			 p_tcell_markers, width = 10, height = 14)
	  ggsave(file.path(supplemental_dir, "tcell_nk_marker_dotplot.png"), 
			 p_tcell_markers, width = 10, height = 14, dpi = 300)
	  write.csv(main_list, 
				file.path(supplemental_dir, "tcell_nk_top_markers.csv"), 
				row.names = FALSE)
	  
	  print("=== T CELL/NK MARKER ANALYSIS COMPLETE ===")
	  print(paste("Marker dot plot saved to:", file.path(supplemental_dir, "tcell_nk_marker_dotplot.png")))
	  print(paste("Marker gene list saved to:", file.path(supplemental_dir, "tcell_nk_top_markers.csv")))
	  
	} else {
	  print("No FDR-significant marker genes found for T cells/NK cells")
	}

	# ============================================================================
	# SUPPLEMENTAL FIGURE: NAIVE AND MEMORY MARKER GENES
	# ============================================================================

	print("=== CREATING NAIVE AND MEMORY MARKER GENES VIOLIN PLOT ===")

	# Define genes of interest
	migration_genes <- c("Sell", "Ccr7", "S1pr1", "Itga1", "Cdh1")

	# Define T cell types only (no NK)
	tcell_types <- c("CD4+ Naive", "CD4+ CM", "CD4+ Treg", "CD8+ CM", "CD8+ TRM")

	# Subset to T cells only
	SSM2sc_tcells_violin <- subset(SSM2sc, subset = cell_type_secondary %in% tcell_types)

	# Define colors for each cluster
	tcell_colors_violin <- c(
	  "CD4+ Naive" = "#F8A19F",
	  "CD4+ CM" = "#FA0087", 
	  "CD4+ Treg" = "#5A5156",
	  "CD8+ CM" = "#2ED9FF",
	  "CD8+ TRM" = "#E4E1E3"
	)

	# Set identity
	Idents(SSM2sc_tcells_violin) <- "cell_type_secondary"

	# Create individual violin plots and combine them
	plot_list <- list()
	for(gene in migration_genes) {
	  p <- VlnPlot(SSM2sc_tcells_violin, 
				   features = gene,
				   pt.size = 0,
				   cols = tcell_colors_violin) +
		theme(legend.position = "none",
			  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
			  axis.title.x = element_blank(),
			  plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
	  
	  plot_list[[gene]] <- p
	}

	# Combine plots vertically
	p_migration_violin <- wrap_plots(plot_list, ncol = 1)

	# Display plot
	print(p_migration_violin)

	# Save plot to supplemental figures directory
	ggsave(file.path(supplemental_dir, "naive_memory_marker_genes_violin.pdf"), 
		   p_migration_violin, width = 10, height = 14)
	ggsave(file.path(supplemental_dir, "naive_memory_marker_genes_violin.png"), 
		   p_migration_violin, width = 10, height = 14, dpi = 300)

	print("=== SUPPLEMENTAL NAIVE AND MEMORY MARKER GENES COMPLETE ===")
	print(paste("Violin plot saved to:", file.path(supplemental_dir, "naive_memory_marker_genes_violin.png")))


	# ============================================================================
	# Figure 5: T CELL PROPORTION STACKED BAR PLOT
	# ============================================================================

	print("Creating T cell proportion stacked bar plot...")

	# Define T cell and NK types
	tcell_nk_types <- c("CD4+ Naive", "CD4+ CM", "CD4+ Treg", "CD8+ CM", "CD8+ TRM", "NK")

	# Extract metadata and filter for specific T cell and NK clusters
	subset_data <- SSM2sc@meta.data[SSM2sc@meta.data$cell_type_secondary %in% tcell_nk_types, ]

	# Create summary table with proportions by treatment
	cluster_treatment_counts <- subset_data %>%
	  group_by(treatment, cell_type_secondary) %>%
	  summarise(count = n(), .groups = 'drop') %>%
	  group_by(treatment) %>%
	  mutate(proportion = count / sum(count))

	# Define custom colors for T cell subtypes and NK cells
	tcell_colors <- c(
	  "CD4+ Naive" = "#F8A19F",
	  "CD4+ CM" = "#FA0087", 
	  "CD4+ Treg" = "#5A5156",
	  "CD8+ CM" = "#2ED9FF",
	  "CD8+ TRM" = "#E4E1E3",
	  "NK" = "#FEAF16"
	)

	# Create stacked bar plot
	p_tcell_proportion <- ggplot(cluster_treatment_counts, 
								 aes(x = treatment, y = proportion, fill = cell_type_secondary)) +
	  geom_bar(stat = "identity", position = "stack", width = 0.6) +
	  scale_fill_manual(values = tcell_colors) +
	  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00), 
						 labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
	  labs(x = "", y = "", fill = "") +
	  theme_minimal() +
	  theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "white", color = NA),
		plot.background = element_rect(fill = "white", color = NA),
		axis.text.x = element_text(size = 14, color = "black"),
		axis.text.y = element_text(size = 12, color = "black"),
		axis.ticks = element_line(color = "black"),
		legend.position = "right",
		legend.text = element_text(size = 11),
		legend.margin = margin(l = 20),
		plot.margin = margin(t = 10, r = 20, b = 10, l = 50)
	  ) +
	  coord_flip()

	print(p_tcell_proportion)
	ggsave(file.path(Figure5_dir, "tcell_proportion_stacked.pdf"), 
		   p_tcell_proportion, width = 10, height = 4)
	ggsave(file.path(Figure5_dir, "tcell_proportion_stacked.png"), 
		   p_tcell_proportion, width = 10, height = 4, dpi = 300)

	# ============================================================================
	# Figure 5: STATISTICAL COMPARISON OF T CELL PROPORTIONS
	# ============================================================================

	print("Performing statistical analysis of T cell and NK proportions...")

	# Get total cells per treatment
	total_placebo <- sum(subset_data$treatment == "placebo")
	total_progesterone <- sum(subset_data$treatment == "progesterone")

	cat("\nTotal T cells and NK cells per treatment:\n")
	cat("Placebo:", total_placebo, "\n")
	cat("Progesterone:", total_progesterone, "\n")

	# Initialize results dataframe
	results <- data.frame(
	  cell_type = character(),
	  placebo_count = numeric(),
	  progesterone_count = numeric(),
	  placebo_prop = numeric(),
	  progesterone_prop = numeric(),
	  fold_change = numeric(),
	  p_value = numeric(),
	  stringsAsFactors = FALSE
	)

	# Perform Fisher's exact test for each cell type
	for(cell_type in tcell_nk_types) {
	  
	  # Get counts for this cell type
	  placebo_this_type <- sum(subset_data$treatment == "placebo" & 
								 subset_data$cell_type_secondary == cell_type)
	  progesterone_this_type <- sum(subset_data$treatment == "progesterone" & 
									  subset_data$cell_type_secondary == cell_type)
	  
	  # Get counts for all other T cell types combined
	  placebo_other_types <- total_placebo - placebo_this_type
	  progesterone_other_types <- total_progesterone - progesterone_this_type
	  
	  # Create contingency table
	  contingency_table <- matrix(c(placebo_this_type, placebo_other_types,
									progesterone_this_type, progesterone_other_types),
								  nrow = 2, byrow = TRUE)
	  
	  # Perform Fisher's exact test
	  fisher_test <- fisher.test(contingency_table)
	  
	  # Calculate proportions
	  placebo_prop <- placebo_this_type / total_placebo
	  progesterone_prop <- progesterone_this_type / total_progesterone
	  
	  # Calculate fold change
	  fold_change <- ifelse(placebo_prop > 0, progesterone_prop / placebo_prop, NA)
	  
	  # Add to results
	  results <- rbind(results, data.frame(
		cell_type = cell_type,
		placebo_count = placebo_this_type,
		progesterone_count = progesterone_this_type,
		placebo_prop = placebo_prop,
		progesterone_prop = progesterone_prop,
		fold_change = fold_change,
		p_value = fisher_test$p.value
	  ))
	}

	# Apply FDR correction
	results$p_adj <- p.adjust(results$p_value, method = "fdr")

	# Add significance indicators
	results$significance <- ifelse(results$p_adj < 0.001, "***",
								   ifelse(results$p_adj < 0.01, "**",
										  ifelse(results$p_adj < 0.05, "*", "ns")))

	# Sort by adjusted p-value
	results <- results[order(results$p_adj), ]

	# Save results as CSV
	write.csv(results, file.path(Figure5_dir, "tcell_statistical_comparison.csv"), row.names = FALSE)

	# Save results as Excel file (requires writexl or openxlsx package)
	if(require(writexl, quietly = TRUE)) {
	  write_xlsx(results, file.path(Figure5_dir, "tcell_statistical_comparison.xlsx"))
	  print("Statistical results saved as Excel file")
	} else if(require(openxlsx, quietly = TRUE)) {
	  write.xlsx(results, file.path(Figure5_dir, "tcell_statistical_comparison.xlsx"))
	  print("Statistical results saved as Excel file")
	} else {
	  print("Note: Install 'writexl' or 'openxlsx' package to save as Excel file")
	}

	cat("\n=== STATISTICAL RESULTS ===\n")
	print(results)

	# ============================================================================
	# Figure 5: T CELL UMAP
	# ============================================================================

	print("Creating T cell and NK UMAP...")

	# Subset to T cells and NK cells
	SSM2sc_tcells <- subset(SSM2sc, subset = cell_type_secondary %in% tcell_nk_types)

	# Define colors including NK
	tcell_colors <- c(
	  "CD4+ Naive" = "#F8A19F",
	  "CD4+ CM" = "#FA0087", 
	  "CD4+ Treg" = "#5A5156",
	  "CD8+ CM" = "#2ED9FF",
	  "CD8+ TRM" = "#E4E1E3",
	  "NK" = "#FEAF16"
	)

	# Create clean UMAP plot (all cells)
	p_tcell_umap <- DimPlot(SSM2sc_tcells, 
							reduction = "umap", 
							group.by = "cell_type_secondary",
							cols = tcell_colors,
							pt.size = 0.5) +
	  labs(title = "", x = "", y = "") +
	  theme_void() +
	  theme(
		legend.title = element_blank(),
		legend.text = element_text(size = 12),
		legend.position = "right",
		plot.background = element_rect(fill = "white", color = NA),
		panel.background = element_rect(fill = "white", color = NA)
	  ) +
	  guides(color = guide_legend(override.aes = list(size = 4)))

	print(p_tcell_umap)
	ggsave(file.path(Figure5_dir, "tcell_umap.pdf"), 
		   p_tcell_umap, width = 10, height = 8)
	ggsave(file.path(Figure5_dir, "tcell_umap.png"), 
		   p_tcell_umap, width = 10, height = 8, dpi = 300)

	# Create UMAP for placebo treatment only
	print("Creating placebo T cell and NK UMAP...")
	SSM2sc_tcells_placebo <- subset(SSM2sc_tcells, subset = treatment == "placebo")

	p_tcell_umap_placebo <- DimPlot(SSM2sc_tcells_placebo, 
									reduction = "umap", 
									group.by = "cell_type_secondary",
									cols = tcell_colors,
									pt.size = 0.5) +
	  labs(title = "Placebo", x = "", y = "") +
	  theme_void() +
	  theme(
		legend.title = element_blank(),
		legend.text = element_text(size = 12),
		legend.position = "right",
		plot.background = element_rect(fill = "white", color = NA),
		panel.background = element_rect(fill = "white", color = NA),
		plot.title = element_text(size = 16, hjust = 0.5)
	  ) +
	  guides(color = guide_legend(override.aes = list(size = 4)))

	print(p_tcell_umap_placebo)
	ggsave(file.path(Figure5_dir, "tcell_umap_placebo.pdf"), 
		   p_tcell_umap_placebo, width = 10, height = 8)
	ggsave(file.path(Figure5_dir, "tcell_umap_placebo.png"), 
		   p_tcell_umap_placebo, width = 10, height = 8, dpi = 300)

	# Create UMAP for progesterone treatment only
	print("Creating progesterone T cell and NK UMAP...")
	SSM2sc_tcells_progesterone <- subset(SSM2sc_tcells, subset = treatment == "progesterone")

	p_tcell_umap_progesterone <- DimPlot(SSM2sc_tcells_progesterone, 
										 reduction = "umap", 
										 group.by = "cell_type_secondary",
										 cols = tcell_colors,
										 pt.size = 0.5) +
	  labs(title = "Progesterone", x = "", y = "") +
	  theme_void() +
	  theme(
		legend.title = element_blank(),
		legend.text = element_text(size = 12),
		legend.position = "right",
		plot.background = element_rect(fill = "white", color = NA),
		panel.background = element_rect(fill = "white", color = NA),
		plot.title = element_text(size = 16, hjust = 0.5)
	  ) +
	  guides(color = guide_legend(override.aes = list(size = 4)))

	print(p_tcell_umap_progesterone)
	ggsave(file.path(Figure5_dir, "tcell_umap_progesterone.pdf"), 
		   p_tcell_umap_progesterone, width = 10, height = 8)
	ggsave(file.path(Figure5_dir, "tcell_umap_progesterone.png"), 
		   p_tcell_umap_progesterone, width = 10, height = 8, dpi = 300)

}

Figure6 <- function(){
	print("Generating Figure6")
		getData("SSM2sc_with_celltypes.RDS","SSM2sc")




	# ============================================================================
	# Figure 6: NEUTROPHIL PERCENTAGE OF IMMUNE CELLS BY TREATMENT
	# ============================================================================
	supplemental_dir <- file.path(base_dir, "Figures", "Supplemental_figures")
	if (!dir.exists(supplemental_dir)) dir.create(supplemental_dir, recursive = TRUE)
	
	# Create Figure 6 subdirectory
	Figure6_dir <- file.path(base_dir, "Figures", "Figure_6")
	if (!dir.exists(Figure6_dir)) dir.create(Figure6_dir, recursive = TRUE)

	# Verify the directory was created successfully
	if (!dir.exists(Figure6_dir)) {
	  stop(paste("ERROR: Could not create directory:", Figure6_dir,
				 "\nCheck that you have write permissions for this location."))
	}

	print("=== Figure 6: NEUTROPHIL PERCENTAGE OF IMMUNE CELLS BY TREATMENT ===")

	# Subset to immune cells only
	SSM2sc_immune <- subset(SSM2sc, subset = cell_type_primary_final == "Immune")
	print(paste("Total immune cells:", ncol(SSM2sc_immune)))

	# Subset to neutrophils only
	SSM2sc_neutrophils <- subset(SSM2sc, subset = cell_type_secondary == "Neutrophil")
	print(paste("Total neutrophils:", ncol(SSM2sc_neutrophils)))
	print("Neutrophils per treatment:")
	print(table(SSM2sc_neutrophils$treatment))

	# Calculate total immune cells per treatment
	immune_counts <- table(SSM2sc_immune$treatment)

	# Count neutrophils by treatment
	neutrophil_counts <- table(SSM2sc_neutrophils$treatment)

	# Create data frame with counts and percentages
	count_df_neut <- data.frame(
	  Treatment = names(neutrophil_counts),
	  Count = as.numeric(neutrophil_counts),
	  Total_Immune = as.numeric(immune_counts[names(neutrophil_counts)]),
	  stringsAsFactors = FALSE
	)

	# Calculate percentage
	count_df_neut$Percentage <- (count_df_neut$Count / count_df_neut$Total_Immune) * 100

	# Capitalize treatment names for plot
	count_df_neut$Treatment <- factor(count_df_neut$Treatment, 
									  levels = c("placebo", "progesterone"),
									  labels = c("Placebo", "Progesterone"))

	# Print summary
	print("Neutrophil percentages:")
	print(count_df_neut)

	# Define treatment colors (same as T cell plot)
	treatment_colors <- c("Placebo" = "#4292C6", "Progesterone" = "#EF3B2C")

	# Create bar plot of neutrophil percentages
	p_neutrophil_percent <- ggplot(count_df_neut, aes(x = Treatment, y = Percentage, fill = Treatment)) +
	  geom_bar(stat = "identity", width = 0.7) +
	  geom_text(aes(label = Count), vjust = -0.5, size = 8) +
	  scale_fill_manual(values = treatment_colors) +
	  labs(title = "Neutrophils as % of Total Immune Cells",
		   x = "Treatment",
		   y = "Percentage of Immune Cells (%)") +
	  theme_classic() +
	  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
			axis.title.x = element_text(size = 16),
			axis.title.y = element_text(size = 16),
			axis.text.x = element_text(size = 14),
			axis.text.y = element_text(size = 14),
			legend.title = element_text(size = 14),
			legend.text = element_text(size = 12),
			legend.position = "right") +
	  ylim(0, max(count_df_neut$Percentage) * 1.15)

	print(p_neutrophil_percent)

	# Save plot
	ggsave(file.path(Figure6_dir, "neutrophil_percentage_by_treatment.pdf"),
		   p_neutrophil_percent, width = 8, height = 6)
	ggsave(file.path(Figure6_dir, "neutrophil_percentage_by_treatment.png"),
		   p_neutrophil_percent, width = 8, height = 6, dpi = 300)

	# Save the data table
	write.csv(count_df_neut, 
			  file.path(Figure6_dir, "neutrophil_percentage_data.csv"), 
			  row.names = FALSE)

	print("=== Figure 6 COMPLETE ===")
	print(paste("Neutrophil percentage bar plot saved to:", 
				file.path(Figure6_dir, "neutrophil_percentage_by_treatment.png")))
	# ============================================================================
	# Figure 6: MACROPHAGE PERCENTAGE OF IMMUNE CELLS BY TREATMENT
	# ============================================================================

	print("=== CREATING MACROPHAGE PERCENTAGE BAR PLOT ===")

	# Subset to immune cells only
	SSM2sc_immune <- subset(SSM2sc, subset = cell_type_primary_final == "Immune")
	print(paste("Total immune cells:", ncol(SSM2sc_immune)))


	# Define macrophage types
	macrophage_types <- c("Macrophage 1", "Macrophage 2", "Macrophage 3",
                      "Macrophage 4", "Macrophage 5", "Macrophage 6")

	# Subset to macrophages only
	SSM2sc_macrophages <- subset(SSM2sc, subset = cell_type_secondary %in% macrophage_types)

	print(paste("Total macrophages:", ncol(SSM2sc_macrophages)))
	print("Macrophages per treatment:")
	print(table(SSM2sc_macrophages$treatment))

	# Calculate total immune cells per treatment
	immune_counts <- table(SSM2sc_immune$treatment)

	# Count macrophages by treatment
	macrophage_counts <- table(SSM2sc_macrophages$treatment)

	# Create data frame with counts and percentages
	count_df_macro <- data.frame(
	  Treatment = names(macrophage_counts),
	  Count = as.numeric(macrophage_counts),
	  Total_Immune = as.numeric(immune_counts[names(macrophage_counts)]),
	  stringsAsFactors = FALSE
	)

	# Calculate percentage
	count_df_macro$Percentage <- (count_df_macro$Count / count_df_macro$Total_Immune) * 100

	# Capitalize treatment names for plot
	count_df_macro$Treatment <- factor(count_df_macro$Treatment, 
									   levels = c("placebo", "progesterone"),
									   labels = c("Placebo", "Progesterone"))

	# Print summary
	print("Macrophage percentages:")
	print(count_df_macro)

	# Define treatment colors (same as previous plots)
	treatment_colors <- c("Placebo" = "#4292C6", "Progesterone" = "#EF3B2C")

	# Create bar plot of macrophage percentages
	p_macrophage_percent <- ggplot(count_df_macro, aes(x = Treatment, y = Percentage, fill = Treatment)) +
	  geom_bar(stat = "identity", width = 0.7) +
	  geom_text(aes(label = Count), vjust = -0.5, size = 8) +
	  scale_fill_manual(values = treatment_colors) +
	  labs(title = "Macrophages as % of Total Immune Cells",
		   x = "Treatment",
		   y = "Percentage of Immune Cells (%)") +
	  theme_classic() +
	  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
			axis.title.x = element_text(size = 16),
			axis.title.y = element_text(size = 16),
			axis.text.x = element_text(size = 14),
			axis.text.y = element_text(size = 14),
			legend.title = element_text(size = 14),
			legend.text = element_text(size = 12),
			legend.position = "right") +
	  ylim(0, max(count_df_macro$Percentage) * 1.15)

	print(p_macrophage_percent)

	# Save plot
	ggsave(file.path(Figure6_dir, "macrophage_percentage_by_treatment.pdf"),
		   p_macrophage_percent, width = 8, height = 6)
	ggsave(file.path(Figure6_dir, "macrophage_percentage_by_treatment.png"),
		   p_macrophage_percent, width = 8, height = 6, dpi = 300)

	# Save the data table
	write.csv(count_df_macro, 
			  file.path(Figure6_dir, "macrophage_percentage_data.csv"), 
			  row.names = FALSE)

	print("=== MACROPHAGE PERCENTAGE BAR PLOT COMPLETE ===")
	print(paste("Macrophage percentage bar plot saved to:", 
				file.path(Figure6_dir, "macrophage_percentage_by_treatment.png")))

	# ============================================================================
	# Supplemental Figure: MACROPHAGE MARKER GENES
	# ============================================================================

	print("=== MACROPHAGE MARKER GENE ANALYSIS ===")

	# Rename macrophage clusters in metadata
	print("Renaming macrophage subtypes in metadata...")
	SSM2sc$cell_type_secondary <- as.character(SSM2sc$cell_type_secondary)
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "Macrophage 1"] <- "LAM"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "Macrophage 2"] <- "Inflammatory"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "Macrophage 3"] <- "AA"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "Macrophage 4"] <- "Undefined"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "Macrophage 5"] <- "TA"
	SSM2sc$cell_type_secondary[SSM2sc$cell_type_secondary == "Macrophage 6"] <- "Proliferative"

	# Save the updated Seurat object with renamed macrophages
	saveRDS(SSM2sc, file = file.path(base_dir, "SSM2sc_with_celltypes.RDS"))
	print("Updated Seurat object saved with renamed macrophage labels")

	# Define macrophage types (now using renamed labels)
	macrophage_types_renamed <- c("LAM", "Inflammatory", "AA", "Undefined", "TA", "Proliferative")

	# Subset to macrophages
	SSM2sc_macrophages_markers <- subset(SSM2sc, subset = cell_type_secondary %in% macrophage_types_renamed)

	print(paste("Filtered from", ncol(SSM2sc), "to", ncol(SSM2sc_macrophages_markers), "macrophage cells"))
	print("Cell counts by type:")
	print(table(SSM2sc_macrophages_markers$cell_type_secondary))

	# Filter mitochondrial and ribosomal genes
	print("Filtering genes...")
	rm.ind <- grep("^mt-|^Rp[ls]", rownames(SSM2sc_macrophages_markers))
	if(length(rm.ind) > 0) {
	  keep <- rownames(SSM2sc_macrophages_markers)[-rm.ind]
	  SSM2sc_macrophages_markers <- SSM2sc_macrophages_markers[keep,]
	  print(paste("Removed", length(rm.ind), "mitochondrial/ribosomal genes"))
	}

	# Filter lowly expressed genes (present in <20 cells)
	counts_data <- GetAssayData(SSM2sc_macrophages_markers, assay = "RNA", layer = "counts")
	keep.ind <- rowSums(counts_data > 0)
	keep.ind <- names(keep.ind[keep.ind >= 20])
	SSM2sc_macrophages_markers <- SSM2sc_macrophages_markers[keep.ind,]
	print(paste("Kept", length(keep.ind), "genes expressed in >=20 cells"))

	# Set up cell type order
	order_vec <- macrophage_types_renamed

	# Set factor levels for proper ordering
	SSM2sc_macrophages_markers$cell_type_secondary <- factor(SSM2sc_macrophages_markers$cell_type_secondary, 
															 levels = order_vec)
	Idents(SSM2sc_macrophages_markers) <- "cell_type_secondary"

	# Find marker genes for each macrophage type
	print("Finding marker genes for each macrophage type...")
	main_list <- data.frame()

	for (cell_type in order_vec) {
	  print(paste("Processing:", cell_type))
	  
	  # Check if this cell type exists in the data
	  if(sum(SSM2sc_macrophages_markers$cell_type_secondary == cell_type) > 10) {
		
		tryCatch({
		  # Find markers for this cell type vs all others
		  cluster_markers <- FindMarkers(SSM2sc_macrophages_markers, 
										 ident.1 = cell_type,
										 only.pos = TRUE,
										 min.pct = 0.25,
										 logfc.threshold = 0.25,
										 test.use = "wilcox",
										 verbose = FALSE)
		  
		  # Filter for positive fold changes
		  cluster_markers <- cluster_markers[cluster_markers$avg_log2FC > 0, ]
		  
		  if(nrow(cluster_markers) > 0) {
			# Apply FDR correction
			cluster_markers$p_val_fdr <- p.adjust(cluster_markers$p_val, method = "fdr")
			
			# Order by FDR-corrected p-values, then by fold change
			cluster_markers <- cluster_markers[order(cluster_markers$p_val_fdr, 
													 -cluster_markers$avg_log2FC), ]
			
			# Filter for FDR significance (FDR p < 0.05)
			significant_markers <- cluster_markers[cluster_markers$p_val_fdr < 0.05, ]
			
			if(nrow(significant_markers) > 0) {
			  # Take top 10 FDR-significant markers
			  n_markers <- min(10, nrow(significant_markers))
			  cluster_markers <- significant_markers[1:n_markers, ]
			  
			  cluster_markers$Gene <- rownames(cluster_markers)
			  cluster_markers$cell <- cell_type
			  
			  main_list <- rbind(main_list, cluster_markers)
			  print(paste("  Found", n_markers, "FDR-significant markers"))
			} else {
			  print(paste("  No FDR-significant markers found for", cell_type))
			}
		  }
		}, error = function(e) {
		  print(paste("  Error processing", cell_type, ":", e$message))
		})
	  } else {
		print(paste("  Skipping", cell_type, "- too few cells"))
	  }
	}

	# Create outputs if markers were found
	if(nrow(main_list) > 0) {
	  print(paste("Found", nrow(main_list), "total FDR-significant marker genes"))
	  
	  # Create dot plot
	  print("Creating macrophage marker gene dot plot...")
	  
	  p_macrophage_markers <- DotPlot(SSM2sc_macrophages_markers, 
									  features = unique(main_list$Gene),
									  cols = c(low = "#f7f7f7", high = "#d7191c"),
									  dot.min = 0.1) +
		scale_fill_gradient2(low = "#2c7bb6",
							 mid = "#f7f7f7",
							 high = "#d7191c") +
		coord_flip() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
			  axis.text.y = element_text(size = 10),
			  axis.title = element_blank(),
			  plot.title = element_text(hjust = 0.5)) +
		ggtitle("Macrophage Subtype Marker Genes")
	  
	  # Display plot
	  print(p_macrophage_markers)
	  
	  # Save outputs
	  ggsave(file.path(supplemental_dir, "macrophage_marker_dotplot.pdf"), 
			 p_macrophage_markers, width = 12, height = 14)
	  ggsave(file.path(supplemental_dir, "macrophage_marker_dotplot.png"), 
			 p_macrophage_markers, width = 12, height = 14, dpi = 300)
	  write.csv(main_list, 
				file.path(supplemental_dir, "macrophage_top_markers.csv"), 
				row.names = FALSE)
	  
	  print("=== MACROPHAGE MARKER ANALYSIS COMPLETE ===")
	  print(paste("Marker dot plot saved to:", file.path(Figure6_dir, "macrophage_marker_dotplot.png")))
	  print(paste("Marker gene list saved to:", file.path(Figure6_dir, "macrophage_top_markers.csv")))
	  
	} else {
	  print("No FDR-significant marker genes found for macrophages")
	}

	# ============================================================================
	# Figure 6: MACROPHAGE UMAP
	# ============================================================================

	print("=== CREATING MACROPHAGE UMAP ===")

	# Subset to macrophages for visualization (using the full SSM2sc object to preserve UMAP)
	SSM2sc_macrophages_umap <- subset(SSM2sc, subset = cell_type_secondary %in% macrophage_types_renamed)

	# Set factor levels to match the order of colors
	SSM2sc_macrophages_umap$cell_type_secondary <- factor(
	  SSM2sc_macrophages_umap$cell_type_secondary,
	  levels = c("LAM", "Inflammatory", "AA", "Undefined", "TA", "Proliferative")
	)

	# Define custom colors for macrophages
	macrophage_colors <- c(
	  "LAM" = "#16FF32",
	  "Inflammatory" = "#3283FE",
	  "AA" = "#7ED7D1",
	  "Undefined" = "#325A9B",
	  "TA" = "#DEA0FD",
	  "Proliferative" = "#C4451C"
	)

	# Create UMAP plot using the existing UMAP coordinates
	p_macrophage_umap <- DimPlot(SSM2sc_macrophages_umap, 
								 reduction = "umap", 
								 group.by = "cell_type_secondary",
								 cols = macrophage_colors,
								 pt.size = 0.5) +
	  labs(title = "Macrophage Subtypes", x = "UMAP 1", y = "UMAP 2") +
	  theme_classic() +
	  theme(
		plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
		legend.title = element_blank(),
		legend.text = element_text(size = 12),
		legend.position = "right",
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12)
	  ) +
	  guides(color = guide_legend(override.aes = list(size = 4)))

	print(p_macrophage_umap)

	# Save UMAP
	ggsave(file.path(Figure6_dir, "macrophage_umap.pdf"), 
		   p_macrophage_umap, width = 10, height = 8)
	ggsave(file.path(Figure6_dir, "macrophage_umap.png"), 
		   p_macrophage_umap, width = 10, height = 8, dpi = 300)

	print("=== MACROPHAGE UMAP COMPLETE ===")
	print(paste("Macrophage UMAP saved to:", file.path(Figure6_dir, "macrophage_umap.png")))


	# ============================================================================
	# Figure 6: MACROPHAGE UMAP BY TREATMENT
	# ============================================================================

	print("=== CREATING TREATMENT-SPECIFIC MACROPHAGE UMAPS ===")

	# Create UMAP for placebo treatment only
	print("Creating placebo macrophage UMAP...")
	SSM2sc_macrophages_placebo <- subset(SSM2sc_macrophages_umap, subset = treatment == "placebo")

	p_macrophage_umap_placebo <- DimPlot(SSM2sc_macrophages_placebo, 
										 reduction = "umap", 
										 group.by = "cell_type_secondary",
										 cols = macrophage_colors,
										 pt.size = 0.5) +
	  labs(title = "Placebo", x = "UMAP 1", y = "UMAP 2") +
	  theme_classic() +
	  theme(
		plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
		legend.title = element_blank(),
		legend.text = element_text(size = 12),
		legend.position = "right",
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12)
	  ) +
	  guides(color = guide_legend(override.aes = list(size = 4)))

	print(p_macrophage_umap_placebo)
	ggsave(file.path(Figure6_dir, "macrophage_umap_placebo.pdf"), 
		   p_macrophage_umap_placebo, width = 10, height = 8)
	ggsave(file.path(Figure6_dir, "macrophage_umap_placebo.png"), 
		   p_macrophage_umap_placebo, width = 10, height = 8, dpi = 300)

	# Create UMAP for progesterone treatment only
	print("Creating progesterone macrophage UMAP...")
	SSM2sc_macrophages_progesterone <- subset(SSM2sc_macrophages_umap, subset = treatment == "progesterone")

	p_macrophage_umap_progesterone <- DimPlot(SSM2sc_macrophages_progesterone, 
											  reduction = "umap", 
											  group.by = "cell_type_secondary",
											  cols = macrophage_colors,
											  pt.size = 0.5) +
	  labs(title = "Progesterone", x = "UMAP 1", y = "UMAP 2") +
	  theme_classic() +
	  theme(
		plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
		legend.title = element_blank(),
		legend.text = element_text(size = 12),
		legend.position = "right",
		axis.title = element_text(size = 14),
		axis.text = element_text(size = 12)
	  ) +
	  guides(color = guide_legend(override.aes = list(size = 4)))

	print(p_macrophage_umap_progesterone)
	ggsave(file.path(Figure6_dir, "macrophage_umap_progesterone.pdf"), 
		   p_macrophage_umap_progesterone, width = 10, height = 8)
	ggsave(file.path(Figure6_dir, "macrophage_umap_progesterone.png"), 
		   p_macrophage_umap_progesterone, width = 10, height = 8, dpi = 300)

	print("=== TREATMENT-SPECIFIC MACROPHAGE UMAPS COMPLETE ===")


	# ============================================================================
	# Figure 6: MACROPHAGE STACKED BAR PLOT
	# ============================================================================

	print("=== CREATING MACROPHAGE STACKED BAR PLOT ===")

	# Extract metadata and filter for macrophages
	macrophage_subset_data <- SSM2sc@meta.data[SSM2sc@meta.data$cell_type_secondary %in% macrophage_types_renamed, ]

	# Create summary table with proportions by treatment
	macrophage_cluster_counts <- macrophage_subset_data %>%
	  group_by(treatment, cell_type_secondary) %>%
	  summarise(count = n(), .groups = 'drop') %>%
	  group_by(treatment) %>%
	  mutate(proportion = count / sum(count))

	# Create stacked bar plot
	p_macrophage_proportion <- ggplot(macrophage_cluster_counts, 
									  aes(x = treatment, y = proportion, fill = cell_type_secondary)) +
	  geom_bar(stat = "identity", position = "stack", width = 0.6) +
	  scale_fill_manual(values = macrophage_colors) +
	  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00), 
						 labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
	  labs(x = "", y = "Proportion", fill = "Macrophage Subtype") +
	  theme_minimal() +
	  theme(
		panel.grid.major.x = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_rect(fill = "white", color = NA),
		plot.background = element_rect(fill = "white", color = NA),
		axis.text.x = element_text(size = 14, color = "black"),
		axis.text.y = element_text(size = 12, color = "black"),
		axis.title.y = element_text(size = 14),
		axis.ticks = element_line(color = "black"),
		legend.position = "right",
		legend.text = element_text(size = 11),
		legend.title = element_text(size = 12),
		plot.margin = margin(t = 10, r = 20, b = 10, l = 50)
	  ) +
	  coord_flip()

	#print(p_macrophage_proportion)

	# Save stacked bar plot
	ggsave(file.path(Figure6_dir, "macrophage_proportion_stacked.pdf"), 
		   p_macrophage_proportion, width = 10, height = 4)
	ggsave(file.path(Figure6_dir, "macrophage_proportion_stacked.png"), 
		   p_macrophage_proportion, width = 10, height = 4, dpi = 300)

	print("=== MACROPHAGE STACKED BAR PLOT COMPLETE ===")
	print(paste("Macrophage stacked bar plot saved to:", file.path(Figure6_dir, "macrophage_proportion_stacked.png")))

	# ============================================================================
	# Figure 6: MACROPHAGE STATISTICAL COMPARISON
	# ============================================================================

	print("=== PERFORMING MACROPHAGE STATISTICAL ANALYSIS ===")

	# Get total cells per treatment
	treatment_totals <- table(macrophage_subset_data$treatment)
	treatments <- names(treatment_totals)

	cat("\nTotal macrophages per treatment:\n")
	print(treatment_totals)

	if(length(treatments) == 2) {
	  
	  # Initialize results dataframe
	  results <- data.frame(
		cell_type = character(),
		treatment1_count = numeric(),
		treatment2_count = numeric(),
		treatment1_prop = numeric(),
		treatment2_prop = numeric(),
		fold_change = numeric(),
		p_value = numeric(),
		stringsAsFactors = FALSE
	  )
	  
	  # Name the columns based on actual treatments
	  names(results)[2:3] <- paste0(treatments, "_count")
	  names(results)[4:5] <- paste0(treatments, "_prop")
	  
	  # Perform Fisher's exact test for each macrophage type
	  for(cell_type in macrophage_types_renamed) {
		
		# Get counts for this cell type
		treatment1_this_type <- sum(macrophage_subset_data$treatment == treatments[1] & 
									  macrophage_subset_data$cell_type_secondary == cell_type)
		treatment2_this_type <- sum(macrophage_subset_data$treatment == treatments[2] & 
									  macrophage_subset_data$cell_type_secondary == cell_type)
		
		# Get counts for all other macrophage types combined
		treatment1_other_types <- treatment_totals[treatments[1]] - treatment1_this_type
		treatment2_other_types <- treatment_totals[treatments[2]] - treatment2_this_type
		
		# Create contingency table
		contingency_table <- matrix(c(treatment1_this_type, treatment1_other_types,
									  treatment2_this_type, treatment2_other_types),
									nrow = 2, byrow = TRUE)
		
		# Perform Fisher's exact test
		fisher_test <- fisher.test(contingency_table)
		
		# Calculate proportions
		treatment1_prop <- treatment1_this_type / treatment_totals[treatments[1]]
		treatment2_prop <- treatment2_this_type / treatment_totals[treatments[2]]
		
		# Calculate fold change (treatment2/treatment1)
		fold_change <- ifelse(treatment1_prop > 0, treatment2_prop / treatment1_prop, NA)
		
		# Add to results
		new_row <- data.frame(
		  cell_type = cell_type,
		  treatment1_count = treatment1_this_type,
		  treatment2_count = treatment2_this_type,
		  treatment1_prop = treatment1_prop,
		  treatment2_prop = treatment2_prop,
		  fold_change = fold_change,
		  p_value = fisher_test$p.value
		)
		names(new_row)[2:3] <- paste0(treatments, "_count")
		names(new_row)[4:5] <- paste0(treatments, "_prop")
		
		results <- rbind(results, new_row)
	  }
	  
	  # Apply FDR correction
	  results$p_adj <- p.adjust(results$p_value, method = "fdr")
	  
	  # Add significance indicators
	  results$significance <- ifelse(results$p_adj < 0.001, "***",
									 ifelse(results$p_adj < 0.01, "**",
											ifelse(results$p_adj < 0.05, "*", "ns")))
	  
	  # Sort by adjusted p-value
	  results <- results[order(results$p_adj), ]
	  
	  # Display results
	  cat(paste("\n=== MACROPHAGE STATISTICAL RESULTS ===\n"))
	  cat("Fisher's Exact Test with FDR correction\n")
	  cat(paste("Comparing", treatments[2], "vs", treatments[1], "for each macrophage subset\n\n"))
	  
	  print(results)
	  
	  # Save results as CSV
	  write.csv(results, file.path(Figure6_dir, "macrophage_statistical_comparison.csv"), row.names = FALSE)
	  
	  # Save results as Excel file (requires writexl or openxlsx package)
	  if(require(writexl, quietly = TRUE)) {
		write_xlsx(results, file.path(Figure6_dir, "macrophage_statistical_comparison.xlsx"))
		print("Statistical results saved as Excel file")
	  } else if(require(openxlsx, quietly = TRUE)) {
		write.xlsx(results, file.path(Figure6_dir, "macrophage_statistical_comparison.xlsx"))
		print("Statistical results saved as Excel file")
	  }
	  
	  # Identify significant differences
	  significant_results <- results[results$p_adj < 0.05, ]
	  
	  if(nrow(significant_results) > 0) {
		cat("\n=== SIGNIFICANT DIFFERENCES (FDR p < 0.05) ===\n")
		for(i in 1:nrow(significant_results)) {
		  row <- significant_results[i, ]
		  direction <- ifelse(row$fold_change > 1, "INCREASED", "DECREASED")
		  cat(sprintf("%s: %s in %s vs %s (FC=%.3f, FDR p=%.2e)\n", 
					  row$cell_type, direction, treatments[2], treatments[1], row$fold_change, row$p_adj))
		}
	  } else {
		cat("\n=== NO SIGNIFICANT DIFFERENCES (FDR p < 0.05) ===\n")
	  }
	  
	} else {
	  cat("Error: Expected 2 treatment groups, found", length(treatments), "\n")
	}

	print("=== MACROPHAGE STATISTICAL ANALYSIS COMPLETE ===")

	}

Figure7 <- function(){
	print("Generating Figure7")
	getData("SSM2sc_with_celltypes.RDS","SSM2sc")



	# ============================================================================
	# Figure 7: LAM GENE SIGNATURE ANALYSIS
	# ============================================================================

	# Create Figure 7 subdirectory
	Figure7_dir <- file.path(base_dir, "Figures", "Figure_7")
	if (!dir.exists(Figure7_dir)) dir.create(Figure7_dir, recursive = TRUE)

	# Verify the directory was created successfully
	if (!dir.exists(Figure7_dir)) {
	  stop(paste("ERROR: Could not create directory:", Figure7_dir,
				 "\nCheck that you have write permissions for this location."))
	}

	# Define macrophage types (now using renamed labels)
	  macrophage_types_renamed <- c("LAM", "Inflammatory", "AA", "Undefined", "TA", "Proliferative")

	  # Subset to macrophages for visualization (using the full SSM2sc object to preserve UMAP)
	  SSM2sc_macrophages_umap <- subset(SSM2sc, subset = cell_type_secondary %in% macrophage_types_renamed)

	  # Set factor levels to match the order of colors
	  SSM2sc_macrophages_umap$cell_type_secondary <- factor(
		SSM2sc_macrophages_umap$cell_type_secondary,
		levels = c("LAM", "Inflammatory", "AA", "Undefined", "TA", "Proliferative")
	  )

	print("=== Figure 7: LAM GENE SIGNATURE ANALYSIS ===")

	# Define LAM signature genes
	LAM_genes <- c("Trem2", "Lipa", "Ctsb", "Ctss", "Fabp5", "Lgals1", "Lgals3", "Cd9", "Apoe", "Spp1")

	# Check which genes are available
	available_LAM_genes <- intersect(LAM_genes, rownames(SSM2sc_macrophages_umap))
	missing_LAM_genes <- setdiff(LAM_genes, rownames(SSM2sc_macrophages_umap))

	cat("Available LAM genes:", paste(available_LAM_genes, collapse = ", "), "\n")
	if(length(missing_LAM_genes) > 0) {
	  cat("Missing genes:", paste(missing_LAM_genes, collapse = ", "), "\n")
	}

	# Calculate LAM signature score (needed for violin plot)
	print("Calculating LAM signature score...")
	SSM2sc_macrophages_umap <- AddModuleScore(
	  object = SSM2sc_macrophages_umap,
	  features = list(available_LAM_genes),
	  name = "LAM_Score",
	  ctrl = 50
	)

	# ============================================================================
	# Figure 7: INDIVIDUAL LAM GENE FEATURE PLOTS
	# ============================================================================

	print("Creating individual LAM gene feature plots...")

	# Define purple color palette
	purple_colors <- c("grey90", "#E6E6FA", "#DDA0DD", "#BA55D3", "#9932CC", "#8B008B", "#4B0082")

	# Create individual feature plots for each available LAM gene
	gene_plots <- list()

	for(gene in available_LAM_genes) {
	  cat("Creating feature plot for", gene, "...\n")
	  
	  p <- FeaturePlot(SSM2sc_macrophages_umap, 
					   features = gene,
					   reduction = "umap",
					   pt.size = 0.8,
					   cols = purple_colors,
					   order = TRUE) +
		ggtitle(gene) +
		labs(x = "UMAP 1", y = "UMAP 2") +
		theme_classic() +
		theme(
		  plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
		  legend.position = "right",
		  legend.title = element_blank(),
		  legend.text = element_text(size = 10),
		  axis.title = element_text(size = 12),
		  axis.text = element_text(size = 10),
		  plot.margin = margin(5, 5, 5, 5),
		  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
		)
	  
	  gene_plots[[gene]] <- p
	}

	# Arrange gene plots in grid (4 columns)
	if(length(gene_plots) > 0) {
	  n_plots <- length(gene_plots)
	  ncols <- 4
	  nrows <- ceiling(n_plots / ncols)
	  
	  gene_grid <- wrap_plots(gene_plots, ncol = ncols, nrow = nrows)
	  
	  # Save individual gene plots
	  ggsave(file.path(Figure7_dir, "LAM_genes_feature_plots.pdf"),
			 gene_grid, width = 16, height = 4 * nrows, dpi = 300)
	  ggsave(file.path(Figure7_dir, "LAM_genes_feature_plots.png"),
			 gene_grid, width = 16, height = 4 * nrows, dpi = 300)
	  
	  cat("Individual LAM gene feature plots saved!\n")
	}

	print("=== INDIVIDUAL LAM GENE FEATURE PLOTS COMPLETE ===")

	# ============================================================================
	# Figure 7: LAM SCORE BY TREATMENT
	# ============================================================================

	print("Creating LAM score violin plot by treatment...")

	# Define treatment colors (matching previous plots)
	treatment_colors <- c("placebo" = "#4292C6", "progesterone" = "#EF3B2C")

	# Prepare data for violin plot
	plot_data <- data.frame(
	  treatment = SSM2sc_macrophages_umap@meta.data$treatment,
	  LAM_Score = SSM2sc_macrophages_umap@meta.data$LAM_Score1
	)

	# Calculate statistics
	placebo_scores <- plot_data$LAM_Score[plot_data$treatment == "placebo"]
	progesterone_scores <- plot_data$LAM_Score[plot_data$treatment == "progesterone"]

	cat("Placebo macrophages:", length(placebo_scores), "\n")
	cat("Progesterone macrophages:", length(progesterone_scores), "\n")
	cat("Placebo median LAM score:", round(median(placebo_scores), 3), "\n")
	cat("Progesterone median LAM score:", round(median(progesterone_scores), 3), "\n")

	# Statistical test (Wilcoxon rank-sum test)
	wilcox_result <- wilcox.test(progesterone_scores, placebo_scores)

	# Determine significance stars
	get_significance_stars <- function(p_value) {
	  if(p_value < 0.0001) return("****")
	  else if(p_value < 0.001) return("***")
	  else if(p_value < 0.01) return("**")
	  else if(p_value < 0.05) return("*")
	  else return("ns")
	}

	significance_stars <- get_significance_stars(wilcox_result$p.value)

	cat("Wilcoxon test p-value:", wilcox_result$p.value, "\n")
	cat("Significance:", significance_stars, "\n")

	# Create violin plot with boxplot overlay
	p_lam_violin <- ggplot(plot_data, aes(x = treatment, y = LAM_Score, fill = treatment)) +
	  geom_violin(alpha = 0.8, color = "black", linewidth = 0.5, trim = FALSE) +
	  geom_boxplot(width = 0.15, alpha = 0.9, color = "black", linewidth = 0.5,
				   outlier.shape = NA) +
	  scale_fill_manual(values = treatment_colors) +
	  scale_x_discrete(labels = c("placebo" = "Placebo", "progesterone" = "Progesterone")) +
	  labs(
		title = "LAM Score by Treatment",
		x = "Treatment",
		y = "LAM Signature Score"
	  ) +
	  theme_minimal() +
	  theme(
		plot.title = element_text(size = 16, hjust = 0.5, face = "bold", margin = margin(b = 20)),
		axis.text.x = element_text(size = 12, color = "black"),
		axis.text.y = element_text(size = 11, color = "black"),
		axis.title.x = element_text(size = 14, color = "black", margin = margin(t = 10)),
		axis.title.y = element_text(size = 14, color = "black", margin = margin(r = 10)),
		legend.position = "none",
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		plot.background = element_blank(),
		axis.line = element_line(color = "black", linewidth = 0.5),
		axis.ticks = element_line(color = "black", linewidth = 0.5),
		plot.margin = margin(20, 20, 20, 20)
	  )

	# Add significance bar and stars if significant
	if(wilcox_result$p.value < 0.05) {
	  # Calculate positions for significance bar
	  max_y <- max(plot_data$LAM_Score, na.rm = TRUE)
	  min_y <- min(plot_data$LAM_Score, na.rm = TRUE)
	  y_range <- max_y - min_y
	  bar_y <- max_y + 0.1 * y_range
	  star_y <- bar_y + 0.02 * y_range
	  
	  p_lam_violin <- p_lam_violin +
		# Add horizontal line
		annotate("segment", x = 1, xend = 2, y = bar_y, yend = bar_y, 
				 color = "black", linewidth = 0.5) +
		# Add vertical lines at ends
		annotate("segment", x = 1, xend = 1, y = bar_y, yend = bar_y - 0.02 * y_range, 
				 color = "black", linewidth = 0.5) +
		annotate("segment", x = 2, xend = 2, y = bar_y, yend = bar_y - 0.02 * y_range, 
				 color = "black", linewidth = 0.5) +
		# Add significance stars
		annotate("text", x = 1.5, y = star_y, label = significance_stars,
				 size = 5, fontface = "bold") +
		# Extend y-axis to accommodate the annotation
		scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
	} else {
	  # If not significant, just add "ns" annotation
	  max_y <- max(plot_data$LAM_Score, na.rm = TRUE)
	  min_y <- min(plot_data$LAM_Score, na.rm = TRUE)
	  y_range <- max_y - min_y
	  
	  p_lam_violin <- p_lam_violin +
		annotate("text", x = 1.5, y = max_y + 0.1 * y_range, 
				 label = "ns", size = 4, fontface = "bold") +
		scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
	}

	print(p_lam_violin)

	# Save violin plot
	ggsave(file.path(Figure7_dir, "LAM_score_by_treatment.pdf"), 
		   p_lam_violin, width = 6, height = 8)
	ggsave(file.path(Figure7_dir, "LAM_score_by_treatment.png"), 
		   p_lam_violin, width = 6, height = 8, dpi = 300)

	print("=== LAM SCORE VIOLIN PLOT COMPLETE ===")
	print(paste("All Figure 7 outputs saved to:", Figure7_dir))


}


cellChatAnalysis <- function(){
	print("Generating cellChat Data")

	# ============================================================================
	# CELLCHAT ANALYSIS: PLACEBO VS PROGESTERONE
	# Complete workflow from subsetting to analysis
	# ============================================================================

	rm(list = ls())

	# Load required libraries
	loadLibrary(CellChat)
	loadLibrary(patchwork)
	loadLibrary(Seurat)
	loadLibrary(ggplot2)
	loadLibrary(ComplexHeatmap)
	options(stringsAsFactors = FALSE)

	# ============================================================================
	# SET UP PATHS (MODIFY THESE AS NEEDED)
	# ============================================================================

	# Set base directory - modify this to your project location
	#base_dir <- "path/to/your/SSM2sc"  # e.g., "C:/Users/YourName/Desktop/SSM2sc"

	# Set up folder paths
	figures_dir <- file.path(base_dir, "Figures")
	output_dir <- file.path(figures_dir, "Figure_8")
	dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

	# Input file path
	input_file <- file.path(base_dir, "SSM2sc_with_celltypes.RDS")

	# Define color scheme
	colors_final <- c(
	  "Tumor" = "#F6222E",
	  "NK" = "#FEAF16", 
	  "CD8+ T cell" = "#2ED9FF",
	  "CD4+ T cell" = "#F8A19F",
	  "CD4+ Treg" = "#5A5156",
	  "LAM" = "#16FF32",
	  "Inflammatory" = "#3283FE",
	  "AA" = "#7ED7D1",
	  "Undefined" = "#325A9B",
	  "TA" = "#DEA0FD",
	  "Proliferative" = "#C4451C"
	)

	# ============================================================================
	# SUBSET AND RECATEGORIZE SEURAT OBJECT
	# ============================================================================

	cat("=== CREATING SUBSETTED SEURAT OBJECT ===\n")

	# Load the full Seurat object
	cat("Loading Seurat object from:", input_file, "\n")
	seurat_full <- readRDS(input_file)
	cat("Original object contains", ncol(seurat_full), "cells\n")

	# Check current cell types
	cat("\nOriginal cell_type_secondary distribution:\n")
	print(table(seurat_full@meta.data$cell_type_secondary))

	# Define which cell types to keep (before recategorization)
	cell_types_to_keep <- c(
	  "CD4+ CM", "CD4+ Naive",  # Will become CD4+ T cell
	  "CD8+ CM", "CD8+ TRM",     # Will become CD8+ T cell
	  "Tumor 1", "Tumor 2", "Tumor 3", "Tumor 4", "Tumor 5", 
	  "Tumor 6", "Tumor 7", "Tumor 8", "Tumor 9", "Tumor 10", "Tumor 11",  # Will become Tumor
	  "CD4+ Treg", "Inflammatory", "Undefined", "AA", "TA", 
	  "Proliferative", "NK", "LAM", "Fibroblast 1", "Fibroblast 2"  # Keep as is
	)

	# Subset the Seurat object
	cat("\nSubsetting Seurat object...\n")
	SSM2sc_subset <- subset(seurat_full, 
							subset = cell_type_secondary %in% cell_types_to_keep)
	cat("After subsetting:", ncol(SSM2sc_subset), "cells\n")

	# Create new recategorized cell type column
	cat("\nRecategorizing cell types...\n")
	SSM2sc_subset$cell_type_combined <- SSM2sc_subset$cell_type_secondary

	# Combine CD4+ subtypes into CD4+ T cell
	SSM2sc_subset$cell_type_combined[SSM2sc_subset$cell_type_secondary %in% c("CD4+ CM", "CD4+ Naive")] <- "CD4+ T cell"
	# Combine Fibroblast into one 
	SSM2sc_subset$cell_type_combined[SSM2sc_subset$cell_type_secondary %in% c("Fibroblast 1", "Fibroblast 2")] <- "Fibroblasts"
	
	# Combine CD8+ subtypes into CD8+ T cell
	SSM2sc_subset$cell_type_combined[SSM2sc_subset$cell_type_secondary %in% c("CD8+ CM", "CD8+ TRM")] <- "CD8+ T cell"

	# Combine all Tumor subtypes into Tumor
	tumor_types <- paste0("Tumor ", 1:11)
	SSM2sc_subset$cell_type_combined[SSM2sc_subset$cell_type_secondary %in% tumor_types] <- "Tumor"

	# Show the new distribution
	cat("\nNew cell_type_combined distribution:\n")
	print(table(SSM2sc_subset$cell_type_combined))

	# Save the subset to Figure_8 folder
	subset_output <- file.path(output_dir, "SSM2sc_subset.RDS")
	cat("\nSaving subset to:", subset_output, "\n")
	saveRDS(SSM2sc_subset, subset_output)
	cat("SSM2sc_subset saved\n")

	# ============================================================================
	# CREATE PLACEBO CELLCHAT OBJECT
	# ============================================================================

	cat("\n=== CREATING PLACEBO CELLCHAT OBJECT ===\n")

	cell.use.placebo <- rownames(SSM2sc_subset@meta.data)[SSM2sc_subset@meta.data$treatment == "placebo"]
	data.input.placebo <- GetAssayData(SSM2sc_subset, layer = "data")[, cell.use.placebo]
	meta.placebo <- SSM2sc_subset@meta.data[cell.use.placebo, ]
	meta.placebo$labels <- factor(meta.placebo$cell_type_combined)

	cellchat_placebo <- createCellChat(object = data.input.placebo, meta = meta.placebo, group.by = "labels")
	cellchat_placebo <- addMeta(cellchat_placebo, meta = meta.placebo)
	cellchat_placebo <- setIdent(cellchat_placebo, ident.use = "labels")

	cellchat_placebo@DB <- CellChatDB.mouse
	cellchat_placebo <- subsetData(cellchat_placebo)
	cellchat_placebo <- identifyOverExpressedGenes(cellchat_placebo)
	cellchat_placebo <- identifyOverExpressedInteractions(cellchat_placebo)
	cellchat_placebo <- computeCommunProb(cellchat_placebo)
	cellchat_placebo <- filterCommunication(cellchat_placebo, min.cells = 10)
	cellchat_placebo <- computeCommunProbPathway(cellchat_placebo)
	cellchat_placebo <- aggregateNet(cellchat_placebo)
	cellchat_placebo <- netAnalysis_computeCentrality(cellchat_placebo, slot.name = "netP")

	cat("Placebo CellChat object created\n")

	# ============================================================================
	# CREATE PROGESTERONE CELLCHAT OBJECT
	# ============================================================================

	cat("\n=== CREATING PROGESTERONE CELLCHAT OBJECT ===\n")

	cell.use.progesterone <- rownames(SSM2sc_subset@meta.data)[SSM2sc_subset@meta.data$treatment == "progesterone"]
	data.input.progesterone <- GetAssayData(SSM2sc_subset, layer = "data")[, cell.use.progesterone]
	meta.progesterone <- SSM2sc_subset@meta.data[cell.use.progesterone, ]
	meta.progesterone$labels <- factor(meta.progesterone$cell_type_combined)

	cellchat_progesterone <- createCellChat(object = data.input.progesterone, meta = meta.progesterone, group.by = "labels")
	cellchat_progesterone <- addMeta(cellchat_progesterone, meta = meta.progesterone)
	cellchat_progesterone <- setIdent(cellchat_progesterone, ident.use = "labels")

	cellchat_progesterone@DB <- CellChatDB.mouse
	cellchat_progesterone <- subsetData(cellchat_progesterone)
	cellchat_progesterone <- identifyOverExpressedGenes(cellchat_progesterone)
	cellchat_progesterone <- identifyOverExpressedInteractions(cellchat_progesterone)
	cellchat_progesterone <- computeCommunProb(cellchat_progesterone)
	cellchat_progesterone <- filterCommunication(cellchat_progesterone, min.cells = 10)
	cellchat_progesterone <- computeCommunProbPathway(cellchat_progesterone)
	cellchat_progesterone <- aggregateNet(cellchat_progesterone)
	cellchat_progesterone <- netAnalysis_computeCentrality(cellchat_progesterone, slot.name = "netP")

	cat("Progesterone CellChat object created\n")

	# ============================================================================
	# MERGE AND COMPARE CELLCHAT OBJECTS
	# ============================================================================

	cat("\n=== MERGING CELLCHAT OBJECTS ===\n")

	object.list <- list(Placebo = cellchat_placebo, Progesterone = cellchat_progesterone)
	cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

	future::plan("sequential")

	cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
	cellchat_merged <- netEmbedding(cellchat_merged, type = "functional", umap.method = "uwot")

	tryCatch({
	  cellchat_merged <- netClustering(cellchat_merged, type = "functional")
	}, error = function(e) {
	  cat("Skipping netClustering\n")
	})

	# ============================================================================
	# COMPARE SIGNIFICANT PATHWAYS
	# ============================================================================

	cat("\n=== COMPARING SIGNIFICANT PATHWAYS ===\n")

	gg1 <- rankNet(cellchat_merged, mode = "comparison",
				   stacked = TRUE, do.stat = TRUE, measure = "weight", font.size = 11)

	pdf(file.path(output_dir, "significant_pathways_placebo_vs_progesterone.pdf"),
		width = 5, height = 13)
	print(gg1)
	dev.off()

	write.csv(gg1$data, file.path(output_dir, "significant_pathways_data.csv"), row.names = FALSE)

	# ============================================================================
	# CREATE SIGNALING ROLE HEATMAPS
	# ============================================================================

	cat("\n=== CREATING SIGNALING ROLE HEATMAPS ===\n")

	# Top 30 pathways (15 from each end)
	gg1_data <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
	top_paths <- gg1_data$data
	top_paths <- c(tail(top_paths, 15)$name, head(top_paths, 15)$name)

	ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]],
											 pattern = "all",
											 signaling = top_paths,
											 title = "Placebo",
											 width = 8.5, height = 12,
											 color.heatmap = "GnBu",
											 cluster.rows = FALSE,
											 cluster.cols = FALSE,
											 color.use = colors_final,
											 font.size = 12)

	ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], 
											 pattern = "all", 
											 signaling = top_paths,
											 title = "Progesterone",
											 width = 8.5, height = 12,
											 color.heatmap = "GnBu",
											 cluster.rows = FALSE,
											 cluster.cols = FALSE,
											 color.use = colors_final,
											 font.size = 12)

	pdf(file.path(output_dir, "heatmap_selected_pathways.pdf"),
		width = 20, height = 26)
	ht_list <- ht1 + ht2
	draw(ht_list, ht_gap = unit(0.5, "cm"))
	dev.off()

	# All pathways heatmap
	pathway.union <- union(object.list[[1]]@netP$pathways, 
						   object.list[[2]]@netP$pathways)

	ht1_all <- netAnalysis_signalingRole_heatmap(object.list[[1]],
												 pattern = "all",
												 signaling = pathway.union,
												 title = "Placebo",
												 width = 8, height = 36,
												 color.heatmap = "GnBu",
												 cluster.rows = FALSE,
												 cluster.cols = FALSE,
												 color.use = colors_final,
												 font.size = 12)

	ht2_all <- netAnalysis_signalingRole_heatmap(object.list[[2]], 
												 pattern = "all", 
												 signaling = pathway.union,
												 title = "Progesterone",
												 width = 8, height = 36,
												 color.heatmap = "GnBu",
												 cluster.rows = FALSE,
												 cluster.cols = FALSE,
												 color.use = colors_final,
												 font.size = 12)

	pdf(file.path(output_dir, "heatmap_all_pathways.pdf"),
		width = 20, height = 36)
	ht_list_all <- ht1_all + ht2_all
	draw(ht_list_all, ht_gap = unit(0.5, "cm"))
	dev.off()

	# ============================================================================
	# CREATE NETWORK COMPARISON PLOTS
	# ============================================================================

	cat("\n=== CREATING NETWORK COMPARISON PLOTS ===\n")

	pathways.show <- as.character(unique(gg1$data$name))
	weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))

	# Create subdirectory for networks
	network_dir <- file.path(output_dir, "network_plots")
	dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)

	for(j in pathways.show) {
	  safe_name <- gsub("[^A-Za-z0-9]", "_", j)
	  
	  # Placebo network
	  pdf(file.path(network_dir, paste0("Placebo_", safe_name, ".pdf")),
		  width = 9, height = 8)
	  tryCatch({
		netVisual_aggregate(object.list[[1]], signaling = j,
							vertex.label.cex = 2.1,
							color.use = colors_final,
							arrow.size = 1.2,
							remove.isolate = FALSE)
	  }, error = function(e) {
		cat("Error plotting Placebo", j, "\n")
	  })
	  dev.off()
	  
	  # Progesterone network
	  pdf(file.path(network_dir, paste0("Progesterone_", safe_name, ".pdf")),
		  width = 9, height = 8)
	  tryCatch({
		netVisual_aggregate(object.list[[2]], signaling = j,
							vertex.label.cex = 2.1,
							color.use = colors_final,
							arrow.size = 1.2,
							remove.isolate = FALSE)
	  }, error = function(e) {
		cat("Error plotting Progesterone", j, "\n")
	  })
	  dev.off()
	  
	  cat("Saved plots for", j, "\n")
	}

	# ============================================================================
	# SAVE CELLCHAT OBJECTS
	# ============================================================================

	cat("\n=== SAVING CELLCHAT OBJECTS ===\n")

	saveRDS(cellchat_placebo, file.path(output_dir, "cellchat_placebo.rds"))
	saveRDS(cellchat_progesterone, file.path(output_dir, "cellchat_progesterone.rds"))
	saveRDS(cellchat_merged, file.path(output_dir, "cellchat_merged.rds"))

	# ============================================================================
	# SUMMARY
	# ============================================================================

	cat("\n=== CELLCHAT ANALYSIS COMPLETE ===\n")
	cat("All files saved to:", output_dir, "\n")
	cat("\nFiles generated:\n")
	cat("1. SSM2sc_subset.RDS (in Figure_8 folder)\n")
	cat("2. significant_pathways_placebo_vs_progesterone.pdf\n")
	cat("3. significant_pathways_data.csv\n")
	cat("4. heatmap_selected_pathways.pdf\n")
	cat("5. heatmap_all_pathways.pdf\n")
	cat("6. network_plots/ folder with individual pathway comparisons\n")
	cat("7. cellchat_placebo.rds\n")
	cat("8. cellchat_progesterone.rds\n")
	cat("9. cellchat_merged.rds\n")

	cat("\nCell type summary:\n")
	print(table(SSM2sc_subset$cell_type_combined))

}

figure8 <- function(){
	print("Generating figure8")
	output_dir <- file.path(figures_dir, "Figure_8")
	cellchatObjects=c("cellchat_placebo","cellchat_progesterone","cellchat_merged")
	for (obj in cellchatObjects){
		if (!(obj %in% ls(envir = .GlobalEnv))){
			if (!file.exists(paste0(output_dir,"/",obj,".rds"))){
					cellChatAnalysis()
	  		}
			assign(obj,readRDS(paste0(output_dir,"/",obj,".rds")))
	  	}
	}
	getData("SSM2sc_with_celltypes.RDS","SSM2sc")

	# ============================================================================
	# Figure 8: SPP1 EXPRESSION VIOLIN PLOT IN MACROPHAGES
	# ============================================================================

	loadLibrary(ggplot2)
	loadLibrary(dplyr)
	loadLibrary(Seurat)

	# Set Figure 8 directory path
	figure8_dir <- file.path(base_dir, "Figures", "Figure_8")

	print("=== CREATING SPP1 EXPRESSION VIOLIN PLOT ===")

	# Define macrophage types (using your renamed types)
	macrophage_types <- c("LAM", "AA", "TA", "Proliferative", "Inflammatory", "Undefined")

	# Subset to macrophages only
	SSM2sc_macrophages <- subset(SSM2sc, subset = cell_type_secondary %in% macrophage_types)

	print(paste("Total macrophages:", ncol(SSM2sc_macrophages)))

	# Extract SPP1 expression data
	spp1_data <- FetchData(SSM2sc_macrophages, 
						   vars = c("Spp1", "cell_type_secondary", "treatment"))

	# Rename columns for clarity
	colnames(spp1_data) <- c("Expression", "Cell_Type", "Treatment")

	# Set factor levels for cell types (to control order on x-axis)
	spp1_data$Cell_Type <- factor(spp1_data$Cell_Type, 
								  levels = c("Inflammatory", "LAM", "Proliferative", 
											 "TA", "AA", "Undefined"))

	# Capitalize treatment names
	spp1_data$Treatment <- factor(spp1_data$Treatment,
								  levels = c("placebo", "progesterone"),
								  labels = c("placebo", "progesterone"))

	# Perform Wilcoxon rank-sum test for each cell type
	print("\nPerforming Wilcoxon rank-sum tests (Mann-Whitney U)...")
	wilcox_results <- list()

	for(cell_type in levels(spp1_data$Cell_Type)) {
	  subset_data <- spp1_data[spp1_data$Cell_Type == cell_type, ]
	  
	  placebo_expr <- subset_data$Expression[subset_data$Treatment == "placebo"]
	  prog_expr <- subset_data$Expression[subset_data$Treatment == "progesterone"]
	  
	  if(length(placebo_expr) > 0 & length(prog_expr) > 0) {
		test_result <- wilcox.test(prog_expr, placebo_expr)
		wilcox_results[[cell_type]] <- test_result$p.value
		
		cat(cell_type, ": p-value =", format(test_result$p.value, scientific = TRUE), "\n")
	  }
	}

	# Define treatment colors
	treatment_colors <- c("placebo" = "#4292C6", "progesterone" = "#EF3B2C")

	# Create the violin plot with corrected boxplots
	p_spp1 <- ggplot(spp1_data, aes(x = Cell_Type, y = Expression, fill = Treatment)) +
	  geom_violin(scale = "width", trim = FALSE, position = position_dodge(0.9)) +
	  geom_boxplot(width = 0.15, position = position_dodge(0.9), 
				   outlier.shape = NA, color = "black", alpha = 0.7) +
	  scale_fill_manual(values = treatment_colors) +
	  labs(title = "SPP1 Expression in Macrophages",
		   x = "Cell Type",
		   y = "Log-normalized Gene Expression") +
	  theme_classic() +
	  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
			axis.title.x = element_text(size = 16),
			axis.title.y = element_text(size = 16),
			axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
			axis.text.y = element_text(size = 14),
			legend.title = element_text(size = 14),
			legend.text = element_text(size = 12),
			legend.position = "right")

	# Add significance stars for p < 0.001
	y_max <- max(spp1_data$Expression) + 0.5

	for(i in 1:length(wilcox_results)) {
	  if(wilcox_results[[i]] < 0.001) {
		cell_type_pos <- i
		p_spp1 <- p_spp1 + 
		  annotate("text", x = cell_type_pos, y = y_max + 0.5, 
				   label = "***", size = 8)
	  }
	}

	print(p_spp1)

	# Save plot to Figure_8 directory
	ggsave(file.path(figure8_dir, "SPP1_expression_macrophages.pdf"),
		   p_spp1, width = 12, height = 8)
	ggsave(file.path(figure8_dir, "SPP1_expression_macrophages.png"),
		   p_spp1, width = 12, height = 8, dpi = 300)

	# Save statistical results
	wilcox_df <- data.frame(
	  Cell_Type = names(wilcox_results),
	  P_value = unlist(wilcox_results),
	  Significance = ifelse(unlist(wilcox_results) < 0.001, "***",
							ifelse(unlist(wilcox_results) < 0.01, "**",
								   ifelse(unlist(wilcox_results) < 0.05, "*", "ns")))
	)

	write.csv(wilcox_df, 
			  file.path(figure8_dir, "SPP1_expression_statistics.csv"), 
			  row.names = FALSE)

	print("=== Figure 8 COMPLETE ===")
	print(paste("SPP1 plot saved to:", file.path(figure8_dir, "SPP1_expression_macrophages.png")))
	print(paste("Statistics saved to:", file.path(figure8_dir, "SPP1_expression_statistics.csv")))

	cat("\nStatistical test used: Wilcoxon rank-sum test (Mann-Whitney U test)")
	cat("\nComparing progesterone vs placebo for each cell type independently\n")

	# ============================================================================
	# Figure 8: CD44 EXPRESSION VIOLIN PLOT IN T CELLS, NK CELLS, AND TUMOR CELLS
	# ============================================================================

	loadLibrary(ggplot2)
	loadLibrary(dplyr)
	loadLibrary(Seurat)

	# Set Figure 8 directory path
	figure8_dir <- file.path(base_dir, "Figures", "Figure_8")

	print("=== CREATING CD44 EXPRESSION VIOLIN PLOT ===")

	# Load the subset object
	SSM2sc_subset <- readRDS(file.path(figure8_dir, "SSM2sc_subset.RDS"))

	print(paste("Total cells in subset:", ncol(SSM2sc_subset)))

	# Extract CD44 expression data using cell_type_combined
	cd44_data <- FetchData(SSM2sc_subset, 
						   vars = c("Cd44", "cell_type_combined", "treatment"))

	# Rename columns for clarity
	colnames(cd44_data) <- c("Expression", "Cell_Type", "Treatment")

	# Define cell types to include (using the combined names)
	cell_types_to_plot <- c("CD4+ Treg", "CD8+ T cell", "NK")

	# Filter for only these cell types
	cd44_data <- cd44_data[cd44_data$Cell_Type %in% cell_types_to_plot, ]

	# Set factor levels for cell types (to control order on x-axis)
	cd44_data$Cell_Type <- factor(cd44_data$Cell_Type, 
								  levels = c( "CD4+ Treg", "CD8+ T cell", "NK"))

	# Treatment as factor
	cd44_data$Treatment <- factor(cd44_data$Treatment,
								  levels = c("placebo", "progesterone"),
								  labels = c("placebo", "progesterone"))

	print(paste("Total cells in plot:", nrow(cd44_data)))
	print("Cell type distribution:")
	print(table(cd44_data$Cell_Type))

	# Perform Wilcoxon rank-sum test for each cell type
	print("\nPerforming Wilcoxon rank-sum tests (Mann-Whitney U)...")
	wilcox_results <- list()

	for(cell_type in levels(cd44_data$Cell_Type)) {
	  subset_data <- cd44_data[cd44_data$Cell_Type == cell_type, ]
	  
	  placebo_expr <- subset_data$Expression[subset_data$Treatment == "placebo"]
	  prog_expr <- subset_data$Expression[subset_data$Treatment == "progesterone"]
	  
	  if(length(placebo_expr) > 0 & length(prog_expr) > 0) {
		test_result <- wilcox.test(prog_expr, placebo_expr)
		wilcox_results[[cell_type]] <- test_result$p.value
		
		cat(cell_type, ": p-value =", format(test_result$p.value, scientific = TRUE), "\n")
	  }
	}

	# Define treatment colors
	treatment_colors <- c("placebo" = "#4292C6", "progesterone" = "#EF3B2C")

	# Create the violin plot
	p_cd44 <- ggplot(cd44_data, aes(x = Cell_Type, y = Expression, fill = Treatment)) +
	  geom_violin(scale = "width", trim = FALSE, position = position_dodge(0.9)) +
	  geom_boxplot(width = 0.15, position = position_dodge(0.9), 
				   outlier.shape = NA, color = "black", alpha = 0.7) +
	  scale_fill_manual(values = treatment_colors) +
	  labs(title = "CD44 Expression in T cells, NK cells, and Tumor cells",
		   x = "Cell Type",
		   y = "Log-normalized Gene Expression") +
	  theme_classic() +
	  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
			axis.title.x = element_text(size = 16),
			axis.title.y = element_text(size = 16),
			axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
			axis.text.y = element_text(size = 14),
			legend.title = element_text(size = 14),
			legend.text = element_text(size = 12),
			legend.position = "right")

	# Add significance stars based on p-values
	y_max <- max(cd44_data$Expression) + 0.5

	for(i in 1:length(wilcox_results)) {
	  p_val <- wilcox_results[[i]]
	  cell_type_pos <- i
	  
	  if(p_val < 0.001) {
		# p < 0.001: ***
		p_cd44 <- p_cd44 + 
		  annotate("text", x = cell_type_pos, y = y_max + 0.5, 
				   label = "***", size = 8)
	  } else if(p_val < 0.01) {
		# p < 0.01: **
		p_cd44 <- p_cd44 + 
		  annotate("text", x = cell_type_pos, y = y_max + 0.5, 
				   label = "**", size = 8)
	  } else if(p_val < 0.05) {
		# p < 0.05: *
		p_cd44 <- p_cd44 + 
		  annotate("text", x = cell_type_pos, y = y_max + 0.5, 
				   label = "*", size = 8)
	  }
	}

	print(p_cd44)

	# Save plot to Figure_8 directory
	ggsave(file.path(figure8_dir, "CD44_expression_T_NK_Tumor.pdf"),
		   p_cd44, width = 12, height = 8)
	ggsave(file.path(figure8_dir, "CD44_expression_T_NK_Tumor.png"),
		   p_cd44, width = 12, height = 8, dpi = 300)

	# Save statistical results
	wilcox_df_cd44 <- data.frame(
	  Cell_Type = names(wilcox_results),
	  P_value = unlist(wilcox_results),
	  Significance = ifelse(unlist(wilcox_results) < 0.001, "***",
							ifelse(unlist(wilcox_results) < 0.01, "**",
								   ifelse(unlist(wilcox_results) < 0.05, "*", "ns")))
	)

	write.csv(wilcox_df_cd44, 
			  file.path(figure8_dir, "CD44_expression_statistics.csv"), 
			  row.names = FALSE)

	print("=== CD44 VIOLIN PLOT COMPLETE ===")
	print(paste("CD44 plot saved to:", file.path(figure8_dir, "CD44_expression_T_NK_Tumor.png")))
	print(paste("Statistics saved to:", file.path(figure8_dir, "CD44_expression_statistics.csv")))

	cat("\nStatistical test used: Wilcoxon rank-sum test (Mann-Whitney U test)")
	cat("\nComparing progesterone vs placebo for each cell type independently\n")

}


Make_All_Figures <- function(){
	print("I am making all the figures!")
	#setup()
	figure2()
	figure3()
	Figure4()
	Figure5()
	Figure6()
	Figure7()
	cellChatAnalysis()
	figure8()
}


helpMessage = function(){
	args <- commandArgs(FALSE)
	file_arg <- args[grepl("^--file=", args)]
	script_path <- basename(sub("^--file=", "", file_arg))
	message("
Usage:
Rscript ",script_path," [-h 2 3 4 5 6 7]


ARGUMENTS:
default			If no value is given, generates all figures images
-h/--help		Prints this help message
2-7			Generates that figures images. Multiple numbers can be entered.

EXAMPLE:
To generate images for figures 3 and 7

Rscript ",script_path," 3 7

")

}




if ( identical(parent.frame(), .GlobalEnv) && !interactive()) {
  functionList=c(helpMessage,figure2,figure3,Figure4,Figure5,Figure6,Figure7,figure8)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0){
  	print("Generating all figures")
  	Make_All_Figures()
  }else{
    for (x in args){
      functionList[[as.numeric(x)]]()
	}
  }
}


	# ============================================================================
	# END OF CODE
	# ============================================================================

