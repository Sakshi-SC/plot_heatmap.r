#####HEATMAP##########
getwd()
setwd("/home/rohan/Documents/Sakshi/plots/heatmap")
library(pheatmap)
#update.packages(ask = FALSE)
getwd()
# Extract normalized counts from the dds object
###Pattern Recognition: Heatmaps make it easy to spot patterns or relationships between variables (such as gene expression across different conditions).
###Identifying Clusters: Clustering can help identify groups of genes or samples that behave similarly under different conditions.
normalized_counts <- counts(dds1, normalized=TRUE)
head(normalized_counts)
#Step 3: filter out the DEGs genew from DESeq2 normalised result ;to include only these genes
#  DESeq2 results object, e.g., `res`
# Filtering DEGs based on adjusted p-value (e.g., < 0.05)
deg_list <- rownames(res[which(res$padj < 0.05), ])
head(deg_list)
# Subset the normalized counts matrix for DEGs
deg_matrix <- normalized_counts[deg_list, ]
head(deg_matrix)
#Step 4: Scaling the Data (Z-Score Scaling)
# Z-score scaling (standardization) across samples
scaled_deg_matrix <- t(scale(t(deg_matrix)))
head(scaled_deg_matrix)
# Create an annotation data frame
sample_info <- data.frame(
   Group = factor(c(rep("control", 3), rep("test", 3))),  # as i am having 3 control and 3 test samples
   row.names = colnames(scaled_deg_matrix)  # Making sure row names match column names in scaled_deg_matrix
)
head(sample_info)
# Define annotation colors
annotation_colors <- list(
   Group = c(
      control = "#708090",  
      test = "grey"      
   )
)


# Generate the heatmap with annotation colors
heatmap <- pheatmap(scaled_deg_matrix,
                    scale = "none",  # Data is already scaled, so no further scaling is applied
                    clustering_distance_rows = "euclidean",  
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",  d
                    show_rownames = FALSE,  # Display row names (genes)
                    show_colnames = TRUE,  # Display column names (samples)
                    color = colorRampPalette(c("lightblue", "skyblue", "deepskyblue", "dodgerblue", "blue"))(50),  
                    annotation_col = sample_info,  
                    annotation_colors = annotation_colors,  
                    fontsize_col = 10,  
                    annotation_names_col = FALSE,  # this is for Hiding annotation names on the horizontal bar
                    main = "")  # Title of the heatmap, keeping it none for now

# Save as PDF
pdf("heatmap_scaled_degs.pdf", width = 10, height = 8) 
pheatmap(scaled_deg_matrix,
         scale = "none",  # Data is already scaled, so no further scaling is applied
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean",  
         clustering_method = "complete",  
         show_rownames = FALSE,  
         show_colnames = TRUE,  
         color = colorRampPalette(c("lightblue", "skyblue", "deepskyblue", "dodgerblue", "blue"))(50),  
         annotation_col = sample_info,  
         annotation_colors = annotation_colors,  
         fontsize_col = 10,  
         annotation_names_col = FALSE,  
         main = "")  
dev.off()  # Close the PDF device
# Save as PNG
png("heatmap_scaled_degs.png", width = 1000, height = 800, res = 150)  
pheatmap(scaled_deg_matrix,
         scale = "none",  # Data is already scaled, so no further scaling is applied
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",  
         show_rownames = FALSE,  
         show_colnames = TRUE,  #
         color = colorRampPalette(c("lightblue", "skyblue", "deepskyblue", "dodgerblue", "blue"))(50), 
         annotation_col = sample_info, 
         annotation_colors = annotation_colors,  
         fontsize_col = 10, 
         annotation_names_col = FALSE, 
         main = "")  # Title of the heatmap
dev.off()  # Close the PNG device




# DE-ANALYSIS
# 
# library(DESeq2)
# setwd("/home/rohan/Documents/Sakshi/mh_rna_last/output/new_deseqanalysis")
# 
# # Loading the data in R
# countData <- read.table("unique_GC_matrix_for_new.csv", header=TRUE, sep=",")
# head(countData)
# rownames(countData) <- countData[,1] # Rename the row names with gene names
# countData <- countData[,-1] # Remove the gene column
# 
# # Loading the sample information
# colData <- read.csv("/home/rohan/Downloads/phenoinfo.csv", sep=",", row.names=1)
# head(colData)
# # Ensure the order of samples in colData matches the columns in countData
# countData <- countData[,rownames(colData)]
# 
# # Print column names of countData
# colnames(countData)
# 
# # Print row names of colData
# rownames(colData)
# 
# # Check if all row names in colData are in the columns of countData
# all(rownames(colData) %in% colnames(countData))
# 
# # Check if all columns in countData are in row names of colData
# all(colnames(countData) %in% rownames(colData))
# 
# 
# 
# # Check if row names of colData match the column names of countData
# if (!all(rownames(colData) %in% colnames(countData))) {
#    stop("Mismatch between sample names in colData and countData. Check the sample names in both files.")
# }
# 
# if (!all(rownames(colData) == colnames(countData))) {
#    stop("Order of samples in colData does not match the columns in countData. Check the order of samples.")
# }
# 
# # Convert 'condition' and 'batch' columns to factors
# colData$condition <- as.factor(colData$condition)
# colData$batch <- as.factor(colData$batch)
# 
# # Check and ensure colData contains 'batch'
# if (!"batch" %in% colnames(colData)) {
#    stop("The 'batch' column is missing in colData.")
# }
# 
# 
# # Filter out genes with a total count of less than 50 across all samples
# countData <- countData[rowSums(countData) >= 50, ]
# cat("Number of genes after filtering low counts:", nrow(countData), "\n")
# 
# 
# # Create DESeq2 dataset object with batch correction
# dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ batch + condition)
# 
# # Normalize counts
# dds1 <- estimateSizeFactors(dds)
# normalized_counts <- counts(dds1, normalized=TRUE)
# head(normalized_counts)
# 
# # Save normalized counts to a CSV file
# write.csv(normalized_counts, file="DEseq2_normalized_counts.csv", quote=F)
# 
# 
# # Run DESeq2 analysis
# dds <- DESeq(dds)
# 
# # Extract results for the condition of interest (Test vs. Control)
# res <- results(dds, contrast = c("condition", "Test", "Control"), alpha = 0.05)
# summary(res)
