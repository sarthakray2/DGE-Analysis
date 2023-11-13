### Create ID vs Sample dataset
# Query GDC metadata without UUID filter
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification")

# Extract sample submitter IDs and UUIDs from the query result
sample_submitter_ids <- query[[1]][[1]]$sample.submitter_id
uuids <- query[[1]][[1]]$id

# Create a dataframe with the extracted information
sample_uuid_df <- data.frame(Sample_Submitter_ID = sample_submitter_ids, UUID = uuids)

write_xlsx(sample_uuid_df, "C:/Users/sarth/Downloads/idvssample.xlsx")

#################################################################################################
### Custom filter 
# Specify the UUID you want to query for
target_uuid <- "6b62d6dd-bb6f-4f39-9041-992c47b875f3"  # Replace with your UUID

# Filter the sample IDs based on the target UUID
target_sample_id <- sample_uuid_df$Sample_Submitter_ID[sample_uuid_df$UUID == target_uuid]

###########################################################################################################
### Read all tsv files and create the combined dataset
# Specify the main folder containing subfolders which contain TSV files
main_folder <- "D:/btp/TCGA/GDCdata/TCGA-BRCA/Transcriptome_Profiling/Gene_Expression_Quantification"

# List subfolders in the main folder
subfolders <- list.files(main_folder, full.names = TRUE)

# Initialize gene_ids and gene_names using the first TSV file
first_tsv_file <- list.files(subfolders[1], pattern = "*.tsv", full.names = TRUE)[1]
gene_ids <- read.table(first_tsv_file, header = TRUE, sep = "\t")$gene_id
gene_names <- read.table(first_tsv_file, header = TRUE, sep = "\t")$gene_name

# Remove irrelevant rows
gene_ids <- gene_ids[-(1:4)]
gene_names <- gene_names[-(1:4)]

# Set gene_ids as row names for the final_data dataframe
final_data <- data.frame(row.names = gene_ids)
final_data$GeneNames <- gene_names

# Loop through each subfolder
for (subfolder in subfolders) {
  # Extract the subfolder name from the path
  subfolder_name <- tools::file_path_sans_ext(basename(subfolder))
  
  # Find the corresponding sample ID from the subfolder name(id) in sample_uuid_df
  sample_id <- sample_uuid_df$Sample_Submitter_ID[sample_uuid_df$UUID == subfolder_name]
  
  # List the TSV file in the subfolder
  tsv_file <- list.files(subfolder, pattern = "*.tsv", full.names = TRUE)
  
  # Select the count data of the unstranded column from the tsv file
  unstranded_column <- read.table(tsv_file, header = TRUE, sep = "\t")$unstranded
  unstranded_column <- unstranded_column[-(1:4)]
  
  # Add this data under the corresponding sample ID
  final_data[sample_id] <- unstranded_column
}

write.xlsx(final_data, "C:/Users/sarth/Downloads/final_data.xlsx", row.names = TRUE)

##########################################################################################
# Clean the environment
remove(first_tsv_file)
remove(gene_ids)
remove(gene_names)
remove(main_folder)
remove(sample_id)
remove(sample_submitter_ids)
remove(subfolder)
remove(subfolder_name)
remove(subfolders)
remove(tsv_file)
remove(unstranded_column)
remove(uuids)
##################################################################################
### Obtain ER_status values from clinical data
# Import clinical data
clinicaldata <- TCGA.BRCA
remove(TCGA.BRCA)

# Extract sample IDs without the last character(representing vial)
sample_ids <- substr(sample_uuid_df$Sample_Submitter_ID, 1, 
                     nchar(sample_uuid_df$Sample_Submitter_ID) - 1)

# Find corresponding ER_status values from clinical data using sample_ids
er_status_values <- clinicaldata$ER_Status_nature2012[match(
  sample_ids, clinicaldata$sampleID)]

# Add the ER_status column to sample_uuid_df
sample_uuid_df$ER_status <- er_status_values

# Rename according to purpose
er_status <- sample_uuid_df
remove(sample_uuid_df)

write_xlsx(er_status, "C:/Users/sarth/Downloads/erstatus.xlsx")

############################################################################################################
### Filter and finetune data
# Find samples with empty ER status
samples_with_empty_er <- er_status$Sample_Submitter_ID[er_status$ER_status == "" 
                                                       | er_status$ER_status == "Indeterminate"]

# Remove columns from final_data where samples have empty or "Indeterminate" ER status
final_data_filtered <- final_data[, !(colnames(final_data) %in% samples_with_empty_er)]
final_data_filtered <- final_data_filtered[, -1]


# Reorder er_status based on the order of samples in final_data_filtered
er_status_filtered <- er_status[match(colnames(final_data_filtered),
                                      er_status$Sample_Submitter_ID), ]
er_status_filtered <- er_status_filtered[, -2]

# Convert ER_status to factor form
er_status_filtered$ER_status <- as.factor(er_status_filtered$ER_status)

write_xlsx(er_status_filtered, "C:/Users/sarth/Downloads/coldata.csv")
write_xlsx(final_data_filtered, "C:/Users/sarth/Downloads/countdata.csv")

###############################################################################################
### DESeq Analysis
# Create the dds object
dds <- DESeqDataSetFromMatrix(countData = final_data_filtered, colData = er_status_filtered, 
                              design = ~ ER_status)

# Filtering genes; different ways
# Calculate the 50th percentile of the sum of counts
threshold <- quantile(rowSums(counts(dds)), probs = 0.5)

# Filter genes based on the threshold
dds <- dds[rowSums(counts(dds)) > threshold, ]

# Normalization to account for different sequencing depth
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

contrast_results <- results(dds, contrast=c("ER_status", "Positive", "Negative"))
contrast_results_df$gene_name <- final_data[row.names(contrast_results_df), 1]

# MA Plot
plotMA(contrast_results)
##################################################################################################
### Filtering results based on custom criteria
# Apply filtering criteria
filtered_results <- contrast_results[abs(contrast_results$log2FoldChange) >= 1 
                                     & contrast_results$padj < 0.05, ]

# Add a new column "gene_name" to filtered_results with gene names from final_data
filtered_results$gene_name <- final_data[row.names(filtered_results), 1]

# View the filtered results
filtered_results

# Calculate the number of genes with a positive log2FoldChange
positive_log2fold_genes <- sum(filtered_results$log2FoldChange > 0)

# Subset the rows where gene names contain "MMP"
mmp_genes <- filtered_results[grepl("HOX", filtered_results$gene_name), ]

print(mmp_genes)

write_xlsx(as.data.frame(contrast_results_df), "C:/Users/sarth/Downloads/all_genes.csv")
###################################################################################
### Create a volcano plot to visualize differential expression results
# Convert DESeqResults to a data frame
contrast_results_df <- as.data.frame(contrast_results)
# Plot
volcano_plot <- ggplot(contrast_results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < 0.01 & abs(log2FoldChange) >= 3, "Significant", "Not Significant")))+
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-3, 3), linetype = "dashed", color = "red") +
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted p-value)") +
  theme_minimal() +
  ggtitle("Volcano Plot") +
  scale_color_manual(name = "Genes", values = c("Significant" = "blue", "Not Significant" = "black"))

print(volcano_plot)
###################################################################################
### Creating a heatmap for the significant genes
# Use variance stabilizing transformation
dds_vst <- vst(dds, blind = FALSE)
dds_vst <- assay(dds_vst)
dds_vst <- as.data.frame(dds_vst)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(filtered_results)
dds_vst <- dds_vst[rownames(dds_vst) %in% sigGenes,]

# Create a new data frame with gene names as the first column
genes_column <- data.frame(genes = filtered_results$gene_name)
# Add the remaining columns from dds_vst
dds_vst <- cbind(genes_column, dds_vst)

# Convert it to long format using melt
dds_vst_melted <- melt(dds_vst, id.vars = "genes") 

# Plot
heatmap_plot <- ggplot(dds_vst_melted, aes(x = variable, y = genes, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0) +
  labs(x = "Samples", y = "Genes", fill = "Expression") +  # Set axis labels
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove text on x-axis
        axis.text.y = element_blank()) + # Remove text on y-axis
 ggtitle("Heatmap of Significantly Differentiated Genes")

# Print the heatmap
print(heatmap_plot)


##################################################################################################
# Clean the environment
remove(er_status_values)
remove(sample_ids)
remove(samples_with_empty_er)
#############################################################