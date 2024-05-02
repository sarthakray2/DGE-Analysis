##########BASICS################
hox_genes <- final_data_filtered[grepl("^HOX", final_data_filtered$GeneNames), ]

transposed_data <- t(hox_genes[, -1])  # Exclude the column with gene names

# Set the transposed gene names as column names
colnames(transposed_data) <- hox_genes$GeneNames

transposed_data <- data.frame(transposed_data)
transposed_data$status <- er_status_filtered$ER_status

# Reorder columns to place 'status' at the beginning
transposed_data <- transposed_data[, c('status', setdiff(names(transposed_data), 'status'))]

hox_data <- transposed_data

##########################################################################################################

mmp_genes <- final_data_filtered[grepl("^MMP", final_data_filtered$GeneNames), ]

transposed_data <- t(mmp_genes[, -1])  # Exclude the column with gene names

# Set the transposed gene names as column names
colnames(transposed_data) <- mmp_genes$GeneNames

transposed_data <- data.frame(transposed_data)
transposed_data$status <- er_status_filtered$ER_status

# Reorder columns to place 'status' at the beginning
transposed_data <- transposed_data[, c('status', setdiff(names(transposed_data), 'status'))]

mmp_data <- transposed_data

remove(transposed_data)
###############################################################################################

all_genes <- final_data_filtered[, 1]
all_gene_expression_data <- t(final_data_filtered)
colnames(all_gene_expression_data) <- all_genes
all_gene_expression_data <- all_gene_expression_data[-1,]
all_gene_expression_data <- data.frame(all_gene_expression_data)
all_gene_expression_data <- as.data.frame(lapply(all_gene_expression_data, function(x) as.numeric(as.character(x))))

#########################################################################################################

# Assuming gene_expression_data is your dataframe with ER status as the first column
# and gene expression values in the remaining columns

# Extract expression matrix (exclude the ER status column)
expression_matrix <- as.matrix(hox_data[, -1])

# Perform variance stabilizing transformation directly
vst_data <- varianceStabilizingTransformation(expression_matrix)

# If you want to add ER status column back to the vst_data
hox_vst_data <- cbind(hox_data[, 1, drop = FALSE], vst_data)

###########################################################################################################

# Extract expression matrix (exclude the ER status column)
expression_matrix <- as.matrix(mmp_data[, -1])

# Perform variance stabilizing transformation directly
vst_data <- varianceStabilizingTransformation(expression_matrix)

# If you want to add ER status column back to the vst_data
mmp_vst_data <- cbind(mmp_data[, 1, drop = FALSE], vst_data)

#############################################################################################################
#######################HOX genes###############################
# Load required libraries
library(ggplot2)

# Assuming 'hox_data' is your dataframe with the described structure

# Extracting status column and gene expression values
status_column <- hox_vst_data[, 1]
gene_expression_data <- hox_vst_data[, -1]  # Exclude the first column (status)

# Initialize a list to store plots
gene_plots <- list()

# Iterate over each gene and create separate plots
for (gene in colnames(gene_expression_data)) {
  # Create a dataframe for the current gene
  gene_data <- data.frame(Sample = rownames(hox_vst_data), Gene = rep(gene, nrow(hox_vst_data)), 
                          Expression = gene_expression_data[, gene], Status = status_column)
  
  # Plot using ggplot2
  gene_plot <- ggplot(gene_data, aes(x = Sample, y = Expression, color = Status)) +
    geom_point() +
    labs(title = paste("Expression of", gene),
         x = "Sample", y = "Expression") +
    theme_minimal()
  
  # Save the plot to the list
  gene_plots[[gene]] <- gene_plot
  
}

# Access individual plots from the list (e.g., the first plot)
print(gene_plots[[30]])

##################################################################################################################

# Initialize a list to store box plots
box_plots <- list()

# Iterate over each gene and create separate box plots
for (gene in colnames(gene_expression_data)) {
  # Create a dataframe for the current gene
  gene_data <- data.frame(Gene = rep(gene, nrow(gene_expression_data)),
                          Expression = gene_expression_data[, gene], 
                          Status = status_column)
  
  # Plot using ggplot2
  box_plot <- ggplot(gene_data, aes(x = Status, y = Expression, fill = Status)) +
    geom_boxplot() +
    labs(title = paste("Expression of", gene),
         x = "Status", y = "Expression") +
    theme_minimal()
  
  # Save the box plot to the list
  box_plots[[gene]] <- box_plot
  
}

# Access individual box plots from the list (e.g., the first plot)
print(box_plots[[30]])
##################################################################################################################
hoxb2_data <- data.frame(Gene = rep("HOXB2", nrow(gene_expression_data)),
                         Expression = gene_expression_data[, "HOXB2"], 
                         Status = status_column)

# Extract expression data for positive and negative status
positive_data <- hoxb2_data$Expression[hoxb2_data$Status == "Positive"]
negative_data <- hoxb2_data$Expression[hoxb2_data$Status == "Negative"]

# Perform Shapiro-Wilk test for normality
shapiro_positive <- shapiro.test(positive_data)
shapiro_negative <- shapiro.test(negative_data)

# Print the test results
print("Shapiro-Wilk Test for Normality - Positive Status:")
print(shapiro_positive)
print("Shapiro-Wilk Test for Normality - Negative Status:")
print(shapiro_negative)

# Create histograms
hist_positive <- ggplot(data.frame(Expression = positive_data), aes(x = Expression)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "HOXB2 - Positive Status",
       x = "Expression Level",
       y = "Frequency")

hist_negative <- ggplot(data.frame(Expression = negative_data), aes(x = Expression)) +
  geom_histogram(binwidth = 1, fill = "salmon", color = "black", alpha = 0.7) +
  labs(title = "HOXB2 - Negative Status",
       x = "Expression Level",
       y = "Frequency")

# Display histograms
print(hist_positive)
print(hist_negative)

# Perform Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(positive_data, negative_data)
print(wilcox_test_result)

####################################################################################################################

# Create an empty data frame to store correlation results
correlation_results <- data.frame(Gene = character(), Correlation = numeric(), stringsAsFactors = FALSE)

# Iterate over each gene and compute correlation with status
for (gene in colnames(gene_expression_data)) {
  # Create a dataframe for the current gene
  gene_data <- data.frame(Expression = gene_expression_data[, gene], Status = status_column)
  
  # Convert 'Status' to a binary variable
  gene_data$Status_binary <- ifelse(gene_data$Status == "Positive", 1, 0)
  
  # Compute correlation between status and gene expression
  correlation_value <- cor.test(gene_data$Expression, gene_data$Status_binary)$estimate
  
  # Append results to the data frame
  correlation_results <- rbind(correlation_results, data.frame(Gene = gene, Correlation = correlation_value))
}

correlation_results <- correlation_results[order(-(correlation_results$Correlation)), ]

# Display the correlation table
print(correlation_results)

################################################################################################################

# Create an empty data frame to store correlation results
correlation_results_2 <- data.frame(Gene = character(), HOXB2_Correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each gene and compute correlation with status
for (gene in colnames(all_gene_expression_data)) {
  # Create a dataframe for the current gene
  gene_data <- data.frame(Expression = all_gene_expression_data[, gene], HOXB2 = hox_data[, 31])
  
  # Compute correlation between status and gene expression
  correlation_value <- cor.test(gene_data$Expression, gene_data$HOXB2)$estimate
  p_value <- cor.test(gene_data$Expression, gene_data$HOXB2)$p.value
  
  # Append results to the data frame
  correlation_results_2 <- rbind(correlation_results_2, data.frame(Gene = gene, HOXB2_Correlation = correlation_value, p_value = p_value))
}

correlation_results_2 <- correlation_results_2[order((correlation_results_2$HOXB2_Correlation)), ]

# Display the correlation table
print(correlation_results_2)

#############################################################################################################
# Hierarchical clustering
dist_matrix <- dist(gene_expression_data, method = "euclidean")
hierarchical_clusters <- hclust(dist_matrix, method = "ward.D2")
row_dend <- as.dendrogram(hierarchical_clusters)

# Create a heatmap
heatmap_data <- data.matrix(gene_expression_data)
colnames(heatmap_data) <- NULL  # Remove gene names for simplicity

# Specify colors for the status column
status_colors <- ifelse(hox_data$status == "Positive", "red", "blue")

# Plot the heatmap with hierarchical clustering
pheatmap(
  heatmap_data,
  clustering_distance_rows = dist_matrix,
  clustering_distance_cols = dist(heatmap_data),
  clustering_method = "ward.D2",
  annotation_row = data.frame(Status = hox_data$status),
  annotation_colors = list(Status = status_colors),
  main = "Hierarchical Clustering of hox_data"
)
####################################################################################################################
#################################MMP genes##############################
# Load required libraries
library(ggplot2)

# Assuming 'mmp_data' is your dataframe with the described structure

# Extracting status column and gene expression values
status_column_b <- mmp_vst_data[, 1]
gene_expression_data_b <- mmp_vst_data[, -1]  # Exclude the first column (status)
status_column_b2 <- mmp_data[, 1]
gene_expression_data_b2 <- mmp_data[, -1]

# Initialize a list to store plots
gene_plots_b <- list()

# Iterate over each gene and create separate plots
for (gene in colnames(gene_expression_data_b)) {
  # Create a dataframe for the current gene
  gene_data_b <- data.frame(Sample = rownames(mmp_vst_data), Gene = rep(gene, nrow(mmp_vst_data)), 
                            Expression = gene_expression_data_b[, gene], Status = status_column_b)
  
  # Plot using ggplot2
  gene_plot <- ggplot(gene_data_b, aes(x = Sample, y = Expression, color = Status)) +
    geom_point() +
    labs(title = paste("Expression of", gene),
         x = "Sample", y = "Expression") +
    theme_minimal()
  
  # Save the plot to the list
  gene_plots_b[[gene]] <- gene_plot
  
  # Save or display the plot as needed
  # ggsave(paste("gene_", gene, "_plot.png"), plot = gene_plot, width = 10, height = 6)  # Uncomment to save the plots
}

# Access individual plots from the list (e.g., the first plot)
print(gene_plots_b[[3]])

##################################################################################################################

# Initialize a list to store box plots
box_plots_b <- list()

# Iterate over each gene and create separate box plots
for (gene in colnames(gene_expression_data_b)) {
  # Create a dataframe for the current gene
  gene_data_b <- data.frame(Gene = rep(gene, nrow(gene_expression_data_b)),
                            Expression = gene_expression_data_b[, gene], 
                            Status = status_column_b)
  
  # Plot using ggplot2
  box_plot <- ggplot(gene_data_b, aes(x = Status, y = Expression, fill = Status)) +
    geom_boxplot() +
    labs(title = paste("Expression of", gene),
         x = "Status", y = "Expression") +
    theme_minimal()
  
  # Save the box plot to the list
  box_plots_b[[gene]] <- box_plot
  
  # Save or display the plot as needed
  # ggsave(paste("box_", gene, "_plot.png"), plot = box_plot, width = 10, height = 6)  # Uncomment to save the box plots
}

# Access individual box plots from the list (e.g., the first plot)
print(box_plots_b[[3]])
###############################################################################################################
mmp11_data <- data.frame(Gene = rep("MMP11", nrow(gene_expression_data_b)),
                         Expression = gene_expression_data_b[, "MMP11"], 
                         Status = status_column)

# Extract expression data for positive and negative status
positive_data <- mmp11_data$Expression[mmp11_data$Status == "Positive"]
negative_data <- mmp11_data$Expression[mmp11_data$Status == "Negative"]

# Perform Shapiro-Wilk test for normality
shapiro_positive <- shapiro.test(positive_data)
shapiro_negative <- shapiro.test(negative_data)

# Print the test results
print("Shapiro-Wilk Test for Normality - Positive Status:")
print(shapiro_positive)
print("Shapiro-Wilk Test for Normality - Negative Status:")
print(shapiro_negative)

# Create histograms
hist_positive <- ggplot(data.frame(Expression = positive_data), aes(x = Expression)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "MMP11 - Positive Status",
       x = "Expression Level",
       y = "Frequency")

hist_negative <- ggplot(data.frame(Expression = negative_data), aes(x = Expression)) +
  geom_histogram(binwidth = 1, fill = "salmon", color = "black", alpha = 0.7) +
  labs(title = "MMP11 - Negative Status",
       x = "Expression Level",
       y = "Frequency")

# Display histograms
print(hist_positive)
print(hist_negative)

# Perform Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(positive_data, negative_data)
print(wilcox_test_result)

##################################################################################################################

# Create an empty data frame to store correlation results
correlation_results_b <- data.frame(Gene = character(), Correlation = numeric(), stringsAsFactors = FALSE)

# Iterate over each gene and compute correlation with status
for (gene in colnames(gene_expression_data_b)) {
  # Create a dataframe for the current gene
  gene_data_b <- data.frame(Expression = gene_expression_data_b[, gene], Status = status_column_b)
  
  # Convert 'Status' to a binary variable
  gene_data_b$Status_binary <- ifelse(gene_data_b$Status == "Positive", 1, 0)
  
  # Compute correlation between status and gene expression
  correlation_value <- cor.test(gene_data_b$Expression, gene_data_b$Status_binary)$estimate
  
  # Append results to the data frame
  correlation_results_b <- rbind(correlation_results_b, data.frame(Gene = gene, Correlation = correlation_value))
}

correlation_results_b <- correlation_results_b[order(-(correlation_results_b$Correlation)), ]

# Display the correlation table
print(correlation_results_b)

################################################################################################################

# Create an empty data frame to store correlation results
correlation_results_b2 <- data.frame(Gene = character(), MMP11_Correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Iterate over each gene and compute correlation with status
for (gene in colnames(all_gene_expression_data)) {
  # Create a dataframe for the current gene
  gene_data_b <- data.frame(Expression = all_gene_expression_data[, gene], MMP11 = mmp_data[, 4])
  
  # Compute correlation between status and gene expression
  correlation_value <- cor.test(gene_data_b$Expression, gene_data_b$MMP11)$estimate
  p_value <- cor.test(gene_data_b$Expression, gene_data_b$MMP11)$p.value
  
  # Append results to the data frame
  correlation_results_b2 <- rbind(correlation_results_b2, data.frame(Gene = gene, 
                                                                     MMP11_Correlation = correlation_value, p_value = p_value))
}

correlation_results_b2 <- correlation_results_b2[order((correlation_results_b2$MMP11_Correlation)), ]

# Display the correlation table
print(correlation_results_b2)


##################################################################################################################

# Hierarchical clustering
dist_matrix <- dist(gene_expression_data, method = "euclidean")
hierarchical_clusters <- hclust(dist_matrix, method = "ward.D2")
row_dend <- as.dendrogram(hierarchical_clusters)

# Create a heatmap
heatmap_data <- data.matrix(gene_expression_data)
colnames(heatmap_data) <- NULL  # Remove gene names for simplicity

# Specify colors for the status column
status_colors <- ifelse(hox_data$status == "Positive", "red", "blue")

# Plot the heatmap with hierarchical clustering
pheatmap(
  heatmap_data,
  clustering_distance_rows = dist_matrix,
  clustering_distance_cols = dist(heatmap_data),
  clustering_method = "ward.D2",
  annotation_row = data.frame(Status = hox_data$status),
  annotation_colors = list(Status = status_colors),
  main = "Hierarchical Clustering of hox_data"
)
################################################################################