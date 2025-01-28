# Load necessary libraries
library(readxl)
library(corrplot)
library(Hmisc)  # For correlation test to get p-values

# Set the working directory to the folder where the script is located
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

##### All samples correlation####
# Load the data from the Excel file in the working directory
file_path <- "BL_OT_Combined.xlsx"
data <- read_excel(file_path, sheet = "Sheet1") #specifies sheet in xlsx file

# Remove rows with specific pt id values
data <- data[!data$`pt_id` %in% c("286_0", "250_0", "325_0", "535_0", "580_0", "612_0", "653_0","110_0", "833_0", "127_0", "257_0", "795_0", "497_0"), ]
data <- data[!data$`pt_id` %in% c("286_1", "580_1", "127_1", "497_1"), ]

# Select cytokine and cell cluster columns
cytokine_data <- data[, c("scd40l", "il-1a", "il-1b","il-1ra","il-2","il-3","il-4","il-5","il-6","il-8","il-9","il-10","il-12(p40)","il-12(p70)","il-13","il-15","il-17a","il-25","il-17f","il-18","il-22","g-csf","gm-csf","tnfa","ifng","mcp-1","mip-1a","mip-1b","rantes","mig","ip-10","vegf-a")]
cell_clusters <- data[, c("B", "DNT", "Myeloid", "NK", "TcCM", "TcEFF", "TcEM", "TcN", "Th17", "Th2", "Th2CM", "Th2EM", "ThCTL", "ThN", "Treg", "UA")]

# Remove rows with NA values in either dataset
complete_cases <- complete.cases(cytokine_data, cell_clusters)
cytokine_data <- cytokine_data[complete_cases, ]
cell_clusters <- cell_clusters[complete_cases, ]

# Combine cytokine and cell cluster data for correlation
combined_data <- cbind(cytokine_data, cell_clusters)

# Calculate correlation matrix and p-values for cytokine to cell type correlations
correlation_results <- rcorr(as.matrix(combined_data))

# Extract correlation matrix and p-values
cor_matrix <- correlation_results$r
p_values <- correlation_results$P

# Set up the PDF file to save the plot in the working directory
pdf(file = file.path(script_dir, "correlation_plot_box_BL_OT_Combo.pdf"), width = 10, height = 10)

# Create a plot focusing on cytokine-to-cell type correlations
corrplot(cor_matrix[1:ncol(cytokine_data), (ncol(cytokine_data)+1):ncol(combined_data)],  # Subset matrix for cytokine-cell correlations
         method = "circle", 
         type = "full", 
         tl.col = "black", 
         tl.cex = 0.8,
         addgrid.col = "gray",  # Light grid lines
         cl.lim = c(-1, 1),  # Set color scale limits
         col = colorRampPalette(c("blue", "white", "red"))(200))  # Color scale: blue (negative) to red (positive)
dev.off()

# Save correlation values (cytokine-to-cell type) as a CSV file
write.csv(cor_matrix[1:ncol(cytokine_data), (ncol(cytokine_data)+1):ncol(combined_data)], 
          file = file.path(script_dir, "correlation_values_box_cytokine_celltype.csv"), row.names = TRUE)

