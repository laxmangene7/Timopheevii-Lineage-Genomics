#################### sample names in each graphs ################################
library(tidyverse)
library(ggplot2)
library(data.table)

# Get list of all .txt files in the current directory
bed_files <- list.files(pattern = "\\.txt$")

# Process each file
for (input_file in bed_files) {
  # Read data
  df4 <- fread(input_file, header = FALSE, check.names = TRUE, data.table = FALSE)
  colnames(df4)[1:4] <- c('rawread', 'chr', 'pos', 'nread')
  
  # Create color bins with extended range (up to 6 copies)
  df4$color_bin <- cut(
    df4$nread,
    breaks = c(0, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.5),
    labels = c("null", "1 copy", "2 copies", "3 copies", "4 copies", "5 copies", "6 copies"),
    include.lowest = TRUE
  )
  
  # Calculate chromosome sizes (max end position per chromosome)
  chr_sizes <- aggregate(pos ~ chr, data = df4, max)
  
  # Merge sizes back to the main data
  df4 <- merge(df4, chr_sizes, by = "chr", suffixes = c("", "_max"))
  
  # Optionally filter extreme coverage values (e.g., > 7x if needed)
  df4_subset <- subset(df4, nread < 3.5)
  
  # Create output filename (replace .txt with .svg)
  output_file <- sub("\\.txt$", ".svg", input_file)
  
  # Get sample name without .txt extension for the plot title
  sample_name <- sub("\\.txt$", "", input_file)
  
  # Define extended color palette
  copy_number_colors <- c(
    "null"       = "#D0ECFF",  # very light cyan
    "1 copy"     = "#9DD9FF",  # pale sky blue
    "2 copies"   = "#6BB8FF",  # medium light blue
    "3 copies"   = "#3A8DFF",  # bright azure
    "4 copies"   = "#1C5FD3",  # strong medium-dark blue
    "5 copies"   = "#0E3C99",  # deep royal blue
    "6 copies"   = "#061E5C"   # very dark navy blue
  )
  
  # Plot with free x-axis scales and proportional lengths
  p <- ggplot(df4_subset, aes(x = pos, y = nread, color = color_bin)) +
    geom_point(size = 1.8) +
    facet_wrap(~ chr, ncol = 1, scales = "free_x", strip.position = "right") +
    scale_x_continuous(
      labels = function(x) x / 1e6,
      expand = c(0.01, 0)
    ) +
    scale_color_manual(
      name = NULL,  # Remove "Copy Number" title from legend
      values = copy_number_colors,
      drop = FALSE
    ) +
    xlab("\nGenomic Position (Mb)") +
    ylab("Normalized coverage (x)\n") +
    ggtitle(sample_name) +  # Add sample name as title
    theme(
      panel.grid.major.x = element_line(size = 0.08, linetype = 'solid', colour = "gray68"),
      panel.grid.major.y = element_line(size = 0.08, linetype = 'solid', colour = "gray68"),
      panel.background = element_rect(fill = "white"),
      panel.spacing = unit(1.5, "lines"),
      strip.text.y.right = element_text(angle = 0, face = "bold", size = 14),
      panel.border = element_rect(colour = "black", size = 0.4, fill = NA),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 20),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center the title
      legend.position = "top",
      legend.justification = "center",
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size = 12),  # Slightly smaller text for single row
      legend.title = element_blank(),  # Remove legend title
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.margin = margin(t = 0, r = 0, b = 5, l = 0),
      legend.spacing.x = unit(0.3, "cm")  # Add spacing between legend items
    ) +
    guides(color = guide_legend(
      nrow = 1,  # Force all items in single row
      override.aes = list(size = 3)  # Slightly smaller legend points
    ))
  
  # Save plot as SVG
  ggsave(
    filename = output_file,
    plot = p,
    device = "svg",
    width = 13.5,
    height = 45,
    units = "in",
    limitsize = FALSE
  )
  
  # Print progress message with sample name (without .txt)
  message("Processed: ", sample_name, " -> Saved as: ", output_file)
}

message("All files processed successfully!")





