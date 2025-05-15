library(tidyverse)
library(ggplot2)
library(data.table)

##### Script to generate plot in Figure 2 of Manuscript

bed_files <- list.files(pattern = "\\.bed$") ## The bed files were produced with the script provided https://github.com/sandeshsth/SkimSeq_Method or [Adhikari et al. 2022; https://www.nature.com/articles/s41598-022-19858-2] 

# Process each file
for (input_file in bed_files) {
  # Read data
  df4 <- fread(input_file, header = FALSE, check.names = TRUE, data.table = FALSE)
  colnames(df4)[1:4] <- c('chr', 'start', 'end', 'nread')
  
  # Create color bins with new labels
  df4$color_bin <- cut(df4$nread,
                       breaks = c(0, 0.25, 0.75, 1.25, 1.75, 2.25),
                       labels = c("null", "1 copy", "2 copies", "3 copies", "4 copies"),
                       include.lowest = TRUE)
  
  # Calculate chromosome sizes (max end position per chromosome)
  chr_sizes <- aggregate(end ~ chr, data = df4, max)
  
  # Merge sizes back to the main data
  df4 <- merge(df4, chr_sizes, by = "chr", suffixes = c("", "_max"))
  
  # Filter data
  df4_subset <- subset(df4, nread < 10)
  
  # Create output filename (replace .bed with .svg)
  output_file <- sub("\\.bed$", ".svg", input_file)
  
  # Define color palette
  copy_number_colors <- c(
    "null" = "#7FDBFF",
    "1 copy" = "#6baed3",
    "2 copies" = "#1E99FF",
    "3 copies" = "#0000CD",
    "4 copies" = "#084594"
  )
  
  # Plot with free x-axis scales and proportional lengths
  p <- ggplot(df4_subset, aes(x = start, y = nread, color = color_bin)) +
    geom_point(size = 1.8) +  # , alpha = 0.7 Add slight transparency if needed
    facet_wrap(~ chr, ncol = 1, scales = "free_x", strip.position = "right") +
    scale_x_continuous(
      labels = function(x) x / 1e6,
      expand = c(0.01, 0)
    ) +
    scale_color_manual(
      name = "Copy Number",
      values = copy_number_colors,
      drop = FALSE
    ) +
    xlab("\nGenomic Position (Mb)") +
    ylab("Normalized coverage (x)\n") +
    theme(
      panel.grid.major.x = element_line(size = 0.08, linetype = 'solid', colour = "gray68"),
      panel.grid.major.y = element_line(size = 0.08, linetype = 'solid', colour = "gray68"),
      panel.background = element_rect(fill = "white"),
      panel.spacing = unit(1.5, "lines"),
      strip.text.y.right = element_text(angle = 0, face = "bold", size = 14),
      panel.border = element_rect(colour = "black", size = 0.4, fill = NA),
      axis.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 20),
      legend.position = "top",
      legend.justification = "center",
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.margin = margin(t = 0, r = 0, b = 10, l = 0)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))  # Larger legend dots
  
  # Save plot as SVG (scalable vector graphics)
  ggsave(
    filename = output_file,
    plot = p,
    device = "svg",
    width = 13.5, 
    height = 45,   ## change height if you want to run only few chromosomes
    units = "in", 
    limitsize = FALSE
  )
  
  # Print progress message
  message("Processed: ", input_file, " -> Saved as: ", output_file)
}

message("All files processed successfully!")
