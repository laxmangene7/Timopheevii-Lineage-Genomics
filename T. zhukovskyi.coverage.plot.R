library(tidyverse)
library(ggplot2)
library(data.table)

##### Script to generate plot in Figure 2 of Manuscript

bed_files <- list.files(pattern = "\\.bed$") ## The bed files were produced with the script provided https://github.com/sandeshsth/SkimSeq_Method or [Adhikari et al. 2022; https://www.nature.com/articles/s41598-022-19858-2] 

library(tidyverse)
library(ggplot2)
library(data.table)

bed_files <- list.files(pattern = "\\.bed$")

for (input_file in bed_files) {
  # Read data
  df4 <- fread(input_file, header = FALSE, check.names = TRUE, data.table = FALSE)
  colnames(df4)[1:4] <- c('chr', 'start', 'end', 'nread')
  
  sample_name <- tools::file_path_sans_ext(basename(input_file))
  
  # Calculate chromosome sizes (max end position per chromosome)
  chr_sizes <- aggregate(end ~ chr, data = df4, max)
  
  # Merge sizes back to the main data
  df4 <- merge(df4, chr_sizes, by = "chr", suffixes = c("", "_max"))
  
  # Filter data
  df4_subset <- subset(df4, nread < 10)
  
  # Define base color labels
  copy_levels <- c("0 copies", "1 copy", "2 copies", "3 copies", "4 copies")
  df4_subset$color_bin_base <- cut(
    df4_subset$nread,
    breaks = c(0, 0.25, 0.75, 1.25, 1.75, 2.25),
    labels = copy_levels,
    include.lowest = TRUE
  )
  

  
  
  # Assign conditional color labels
  df4_subset$color_bin <- ifelse(
    grepl("TA2804", df4_subset$chr),
    paste("TA2804", df4_subset$color_bin_base),
    ifelse(
      grepl("TA10622", df4_subset$chr),
      paste("TA10622", df4_subset$color_bin_base),
      paste("OTHER", df4_subset$color_bin_base)
    )
  )
  
  
  
  # Define palettes
  copy_number_colors_normal <- c(
    "0 copies" = "#a9c8e1ff",
    "1 copy" = "#61a5c2ff",
    "2 copies" = "#2c7da0ff",
    "3 copies" = "#2a6f97ff",
    "4 copies" = "#01497cff"
  )
  copy_number_colors_reversed <- rev(copy_number_colors_normal)
  copy_levels_rev <- rev(copy_levels)
  # Combined color map
  copy_number_colors_all <- c(
    setNames(copy_number_colors_normal, paste("TA2804", copy_levels)),
    setNames(copy_number_colors_reversed, paste("TA10622", copy_levels))#,
    #setNames(copy_number_colors_normal, paste0("OTHER_", copy_levels))
  )
  
  print(copy_number_colors_all)
  
  #copy_levels_custom <- c(         
  
  #copy_levels_custom <- c("TA2804 3 copies", "TA2804 2 copies", "TA2804 1 copy", "TA2804 0 copies", "TA2804 4 copies", "TA10622 3 copies", "TA10622 2 copies", "TA10622 1 copy", "TA10622 0 copies", "TA10622 4 copies")
  #copy_number_colors_custom <- c("#a9c8e1ff", "#61a5c2ff", "#2c7da0ff", "#2a6f97ff", "#01497cff", "#01497cff", "#2a6f97ff", "#2c7da0ff", "#61a5c2ff","#a9c8e1ff")
  copy_levels_custom <- c("TA2804 0 copies", "TA2804 1 copy", "TA2804 2 copies", "TA2804 3 copies", "TA2804 4 copies", "TA10622 0 copies", "TA10622 1 copy", "TA10622 2 copies", "TA10622 3 copies", "TA10622 4 copies")
  
  df4_subset$color_bin <- factor(df4_subset$color_bin, levels = names(copy_number_colors_all))
  
  
  sample_name <- tools::file_path_sans_ext(basename(input_file))
  
  # Plot
  p <- ggplot(df4_subset, aes(x = start, y = nread, color = color_bin)) +
    geom_point(size = 1.8) +
    facet_wrap(~ chr, ncol = 1, scales = "free_x", strip.position = "right") +
    scale_x_continuous(labels = function(x) x / 1e6, expand = c(0.01, 0)) +
    scale_color_manual(
      name = "Copy number",
      values = copy_number_colors_all,
      drop = FALSE,
      labels = copy_levels_custom
    ) +
    xlab("\nGenomic position (Mb)") +
    ylab("Normalised coverage (x)\n") +
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
  
  
  p <- p + ggtitle(sample_name) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  

  if (!exists("pdf_started")) {
  pdf("combined_output.pdf", width = 13.5, height = 45)
  assign("pdf_started", TRUE, envir = .GlobalEnv)
  }
  
  # Print the plot to the current PDF device
  print(p)
  
}

#message("All files processed successfully!")
if (exists("pdf_started")) dev.off()



