

################ Script modify for publication and labeling each track #############
# Clear any existing circular layout
circos.clear()

library(ggplot2)

library(circlize)

# Text color for the labels
col_text <- "black"

# Set global parameters for the plot
circos.par("track.height" = 0.4, 
           gap.degree = c(rep(2.5, 13), 10, rep(2.5, 14)),  # Adjust the gap for separation between TA877 and TA2804
           cell.padding = c(0, 0, 0, 0))

# Define chromosome colors: blue for TA877, orange for TA2804
chr_colors <- c("chr1A_TA877" = "blue3", "chr2A_TA877" = "blue3", "chr3A_TA877" = "blue3", "chr4A_TA877" = "blue3",
                "chr5A_TA877" = "blue3", "chr6A_TA877" = "blue3", "chr7A_TA877" = "blue3", "chr1G_TA877" = "#1e97ff",
                "chr2G_TA877" = "#1e97ff", "chr3G_TA877" = "#1e97ff", "chr4G_TA877" = "#1e97ff", "chr5G_TA877" = "#1e97ff",
                "chr6G_TA877" = "#1e97ff", "chr7G_TA877" = "#1e97ff",
                "chr7G_TA2804" = "#FFB668", "chr6G_TA2804" = "#FFB668",
                "chr5G_TA2804" = "#FFB668", "chr4G_TA2804" = "#FFB668", "chr3G_TA2804" = "#FFB668", "chr2G_TA2804" = "#FFB668",
                "chr1G_TA2804" = "#FFB668", "chr7A_TA2804" = "orange3", "chr6A_TA2804" = "orange3", "chr5A_TA2804" = "orange3",
                "chr4A_TA2804" = "orange3", "chr3A_TA2804" = "orange3", "chr2A_TA2804" = "orange3", "chr1A_TA2804" = "orange3")

# Initialize circular layout with chromosome limits and specific colors
circos.initialize(factors = c("chr1A_TA877", "chr2A_TA877", "chr3A_TA877", "chr4A_TA877", "chr5A_TA877", "chr6A_TA877",  "chr7A_TA877", 
                              "chr1G_TA877",  "chr2G_TA877",  "chr3G_TA877","chr4G_TA877",  "chr5G_TA877",  "chr6G_TA877", "chr7G_TA877",
                              
                              "chr7G_TA2804", "chr6G_TA2804", "chr5G_TA2804", "chr4G_TA2804", "chr3G_TA2804", "chr2G_TA2804", "chr1G_TA2804", 
                              "chr7A_TA2804",  "chr6A_TA2804",  "chr5A_TA2804","chr4A_TA2804",  "chr3A_TA2804",  "chr2A_TA2804","chr1A_TA2804"),
                 xlim = matrix(c(rep(0, 28), c(602656826, 771401179, 675871483, 773897100, 691724660,  595162121, 728204638, 
                                               499920424, 683726055, 670590195, 659442802,  662613401, 593054900, 715925247,
                                               
                                               696135035, 599005360, 649451578, 726833302, 667333744, 689748856, 497037353,
                                               746009334,  582428954,  614132190, 772671662,  670699037,  768562081,  613121534)), ncol = 2))

# Track for chromosome names, displaying simplified names
circos.track(ylim = c(0, 3), ## changed from 1 to 3
             panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               display_chr <- sub("_.*", "", chr)
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               circos.text(mean(xlim), mean(ylim), display_chr, cex = 0.95, col = "white",
                           font = 2, # bold font
                           facing = "bending.inside", niceFacing = TRUE)
             },
             bg.col = chr_colors[as.character(get.all.sector.index())],  # Set background color for each chromosome
             bg.border = FALSE,
             track.height = 0.1) 
            

# Set breakpoints for the axis
brk <- c(0, 200, 400, 600, 800) * 10^6 ## changed it from 0, 100, 200 to 800 to this style to increase font size

# Track for chromosome axis, same for all chromosomes
circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
               circos.axis(h = "top",
                           major.at = brk,
                           labels = round(brk / 10^6, 1),
                           labels.cex = 0.8,
                           col = col_text,
                           labels.col = col_text,
                           lwd = 0.9,
                           labels.facing = "clockwise")
             },
             bg.border = FALSE)



### Density plots and genomic tracks

# High-confidence genes density plot
circos.genomicTrack(data = read.table("HC_both.genome.gene.density.txt", header = TRUE), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, 1 * (value), type = "h", area = TRUE, col = "grey10", lwd = 0.7)
                    }, 
                    track.height = 0.06, bg.border = FALSE)

# Gypsy-TE density plot
circos.genomicTrack(data = read.table("both_TA2804.TA877.Gypsy.10Mb.density.txt", header = TRUE),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, 1 * (value), type = "h", area = TRUE, col = "pink4", lwd = 0.7)
                    }, 
                    track.height = 0.06, bg.border = FALSE)

# Copia-TE elements density plot
circos.genomicTrack(data = read.table("both_TA2804.TA877.Copia.10Mb.density.txt", header = TRUE),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, 1 * (value), type = "h", area = TRUE, col = "magenta", lwd = 0.7)
                    }, 
                    track.height = 0.06, bg.border = FALSE)

# CENh3 density plot
circos.genomicTrack(data = read.table("combined.CENh3.TA2804.mapping.for.circos.txt", header = TRUE),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, 1 * (value), type = "h", area = TRUE, col = "brown", lwd = 0.7)
                    }, 
                    track.height = 0.06, bg.border = FALSE)

# Synteny plot between homeologous chromosomes
syntenic_data <- read.table("filtered_matching_homeologoul_translocation.5A.txt", header = TRUE)

# Define color mapping for different chromosome pairs
# Define color mapping for different chromosome pairs
color_mapping <- c(
  "chr1A_TA2804-chr1A_TA877" = "#D01C8B", "chr1G_TA2804-chr1G_TA877" = "#D7191C",
  "chr2A_TA2804-chr2A_TA877" = "#FDAE61", "chr2G_TA2804-chr2G_TA877" = "#F4A582",
  "chr3A_TA2804-chr3A_TA877" = "#FFEDA0", "chr3G_TA2804-chr3G_TA877" = "#FFFFCC",
  "chr4A_TA2804-chr4A_TA877" = "#D9EF8B", "chr4G_TA2804-chr4G_TA877" = "purple4",
  "chr5A_TA2804-chr5A_TA877" = "orange", "chr5G_TA2804-chr5G_TA877" = "#2C7BB6",
  "chr6A_TA2804-chr6A_TA877" = "#ABD9E9", "chr6G_TA2804-chr6G_TA877" = "#67A9CF",
  "chr7A_TA2804-chr7A_TA877" = "blue2", "chr7G_TA2804-chr7G_TA877" = "#0571B0",
  "chr7A_TA2804-chr4G_TA877" = "purple4", "chr4G_TA2804-chr7A_TA877" = "blue2",
  "chr4G_TA2804-chr5A_TA877" = "orange"
)


# Draw links between orthologues
for (i in 1:nrow(syntenic_data)) {
  chr1 <- syntenic_data[i, "chr1"]
  chr2 <- syntenic_data[i, "chr2"]
  region1 <- syntenic_data[i, 2:3]  # Start1 and End1
  region2 <- syntenic_data[i, 5:6]  # Start2 and End2
  
  # Create key for color mapping
  key <- paste(chr1, chr2, sep = "-")
  
  
  # Check if key exists in color_mapping, otherwise assign a default color
  if (key %in% names(color_mapping)) {
    col <- color_mapping[key]
  } else {
    col <- "black"  # Default color if the key is not found
  }
  
  # Draw genomic link with specific color and thin line
  circos.genomicLink(region1 = data.frame(syntenic_data[i, 1:3]), 
                     region2 = data.frame(syntenic_data[i, 4:6]), 
                     col = color_mapping[key], 
                     lwd = 0.15)  # Thin line (can adjust if needed) replaced 0.5 by 0.15
}

# Clear the plot when done
circos.clear()


