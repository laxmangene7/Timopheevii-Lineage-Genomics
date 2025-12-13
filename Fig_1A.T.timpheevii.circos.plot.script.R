library(ggplot2)
library(circlize)

# ------------------------------------------------------------------
# OPEN THE SVG DEVICE FIRST
# ------------------------------------------------------------------
svg("circos_plot_updated.11.24.2025.svg", width = 11, height = 11)

# ------------------------------------------------------------------
# REDUCE OUTER MARGINS FOR GOOD LABEL POSITIONS
# ------------------------------------------------------------------
par(mai = c(0.2, 0.6, 0.3, 0.6)) 
# bottom, left, top, right

# Clear any previous circos
circos.clear()

# Text color for axes/labels
col_text <- "black"

# Chromosome order
chromosomes <- c(
  "chr1A_TA877","chr2A_TA877","chr3A_TA877","chr4A_TA877","chr5A_TA877","chr6A_TA877","chr7A_TA877",
  "chr1G_TA877","chr2G_TA877","chr3G_TA877","chr4G_TA877","chr5G_TA877","chr6G_TA877","chr7G_TA877",
  "chr7G_TA2804","chr6G_TA2804","chr5G_TA2804","chr4G_TA2804","chr3G_TA2804","chr2G_TA2804","chr1G_TA2804",
  "chr7A_TA2804","chr6A_TA2804","chr5A_TA2804","chr4A_TA2804","chr3A_TA2804","chr2A_TA2804","chr1A_TA2804",
  "chr1A","chr2A","chr3A","chr4A","chr5A","chr6A","chr7A",
  "chr1G","chr2G","chr3G","chr4G","chr5G","chr6G","chr7G"
)

# Gaps
default_gap <- 2
gaps <- rep(default_gap, length(chromosomes))

gaps[which(chromosomes == "chr7G_TA877")]   <- 8
gaps[which(chromosomes == "chr1A_TA2804")]  <- 6
gaps[which(chromosomes == "chr7G")]         <- 6

# Chrom lengths (your values)
xlim_vals <- matrix(c(rep(0, 42), c(
  602656826, 771401179, 675871483, 773897100, 691724660, 595162121, 728204638,
  499920424, 683726055, 670590195, 659442802, 662613401, 593054900, 715925247,
  696135035, 599005360, 649451578, 726833302, 667333744, 689748856, 497037353,
  746009334, 582428954, 614132190, 772671662, 670699037, 768562081, 613121534,
  606421902, 766211424, 672715223, 764583328, 684654349, 591560360, 734384476,
  508095237, 680075684, 840473417, 650379239, 479247860, 584502105, 697235649
)), ncol = 2)

# Circos parameters
circos.par(
  track.height = 0.4,
  gap.degree = gaps,
  cell.padding = c(0,0,0,0),
  start.degree = 30
)

# Colors for each chromosome
chr_colors <- c(
  "chr1A_TA877"="orange2","chr2A_TA877"="orange2","chr3A_TA877"="orange2","chr4A_TA877"="orange2",
  "chr5A_TA877"="orange2","chr6A_TA877"="orange2","chr7A_TA877"="orange2",
  "chr1G_TA877"="#FFB668","chr2G_TA877"="#FFB668","chr3G_TA877"="#FFB668","chr4G_TA877"="#FFB668",
  "chr5G_TA877"="#FFB668","chr6G_TA877"="#FFB668","chr7G_TA877"="#FFB668",
  "chr7G_TA2804"="#abc2ff","chr6G_TA2804"="#abc2ff","chr5G_TA2804"="#abc2ff","chr4G_TA2804"="#abc2ff",
  "chr3G_TA2804"="#abc2ff","chr2G_TA2804"="#abc2ff","chr1G_TA2804"="#abc2ff",
  "chr7A_TA2804"="#00b3d8","chr6A_TA2804"="#00b3d8","chr5A_TA2804"="#00b3d8","chr4A_TA2804"="#00b3d8",
  "chr3A_TA2804"="#00b3d8","chr2A_TA2804"="#00b3d8","chr1A_TA2804"="#00b3d8",
  "chr1A"="#266ed9","chr2A"="#266ed9","chr3A"="#266ed9","chr4A"="#266ed9","chr5A"="#266ed9","chr6A"="#266ed9","chr7A"="#266ed9",
  "chr1G"="#1e97ff","chr2G"="#1e97ff","chr3G"="#1e97ff","chr4G"="#1e97ff","chr5G"="#1e97ff","chr6G"="#1e97ff","chr7G"="#1e97ff"
)

# Initialize
circos.initialize(factors = chromosomes, xlim = xlim_vals)

# ------------------------------------------------------------------
# CHROMOSOME LABELS
# ------------------------------------------------------------------
circos.track(
  ylim = c(0, 3),
  bg.col = chr_colors[chromosomes],
  bg.border = FALSE,
  track.height = 0.1,
  panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    label <- gsub("chr|_.*", "", chr)
    circos.text(
      mean(CELL_META$xlim), mean(CELL_META$ylim),
      label,
      cex = 1.3, col = "white", font = 2,
      facing = "bending.inside", niceFacing = TRUE
    )
  }
)

# ------------------------------------------------------------------
# GENOMIC TRACKS
# ------------------------------------------------------------------

# Gene density
circos.genomicTrack(
  read.table("HC_three.genomes.gene.density.txt", header = TRUE),
  track.height = 0.06, bg.border = FALSE,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type = "h", area = TRUE,
                        col = "grey10", lwd = 0.7)
  }
)

# Gypsy
circos.genomicTrack(
  read.table("three.genomes.Gypsy.txt", header = TRUE),
  track.height = 0.06, bg.border = FALSE,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type = "h", area = TRUE,
                        col = "pink4", lwd = 0.7)
  }
)

# Copia
circos.genomicTrack(
  read.table("three.genomes.Copia.txt", header = TRUE),
  track.height = 0.06, bg.border = FALSE,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type = "h", area = TRUE,
                        col = "magenta", lwd = 0.7)
  }
)

# CENH3
circos.genomicTrack(
  read.table("combined.CENh3.TA2804.mapping.for.circos.txt", header = TRUE),
  track.height = 0.06, bg.border = FALSE,
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, type = "h", area = TRUE,
                        col = "brown", lwd = 0.7)
  }
)

# ------------------------------------------------------------------
# SYNTENY LINKS
# ------------------------------------------------------------------

syntenic_data <- read.table("only.translocation.txt", header = TRUE)

color_mapping <- c(
  "chr1A_TA2804-chr1A_TA877" = "#D01C8B", "chr1G_TA2804-chr1G_TA877" = "#D7191C",
  "chr2A_TA2804-chr2A_TA877" = "#FDAE61", "chr2G_TA2804-chr2G_TA877" = "#F4A582",
  "chr3A_TA2804-chr3A_TA877" = "#FFEDA0", "chr3G_TA2804-chr3G_TA877" = "#FFFFCC",
  "chr4A_TA2804-chr4A_TA877" = "#D9EF8B", "chr4G_TA2804-chr4G_TA877" = "purple4",
  "chr5A_TA2804-chr5A_TA877" = "orange", "chr5G_TA2804-chr5G_TA877" = "#2C7BB6",
  "chr6A_TA2804-chr6A_TA877" = "#ABD9E9", "chr6G_TA2804-chr6G_TA877" = "#67A9CF",
  "chr7A_TA2804-chr7A_TA877" = "blue2", "chr7G_TA2804-chr7G_TA877" = "#0571B0",
  
  "chr7A_TA2804-chr4G_TA877" = "#00b3d8", "chr4G_TA2804-chr7A_TA877" = "#abc2ff",
  "chr4G_TA2804-chr5A_TA877" = "#abc2ff",  "chr5A_TA2804-chr4G_TA877" = "#00b3d8",
  
  "chr1A_TA2804-chr1A" = "#D01C2B", "chr1G_TA2804-chr1G" = "#D7118C",
  "chr2A_TA2804-chr2A" = "#FDAE21", "chr2G_TA2804-chr2G" = "#F4A532",
  "chr3A_TA2804-chr3A" = "#FFEDA5", "chr3G_TA2804-chr3G" = "#FFFFCB",
  "chr4A_TA2804-chr4A" = "#D9EF3B", "chr4G_TA2804-chr4G" = "purple2",
  "chr5A_TA2804-chr5A" = "orange2", "chr5G_TA2804-chr5G" = "#2C7BB3",
  "chr6A_TA2804-chr6A" = "#ABD9E9", "chr6G_TA2804-chr6G" = "#67A4CF",
  "chr7A_TA2804-chr7A" = "blue1", "chr7G_TA2804-chr7G" = "#0571B2",
  
  "chr5A_TA2804-chr4G" = "#266ed9", "chr4G_TA2804-chr5A" = "#1e97ff",
  "chr3G_TA2804-chr5G" = "#1e97ff", "chr5G_TA2804-chr3G" = "#1e97ff",
  
  "chr1A_TA877-chr1A" = "#D01C2B", "chr1G_TA877-chr1G" = "#D7118C",
  "chr2A_TA877-chr2A" = "#FDAE21", "chr2G_TA877-chr2G" = "#F4A532",
  "chr3A_TA877-chr3A" = "#FFEDA5", "chr3G_TA877-chr3G" = "#FFFFCB",
  "chr4A_TA877-chr4A" = "#D9EF3B", "chr4G_TA877-chr4G" = "purple2",
  "chr5A_TA877-chr5A" = "orange2", "chr5G_TA877-chr5G" = "#2C7BB3",
  "chr6A_TA877-chr6A" = "#ABD9E9", "chr6G_TA877-chr6G" = "#67A4CF",
  "chr7A_TA877-chr7A" = "blue1", "chr7G_TA877-chr7G" = "#0571B2",
  
  "chr7A_TA877-chr4G" = "orange2", "chr4G_TA877-chr7A" = "#FFB668",
  "chr4G_TA877-chr5A" = "#FFB668", "chr5G_TA877-chr3G" = "#FFB668",
  "chr3G_TA877-chr5G" = "#FFB668"
)



for (i in 1:nrow(syntenic_data)) {
  chr1 <- syntenic_data$chr1[i]
  chr2 <- syntenic_data$chr2[i]
  key <- paste(chr1, chr2, sep = "-")
  col <- ifelse(key %in% names(color_mapping),
                color_mapping[key], "black")
  circos.genomicLink(
    region1 = syntenic_data[i, 1:3],
    region2 = syntenic_data[i, 4:6],
    col = col,
    lwd = 0.15
  )
}

# ------------------------------------------------------------------
# TITLES / LEGENDS
# ------------------------------------------------------------------

#title(main = "KU-1948B (ARA-1)", cex.main = 2, font.main = 2, line = -1)

# Place TA2804 LEFT side
#mtext("TA2804 (Domesticated)", side = 2, line = 0.2, cex = 1.8, font = 2)

# Place TA877 RIGHT side
#mtext("TA877 (ARA-0)", side = 4, line = 0.2, cex = 1.8, font = 2)

title(main = "KU-1948B (ARA-1)", cex.main = 2, font.main = 2, line = -1)

# Place TA2804 LEFT side, smaller font, closer to circos
mtext("TA2804 (Domesticated)", side = 2, line = 0, cex = 1.5, font = 2)

# Place TA877 RIGHT side, smaller font, closer to circos
mtext("TA877 (ARA-0)", side = 4, line = 0, cex = 1.5, font = 2)



# ------------------------------------------------------------------
# End
# ------------------------------------------------------------------
dev.off()
circos.clear()


