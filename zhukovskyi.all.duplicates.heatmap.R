######### dual length allele matching  ####################################


library(ggplot2)
library(reshape2)
library(ggnewscale) 
library(svglite)# For multiple fill scales

data <- read.csv("zhukovskyi_duplicates.arranged.csv", row.names = 1)
matrix_data <- as.matrix(data)
melted_data <- melt(matrix_data)
colnames(melted_data) <- c("Row", "Column", "Value")

# Prepare data subsets
data_low <- subset(melted_data, Value <= 0.98)
data_high <- subset(melted_data, Value > 0.98)

ggplot() +
  # Lower range (0-0.95) with blue gradient
  geom_tile(data = data_low, aes(x = Column, y = Row, fill = Value)) +
  scale_fill_gradientn(
    name = " ",
    colours = c("white", "#F9C74F"),
    limits = c(0, 0.98),
    breaks = seq(0, 0.98, 0.14),
    na.value = "grey"
  ) +
  # New fill scale for high values
  ggnewscale::new_scale_fill() +
  # Upper range (0.95-1) with red-gold gradient
  geom_tile(data = data_high, aes(x = Column, y = Row, fill = Value)) +
  scale_fill_gradientn(
    name = " ", ## put nothing in it title
    colours = c("#61F4DE", "#6E78FF"),
    limits = c(0.98, 1),
    breaks = seq(0.98, 1, 0.005),
    labels = c("0.98", "0.985","0.99", "0.995", "1.00"),
    na.value = "grey"
  ) +
  # Plot styling
  labs(title = "Dual Gradient Heatmap", x = "Columns", y = "Rows") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    legend.box = "vertical"
  )

# Save with adequate dimensions
ggsave("dual_legend_heatmap3.final.svg", width = 14, height = 11, dpi = 300)







