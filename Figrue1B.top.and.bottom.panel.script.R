library(tidyverse)
library(ggplot2)
library(ggExtra)
library(gridExtra)
options(scipen=999)
library(data.table)


########### figure 1B. @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

df4 <- fread("ChIP-Seq_1b.top.10.kb.chr4G.7A.10x.txt", header = T, check.names = T, data.table = F)

colnames(df4)[1:6] = c('raw_count', 'chr', 'pos', 'read.mapped', 'nread', 'sample_type')
head(df4)

df4 <- df4[ which(df4$chr=='chr7A_TA2804'), ] ## this was just to plot chr7A 
head(df4)

df6 <- fread("ReSeq_1b.bottom.10.kb.chr4G.7A.10x.txt", header = F, check.names = T, data.table = F)

head(df6)
colnames(df6)[1:6] = c('raw_count', 'chr', 'pos', 'read.mapped','nread', 'sample_type')
head(df6)

# df7 <- df6[ which(df6$pos > 160000000 & df6$pos < 220000000 & df6$chr=='chr5D'), ]
df6 <- df6[ which(df6$chr=='chr7A_TA2804'), ] ## this was just to plot chr7A 

head(df6)

plot1 <- ggplot(data = subset(df6, df6$nread<35), aes(x=pos, y=nread, colour=read.mapped)) + 
  geom_point(size=0.8) + ## changed from 0.5 to 0.8
  geom_point(data = subset(df4, df4$nread<60), aes(x=pos, y=nread, colour=read.mapped), size=0.8)+ 
  facet_wrap( ~ chr, ncol=1, strip.position = "right" ) +   ## ,scales = "free_x"
  scale_x_continuous(labels=function(x)x/1000000, breaks = seq(0, 9000000000, by = 100000000),expand = c(0.01, 0)) +
  scale_color_manual(values = c("Cenh3" = '#690B22', "TA877_chr7A.reads" = '#0d47a1', "TA877_chr4G.reads" = '#ff6d00', "TA877_chr5A.reads" = '#1976d2')) + ## alpha =1 means fully opaque and alpha =0 means fully transparent
  theme(
    strip.text = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill="white"),
    #panel.spacing = unit(1.5, "lines"),
  )

ggsave("Full.arm.bottom.and.top.combined.png", plot1,
       width = 15.94,
       height = 5.18,
       units = "in",
       dpi = 400,
       bg = "white",
       limitsize = FALSE)





### zoom in figure 1B#####################################################################################################################################


library(data.table)
library(dplyr)
library(ggplot2)


df4_filtered <- df4 %>%
  filter(
    #  (chr == "chr4G_TA2804" & pos >= 355000000 & pos <= 390000000) |
    (chr == "chr7A_TA2804" & pos >= 354000000 & pos <= 384000000)
  )

head(df4_filtered)

df6_filtered <- df6 %>%
  filter(
    #  (chr == "chr4G_TA2804" & pos >= 355000000 & pos <= 390000000) |
    (chr == "chr7A_TA2804" & pos >= 354000000 & pos <= 384000000)
  )

head(df6_filtered)




## zoom in 1B. top.panel..CENH3

plot2 <- ggplot(data = subset(df4_filtered, df4_filtered$nread<60), aes(x=pos, y=nread, colour=read.mapped)) + 
  geom_point(size=3) +
  scale_x_continuous(labels=function(x)x/1000000, breaks = seq(354000000, 384000000, by = 1000000),expand = c(0.01, 0)) +
  scale_color_manual(values = c("Cenh3" = scales::alpha('#690B22', 1), "TA877_chr7A.reads" = scales::alpha('#0d47a1', 0), "TA877_chr4G.reads" = scales::alpha('#ff6d00', 0), "TA877_chr5A.reads" = scales::alpha('#1976d2', 0))) + ## alpha =1 means fully opaque and alpha =0 means fully transparent
  
  ###scales::alpha('red', 0.5)
  
  theme(
    
    strip.text = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill="white"),    
  )

ggsave("Figure.1B.top.cenh3.png", plot2,
       width = 15.94,
       height = 5.18,
       units = "in",
       dpi = 400,
       bg = "white",
       limitsize = FALSE)


########### zoom in 1B. bottom.panel..CENH3 @@@@@@@@@@@@@@@@@@@@@@@@@@@@
plot3 <- ggplot(
  data = subset(df6_filtered, df6_filtered$nread < 35), 
  aes(x = pos, y = nread, colour = read.mapped)
) + 
  geom_point(size = 3) +
  scale_x_continuous(
    labels = function(x) x/1000000,
    breaks = seq(354000000, 384000000, by = 1000000),
    expand = c(0.01, 0)
  ) +
  scale_color_manual(
    values = c(
      "TA877_chr7A.reads" = '#0d47a1',  # Blue (now visible)
      "TA877_chr4G.reads" = '#ff6d00' 
    ),
    # Optional: Add labels for legend clarity
    labels = c(
      "TA877_chr7A.reads" = "chr7A Reads",
      "TA877_chr4G.reads" = "chr4G Reads"
    )
  ) +
  theme(
    strip.text = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "white")
  )

ggsave("figure1b.bottom.reseq.png", plot3, 
       width = 12, height = 5, dpi = 400, bg = "white")












