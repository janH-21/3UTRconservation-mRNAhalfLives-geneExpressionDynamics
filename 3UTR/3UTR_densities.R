## load packages
require(tidyverse)
require(ggpubr)

## load data
hg38_pivot_half <- read.csv("dataR/hg38_pivot_half.csv")

## x scale max.
xlimit = 1250;

#### cumulative distributions of single factors ####
## by conservation
colors = c("all" = "black", "cons" = "green", "non-cons" = "blue")
consCums <- ggplot() +
  stat_ecdf(data = hg38_pivot_half, aes(x=Half.life.min, color = "all")) +
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(threeUTR_conservation > 0.6), aes(x=Half.life.min, color = "cons")) + 
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(threeUTR_conservation <= 0.6), aes(x=Half.life.min, color = "non-cons")) +
  labs(x = "half-life [min]", y = "cumulative gene density", color = "") +
  scale_color_manual(values = colors) + 
  theme(legend.position="bottom", legend.box = "horizontal", legend.text = element_text(size=12),
        plot.title = element_text(size=13, face="bold", hjust=0.5, lineheight=1.2),
        plot.subtitle = element_text(size=10, face="italic", hjust=0.5),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12)) + 
  xlim(0,xlimit)

## by length
colors = c("all" = "black", "long" = "green", "short" = "blue")
lenCums <- ggplot() +
  stat_ecdf(data = hg38_pivot_half, aes(x=Half.life.min, color = "all")) +
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(threeUTR_length > 1000), aes(x=Half.life.min, color = "long")) + 
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(threeUTR_length <= 1000), aes(x=Half.life.min, color = "short")) +
  labs(x = "half-life [min]", y = "cumulative gene density", color = "") +
  scale_color_manual(values = colors) + 
  theme(legend.position="bottom", legend.box = "horizontal", legend.text = element_text(size=12),
        plot.title = element_text(size=13, face="bold", hjust=0.5, lineheight=1.2),
        plot.subtitle = element_text(size=10, face="italic", hjust=0.5),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12)) + 
  xlim(0,xlimit)

## by transcription factors
colors = c("all" = "black", "TF" = "green", "non-TF" = "blue")
TFCums <- ggplot() +
  stat_ecdf(data = hg38_pivot_half, aes(x=Half.life.min, color = "all")) +
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(isTF == TRUE), aes(x=Half.life.min, color = "TF")) + 
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(isTF != TRUE), aes(x=Half.life.min, color = "non-TF")) +
  labs(x = "half-life [min]", y = "cumulative gene density", color = "") +
  scale_color_manual(values = colors) + 
  theme(legend.position="bottom", legend.box = "horizontal", legend.text = element_text(size=12),
        plot.title = element_text(size=13, face="bold", hjust=0.5, lineheight=1.2),
        plot.subtitle = element_text(size=10, face="italic", hjust=0.5),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12)) + 
  xlim(0,xlimit)

ggarrange(consCums, lenCums, TFCums, nrow = 1)

#### cumulative distributions of combined factors ####
colors = c("all" = "black", "long + cons TFs" = "green", "long + cons" = "royalblue4" , "all TFs" = "deepskyblue2")
cumDens <- ggplot() +
  stat_ecdf(data = hg38_pivot_half, aes(x=Half.life.min, color = "all")) +
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(threeUTR_conservation > 0.6 & threeUTR_length > 1000 &  isTF == TRUE), aes(x=Half.life.min, color = "long + cons TFs")) + 
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(threeUTR_conservation > 0.6 & threeUTR_length > 1000), aes(x=Half.life.min, color = "long + cons")) +
  stat_ecdf(data = hg38_pivot_half %>% dplyr::filter(isTF == TRUE), aes(x=Half.life.min, color = "all TFs")) + 
  labs(x = "half-life [min]", y = "cumulative gene density", color = "") +
  scale_color_manual(values = colors) + 
  theme(legend.position="none",
        plot.title = element_text(size=13, face="bold", hjust=0.5, lineheight=1.2),
        plot.subtitle = element_text(size=10, face="italic", hjust=0.5),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12)) + 
  xlim(0,xlimit)
cumDens
