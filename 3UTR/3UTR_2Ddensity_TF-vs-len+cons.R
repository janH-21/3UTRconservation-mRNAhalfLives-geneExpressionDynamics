## load packages
require(tidyverse)
require(ggpubr)

## load data
hg38_halfLife <- read.csv("dataR/hg38_pivot_half.csv")

## prepare data
hg38_halfLife_TF_ConsDens<- hg38_halfLife$threeUTR_conservation %>% subset(hg38_halfLife$isTF) %>% density
hg38_halfLife_TF_ConsDens_maxInd <- hg38_halfLife_TF_ConsDens$x[hg38_halfLife_TF_ConsDens$y %>% which.max()]
hg38_halfLife_TF_LenDens<- hg38_halfLife$threeUTR_length %>% subset(hg38_halfLife$isTF) %>% density
hg38_halfLife_TF_LenDens_maxInd <- hg38_halfLife_TF_LenDens$x[hg38_halfLife_TF_LenDens$y %>% which.max()]
hg38_halfLife_TF_halfDens<- hg38_halfLife$Half.life.min %>% subset(hg38_halfLife$isTF) %>% density
hg38_halfLife_TF_halfDens_maxInd <- hg38_halfLife_TF_halfDens$x[hg38_halfLife_TF_halfDens$y %>% which.max()]
hg38_halfLife_nonTF_ConsDens<- hg38_halfLife$threeUTR_conservation %>% subset(!hg38_halfLife$isTF) %>% density
hg38_halfLife_nonTF_ConsDens_maxInd <- hg38_halfLife_nonTF_ConsDens$x[hg38_halfLife_nonTF_ConsDens$y %>% which.max()]
hg38_halfLife_nonTF_LenDens<- hg38_halfLife$threeUTR_length %>% subset(!hg38_halfLife$isTF) %>% density
hg38_halfLife_nonTF_LenDens_maxInd <- hg38_halfLife_nonTF_LenDens$x[hg38_halfLife_nonTF_LenDens$y %>% which.max()]
hg38_halfLife_nonTF_halfDens<- hg38_halfLife$Half.life.min %>% subset(!hg38_halfLife$isTF) %>% density
hg38_halfLife_nonTF_halfDens_maxInd <- hg38_halfLife_nonTF_halfDens$x[hg38_halfLife_nonTF_halfDens$y %>% which.max()]

## 2D densty: half-life vs. conservation
gConsHalf <- ggplot(hg38_halfLife, aes(x = threeUTR_conservation,  y = Half.life.min , color = isTF)) +
  geom_density_2d(alpha = .75, size = 1)  + 
  scale_color_manual(name="", 
                     values = c("#F8766D", "#00BFC4"),
                     labels = c("non-TFs", "TFs"))+
  geom_segment(aes(x = hg38_halfLife_nonTF_ConsDens_maxInd,
                   y = hg38_halfLife_nonTF_halfDens_maxInd,
                   xend = hg38_halfLife_TF_ConsDens_maxInd,
                   yend = hg38_halfLife_TF_halfDens_maxInd),
               arrow=arrow(length=unit(2, "mm")), size = 0.8, color="black") +
  labs(x = "3'-UTR conservation",y = "Half-life [min]") +
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        legend.text = element_text(size=12),
        plot.title = element_text(size=8, face="bold", hjust=0.5, lineheight=1.2),
        plot.subtitle = element_text(size=6, face="italic", hjust=0.5),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12)) 

## 2D density: half-life vs. length
gLenHalf <- ggplot(hg38_halfLife, aes(x = threeUTR_length,  y = Half.life.min , color = isTF)) +
  geom_density_2d(alpha = .75, size = 1)  + 
  scale_color_manual(name="", 
                     values = c("#F8766D", "#00BFC4"),
                     labels = c("non-TFs", "TFs"))+
  geom_segment(aes(x = hg38_halfLife_nonTF_LenDens_maxInd,
                   y = hg38_halfLife_nonTF_halfDens_maxInd,
                   xend = hg38_halfLife_TF_LenDens_maxInd,
                   yend = hg38_halfLife_TF_halfDens_maxInd),
               arrow=arrow(length=unit(2, "mm")), size = 0.8, color="black") +
  labs(x = "3'-UTR length [bp]", y = "Half-life [min]") +
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        legend.text = element_text(size=12),
        plot.title = element_text(size=8, face="bold", hjust=0.5, lineheight=1.2),
        plot.subtitle = element_text(size=6, face="italic", hjust=0.5),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12)) 

## plot
ggarrange(gConsHalf, gLenHalf, nrow = 1, common.legend = TRUE)
