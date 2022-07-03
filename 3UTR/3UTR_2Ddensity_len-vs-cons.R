## load packages
require("tidyverse")

## load data
hg38_TF <- read.csv("dataR/hg38_pivot_TF.csv")

## prepare data
hg38_TF_TF_ConsDens<- hg38_TF$threeUTR_conservation %>% subset(hg38_TF$isTF) %>% density
hg38_TF_TF_ConsDens_maxInd <- hg38_TF_TF_ConsDens$x[hg38_TF_TF_ConsDens$y %>% which.max()]
hg38_TF_nonTF_ConsDens<- hg38_TF$threeUTR_conservation %>% subset(!hg38_TF$isTF) %>% density
hg38_TF_nonTF_ConsDens_maxInd <- hg38_TF_nonTF_ConsDens$x[hg38_TF_nonTF_ConsDens$y %>% which.max()]
hg38_TF_TF_LenDens<- hg38_TF$threeUTR_length %>% subset(hg38_TF$isTF) %>% density
hg38_TF_TF_LenDens_maxInd <- hg38_TF_TF_LenDens$x[hg38_TF_TF_LenDens$y %>% which.max()]
hg38_TF_nonTF_LenDens<- hg38_TF$threeUTR_length %>% subset(!hg38_TF$isTF) %>% density
hg38_TF_nonTF_LenDens_maxInd <- hg38_TF_nonTF_LenDens$x[hg38_TF_nonTF_LenDens$y %>% which.max()]

## 2D density: length vs. conservation
gAllDensity <- ggplot(hg38_TF, aes(x = threeUTR_conservation,  y = threeUTR_length, color = isTF)) +
  geom_density_2d(alpha = .75, size = 1) + 
  scale_color_manual(name="", 
                     values = c("#F8766D", "#00BFC4"),
                     labels = c("non-TFs", "TFs")) +
  geom_segment(aes(x = hg38_TF_nonTF_ConsDens_maxInd,
                   y = hg38_TF_nonTF_LenDens_maxInd,
                   xend = hg38_TF_TF_ConsDens_maxInd,
                   yend = hg38_TF_TF_LenDens_maxInd),
               arrow=arrow(length=unit(2, "mm")), size = 0.8, color="black") +
  labs(x = "3'-UTR conservation",y = "3'-UTR length [bp]") +
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
gAllDensity
