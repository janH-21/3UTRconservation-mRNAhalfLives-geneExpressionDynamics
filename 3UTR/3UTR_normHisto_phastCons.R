## load packages
require(tidyverse)

## load data
hg38_ordered <- read.csv("dataR/hg38_ordered.csv")

## normalized histogram
hg38_hist_transcript <- hg38_ordered %>% 
  dplyr::filter(sequence_type %in% c("CDS", "3UTR", "5UTR")) %>%
  dplyr::filter(length > 10) %>%
  ggplot(aes(x=conservation, color=sequence_type)) +
  geom_histogram(aes(y = ..density..), position="identity", alpha = 0, binwidth = 0.02, size = 0.33) +
  scale_color_manual(name="",
                     labels = c("3'-UTRs", "5'-UTRs", "CDS"),
                     values = c("red", "green", "blue"))+
  labs(title = "Histogram normalized to group count\n", x = "Conservation",y = "Transcript Density") +
  theme(legend.position="bottom", 
        legend.box = "horizontal",
        legend.text = element_text(size=12),
        plot.title = element_text(size=8, face="bold", hjust=0.5, lineheight=1.2),
        plot.subtitle = element_text(size=6, face="italic", hjust=0.5),
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12))  
hg38_hist_transcript
