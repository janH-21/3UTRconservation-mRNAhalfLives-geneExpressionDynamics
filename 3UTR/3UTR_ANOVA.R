## load packages
require(tidyverse)
require(rstatix)
require(magrittr)

## load data
hg38_pivot_half <- read.csv("dataR/hg38_pivot_half.csv")

## build model
anovaData <- cbind(halfLife = hg38_pivot_half$Half.life.min,
                   isTF = (hg38_pivot_half$isTF == TRUE),
                   isCO = (hg38_pivot_half$threeUTR_conservation >0.6),
                   isLE = (hg38_pivot_half$threeUTR_length > 1000)
) %>% as.data.frame()

ANOVAmodel1 <- (lm((halfLife) ~ isTF+isCO+isLE, data=anovaData)) # without interactions
ANOVAmodel2 <- (lm((halfLife) ~ isTF*isCO*isLE, data=anovaData)) # with interactions

##  statistics
anovaData %>% 
  group_by(isTF, isCO, isLE) %>%
  get_summary_stats(halfLife, type = "median_iqr")

summary(ANOVAmodel1)
summary(ANOVAmodel2) 

## plot results
anova2 <- summary(ANOVAmodel2) 
barANOVA <- as.data.frame(anova2$coefficients) %>% 
  tibble::rownames_to_column('group') %>%
  set_colnames(c("group", "Estimate", "Error", "tValue", "p")) %>%
  ggplot(aes(x=factor(group, c("(Intercept)","isCO","isLE","isTF","isCO:isLE","isTF:isCO","isTF:isLE","isTF:isCO:isLE")))) + 
  geom_bar(aes(y = Estimate), position="dodge",stat="identity") +
  geom_errorbar(aes(ymin = Estimate - Error, ymax = Estimate + Error), width = 0.5, size = 1.2, color = "black") +
  ylab("ANOVA Model Estimate: Influence on half-life [min]") +
  xlab("Group") + 
  theme(axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12),  
        axis.text.x = element_text(size=12),  
        axis.text.y = element_text(size=12))
barANOVA

## AIC
akaike <- AIC(ANOVAmodel1, ANOVAmodel2) 
exp((akaike$AIC[2]-akaike$AIC[1])/2)
