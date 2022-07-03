#### load packages ####
require(tidyverse)
require(seewave)
require(ggpubr)
require(scales)
require(reshape2)

#### function definitions ####
## bode plot of phase
plotPhase <- function(dataModel, variable, colorScale, fontSettings, legendTitle, xlabel) {
  ggPhase <- ggplot(data = dataModel, aes(x=frequency, y=mad/3.1451*180)) +
    geom_line(aes(group = factor(variable), color=factor(variable))) +
    scale_x_settings +
    xlab(xlabel) + ylab("Phase shift [°]") +
    colorScale + fontSettings + guides(col=guide_legend(legendTitle))
  ggPhase
}

## bode plot of energy
plotEnergy <- function(dataModel, variable, colorScale, fontSettings, legendTitle, xlabel) {
  ggPhase <- ggplot(data = dataModel, aes(x=frequency, y=20*log10(energy))) +
    geom_line(aes(group = factor(variable), color=factor(variable))) +
    scale_x_settings +
    xlab(xlabel) + ylab("Magnitude [dB]") +
    colorScale + fontSettings + guides(col=guide_legend(legendTitle))
  ggPhase
}

## bode plot of standardized energy
plotStandEnergy <- function(dataModel, variable, colorScale, fontSettings, legendTitle, xlabel) {
  for(iter in unique(variable)){
    dataModel$normEnergy[variable == iter] <- dataModel$energy[variable == iter] / dataModel$energy[variable == iter & dataModel$frequency == min(dataModel$frequency)]
  }
  
  ggPhase <- ggplot(data = dataModel, aes(x=frequency, y=20*log10(normEnergy))) +
    geom_line(aes(group = factor(variable), color=factor(variable))) +
    scale_x_settings +
    xlab(xlabel) + ylab("Magnitude [dB]") +
    colorScale + fontSettings + guides(col=guide_legend(legendTitle))
  ggPhase
}

#### aesthetics ####
colorScale6 <- scale_color_manual(values=c('#000033','#000066','#000099','#0000CC','#0000FF','#3366CC'))
colorScale7 <- scale_color_manual(values=c('#000033','#000066','#000099','#0000CC','#0000FF','#3366CC','#6699CC')) 
colorScale2 <- scale_color_manual(values=c('#000033','#6699CC')) 
colorScale4 <- scale_color_manual(values=c('black','blue','green','red')) 

fontSettings <- theme(legend.text = element_text(size=12), legend.box = "horizontal",
                      plot.title = element_text(size=13, face="bold", hjust=0.5, lineheight=1.2),
                      plot.subtitle = element_text(size=10, face="italic", hjust=0.5),
                      axis.title.x = element_text(size=12), 
                      axis.title.y = element_text(size=12),  
                      axis.text.x = element_text(size=12),  
                      axis.text.y = element_text(size=12))

scale_x_settings <- scale_x_continuous(trans = log10_trans(), 
                                       breaks = trans_breaks("log10", function(x) 10^x),
                                       labels = trans_format("log10", math_format(10^.x))) 

#### bode plots ####
## bode plot: mRNA degradation rate d0
dataModel <- read.table("modelAnalysisR/result_d.txt",header=T) %>% arrange(frequency) 
ggPhase.d0 <- plotPhase(dataModel, dataModel$d0, colorScale6, fontSettings, expression("Degradation rate d"[0]*":"), "Frequency [Hz]")
ggEnergy.d0 <- plotEnergy(dataModel, dataModel$d0, colorScale6, fontSettings,expression("Degradation rate d"[0]*":"), "Frequency [Hz]")
ggStandEnergy.d0 <- plotStandEnergy(dataModel, dataModel$d0, colorScale6, fontSettings, expression("Degradation rate d"[0]*":"), "Frequency [Hz]")

gg.d0 <- ggarrange(ggEnergy.d0, ggStandEnergy.d0, ggPhase.d0, 
                   nrow = 3, ncol = 1, common.legend = TRUE, legend = "bottom")
gg.d0

## bode plot: mRNA translation rate k1
dataModel <- read.table("modelAnalysisR/result_k.txt",header=T) %>% arrange(frequency) 
ggPhase.k1 <- plotPhase(dataModel, dataModel$k1, colorScale6, fontSettings, expression("Translation rate k"[1]*":"), "Frequency [Hz]")
ggEnergy.k1 <- plotEnergy(dataModel, dataModel$k1, colorScale6, fontSettings, expression("Translation rate k"[1]*":"), "Frequency [Hz]")
ggStandEnergy.k1 <- plotStandEnergy(dataModel, dataModel$k1, colorScale6, fontSettings, expression("Translation rate k"[1]*":"), "Frequency [Hz]")

gg.k1 <- ggarrange(ggEnergy.k1, ggStandEnergy.k1, ggPhase.k1, 
                   nrow = 3, ncol = 1, common.legend = TRUE, legend = "bottom")
gg.k1

## bode plot biological: data
dataModel <- read.table("modelAnalysisR/result_bio.txt",header=T) %>% arrange(frequency) 
ggPhase.bio <- plotPhase(dataModel, dataModel$parameterset, colorScale2, fontSettings, "TF/target:", "Frequency [Hz]")
ggEnergy.bio <- plotEnergy(dataModel, dataModel$parameterset, colorScale2, fontSettings, "TF/target:", "Frequency [Hz]")
ggStandEnergy.bio <- plotStandEnergy(dataModel, dataModel$parameterset, colorScale2, fontSettings, "TF/target:", "Frequency [Hz]")

gg.bio <- ggarrange(ggEnergy.bio, ggStandEnergy.bio, ggPhase.bio, 
                    nrow = 3, ncol = 1, common.legend = TRUE, legend = "bottom")
gg.bio


#### time course example ####
dataTimeCourse <- read.table("modelAnalysisR/timecourse_RelA_VCAM1_f_0.5.txt",header=T)

ggStart <- dataTimeCourse[1:100, c('t','k0', 'TF_mRNA', 'TF_protein', 'target_mRNA')] %>% melt(id.vars = 't', mesure.vars = c('k0', 'TF_mRNA', 'TF_protein', 'target_mRNA')) %>%
  ggplot(aes(x=t, y=value, colour = variable)) +
  geom_path() +
  colorScale4 +
  fontSettings +
  ylab("Concentration [AU]") +
  xlab("Time [h]") + 
  coord_cartesian(ylim=c(-1,4))

ggAll <- dataTimeCourse[1:2500, c('t','k0', 'TF_mRNA', 'TF_protein', 'target_mRNA')] %>% melt(id.vars = 't', mesure.vars = c('k0', 'TF_mRNA', 'TF_protein', 'target_mRNA')) %>%
  ggplot(aes(x=t, y=value, colour = variable)) +
  geom_line() +
  colorScale4 +
  fontSettings +
  ylab("Concentration [AU]") +
  xlab("Time [h]")

ggEqui <- dataTimeCourse[6000:6100, c('t','k0', 'TF_mRNA', 'TF_protein', 'target_mRNA')] %>% melt(id.vars = 't', mesure.vars = c('k0', 'TF_mRNA', 'TF_protein', 'target_mRNA')) %>%
  ggplot(aes(x=t, y=value, colour = variable)) +
  geom_line() +
  colorScale4 +
  fontSettings +
  ylab("Concentration [AU]") +
  xlab("Time [h]")

gg.time <- ggarrange(ggStart, ggAll, ggEqui, 
                     nrow = 3, ncol = 1, common.legend = TRUE, legend = "bottom")
gg.time