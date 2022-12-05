source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
load('~/research/2022/ANTS/AnTracks/data/det.RData')
load('~/research/2022/ANTS/AnTracks/data/sto.RData')

det <- lapply(det, food_trails)
# det <- lapply(det, function(i){
#      i@interactions <- data.frame()
#      i <- interaction_matrix(i)
#      i
# })
det <- lapply(det, interaction_matrix)

DET_trails <- lapply(det, plot_trails)
ggarrange(plotlist = L)


det <- lapply(det, ait)
det <- lapply(det, iit)


detAiT <- lapply(det, plot_AiT)
detIiT <- lapply(det, plot_IiT)

ggarrange(plotlist = detAiT, labels = vapply(det, function(i) i@date, character(1)), common.legend = T,
          hjust = -3)

ggarrange(plotlist = detIiT, labels = vapply(det, function(i) i@date, character(1)), common.legend = T,
          hjust = -3)





### EL MATEIX PER ESTOCÀSTICS... ###
load('~/research/2022/ANTS/AnTracks/data/sto.RData')

sto <- lapply(sto, food_trails)
sto <- lapply(sto, interaction_matrix)

STO_trails <- lapply(sto, plot_trails)
ggarrange(plotlist = STO_trails)


sto <- lapply(sto, ait)
sto <- lapply(sto, iit)


stoAiT <- lapply(sto, plot_AiT)
stoIiT <- lapply(sto, plot_IiT)

ggarrange(plotlist = stoAiT, labels = vapply(sto, function(i) i@date, character(1)), common.legend = T,
          hjust = -3)

ggarrange(plotlist = stoIiT, labels = vapply(sto, function(i) i@date, character(1)), common.legend = T,
          hjust = -3)
