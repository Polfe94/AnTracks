source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
load('~/research/2022/ANTS/AnTracks/data/det.RData')
load('~/research/2022/ANTS/AnTracks/data/sto.RData')

#### +++ OPTIMAL PATHS +++ ####
## DETERMINIST
det <- lapply(det, optimal_food_trails)
det <- lapply(det, interaction_matrix)
det <- lapply(det, ait, TRUE)
det <- lapply(det, iit, TRUE)

det_optimal_trails <- lapply(det, plot_trails)
det_optimal_ait <- lapply(det, plot_AiT)
det_optimal_iit <- lapply(det, plot_IiT, norm = T) # normalized version

ggarrange(plotlist = det_optimal_trails, labels = sapply(det, function(i) i@date))
ggarrange(plotlist = det_optimal_ait, labels = sapply(det, function(i) i@date), common.legend = T, hjust = -3)
ggarrange(plotlist = det_optimal_iit, labels = sapply(det, function(i) i@date), common.legend = T, hjust = -3)

## STOCHASTIC
sto <- lapply(sto, optimal_food_trails)
sto <- lapply(sto, interaction_matrix)
sto <- lapply(sto, ait, TRUE)
sto <- lapply(sto, iit, TRUE)

sto_optimal_trails <- lapply(sto, plot_trails)
sto_optimal_ait <- lapply(sto, plot_AiT)
sto_optimal_iit <- lapply(sto, plot_IiT, norm = T) # normalized version

ggarrange(plotlist = sto_optimal_trails, labels = sapply(sto, function(i) i@date))
ggarrange(plotlist = sto_optimal_ait, labels = sapply(sto, function(i) i@date), common.legend = T, hjust = -3)
ggarrange(plotlist = sto_optimal_iit, labels = sapply(sto, function(i) i@date), common.legend = T, hjust = -3)


#### +++ FOOD TRAILS +++ ####
## DETERMINIST
det <- lapply(det, food_trails)
det <- lapply(det, interaction_matrix)
det <- lapply(det, function(i){
     i@AiT <- data.frame()
     i <- ait(i, TRUE)
     i
})
det <- lapply(det, function(i){
     i@IiT <- data.frame()
     i <- iit(i, TRUE)
     i
})

det_food_trails <- lapply(det, plot_trails)
det_ait <- lapply(det, plot_AiT)
det_iit <- lapply(det, plot_IiT, norm = T) # normalized version

ggarrange(plotlist = det_food_trails, labels = sapply(det, function(i) i@date))
ggarrange(plotlist = det_ait, labels = sapply(det, function(i) i@date), common.legend = T, hjust = -3)
ggarrange(plotlist = det_iit, labels = sapply(det, function(i) i@date), common.legend = T, hjust = -3)

## STOCHASTIC
sto <- lapply(sto, food_trails)
sto <- lapply(sto, interaction_matrix)
sto <- lapply(sto, function(i){
     i@AiT <- data.frame()
     i <- ait(i, TRUE)
     i
})
sto <- lapply(sto, function(i){
     i@IiT <- data.frame()
     i <- iit(i, TRUE)
     i
})

sto_food_trails <- lapply(sto, plot_trails)
sto_ait <- lapply(sto, plot_AiT)
sto_iit <- lapply(sto, plot_IiT, norm = T) # normalized version

ggarrange(plotlist = sto_food_trails, labels = sapply(sto, function(i) i@date))
ggarrange(plotlist = sto_ait, labels = sapply(sto, function(i) i@date), common.legend = T, hjust = -3)
ggarrange(plotlist = sto_iit, labels = sapply(sto, function(i) i@date), common.legend = T, hjust = -3)
