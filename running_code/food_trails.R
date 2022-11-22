source('~/research/2022/ANTS/AnTracks/src/generic.R')
source('~/research/2022/ANTS/AnTracks/src/nodes.R')
source('~/research/2022/ANTS/AnTracks/src/coords.R')

load('/home/polfer/research/2022/ANTS/AnTracks/results/det_coords.RData')
load('/home/polfer/research/2022/ANTS/AnTracks/results/sto_coords.RData')

for(i in seq_along(det)){
     det[[i]]$trails <- food_trails(det[[i]], method = 'post-exploration', prob = 0.9) 
     det[[i]]$AiT <- activity_in_trail(det[[i]])
     sto[[i]]$trails <- food_trails(sto[[i]], method = 'post-exploration', prob = 0.9) 
     sto[[i]]$AiT <- activity_in_trail(sto[[i]])
}

for(i in seq_along(det)){
     det[[i]]$AiT <- activity_in_trail(det[[i]], plot = F)
     sto[[i]]$AiT <- activity_in_trail(sto[[i]], plot = F)
     
     det[[i]]$IiT <- interactions_in_trail(det[[i]])
     sto[[i]]$IiT <- interactions_in_trail(sto[[i]])
}


det_trails <- lapply(det, plot_trails)
det_AiT <- lapply(det, plot_AiT)
det_diff <- lapply(det, function(i){
     y <- i$AiT$value[i$AiT$variable == 'inTrail'] - i$AiT$value[i$AiT$variable == 'outTrail']
     data <- data.frame(x = seq_along(y), y = y)
     ggplot(data = data, aes(x / 120, y)) + geom_line()+
          scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
          scale_y_continuous('Difference in activity', breaks = seq(min(y), max(y), 5))+
          geom_vline(xintercept = range(do.call('rbind', i$food)$t)/120, linetype = 2, size = 1)
})
det_diffXcap <- lapply(det, function(i){
     y <- i$AiT$value[i$AiT$variable == 'inTrail']
     data <- data.frame(x = i$AiT$value[i$AiT$variable == 'N'], y = y)
     ggplot(data = data, aes(x, y)) + geom_point()+
          scale_x_continuous('Total activity', breaks = seq(0, max(data$x), 5))+
          scale_y_continuous('Difference in activity', breaks = seq(min(y), max(y), 5))
})
det_Idiff <- lapply(det, function(i){
     y <- i$IiT$value[i$IiT$variable == 'inTrail'] - i$IiT$value[i$IiT$variable == 'outTrail']
     data <- data.frame(x = seq_along(y), y = y)
     ggplot(data = data, aes(x / 120, y)) + geom_line()+
          scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
          scale_y_continuous('Difference in activity', breaks = seq(min(y), max(y), 1.5))+
          geom_vline(xintercept = range(do.call('rbind', i$food)$t)/120, linetype = 2, size = 1)
})
det_IiT <- lapply(det, plot_IiT)
det_spatial <- lapply(det, plot_interactions)

sto_trails <- lapply(sto, plot_trails)
sto_AiT <- lapply(sto, plot_AiT)
sto_diff <- lapply(sto, function(i){
     y <- i$AiT$value[i$AiT$variable == 'inTrail'] - i$AiT$value[i$AiT$variable == 'outTrail']
     data <- data.frame(x = seq_along(y), y = y)
     ggplot(data = data, aes(x / 120, y)) + geom_line()+
          scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
          scale_y_continuous('Difference in activity', breaks = seq(min(y), max(y), 5))+
          geom_vline(xintercept = range(do.call('rbind', i$food)$t)/120, linetype = 2, size = 1)
})
sto_diffXcap <- lapply(sto, function(i){
     y <- i$AiT$value[i$AiT$variable == 'inTrail']
     data <- data.frame(x = i$AiT$value[i$AiT$variable == 'N'], y = y)
     ggplot(data = data, aes(x, y)) + geom_point()+
          scale_x_continuous('Total activity', breaks = seq(0, max(data$x), 5))+
          scale_y_continuous('Difference in activity', breaks = seq(min(y), max(y), 5))
})
sto_Idiff <- lapply(sto, function(i){
     y <- i$IiT$value[i$IiT$variable == 'inTrail'] - i$IiT$value[i$IiT$variable == 'outTrail']
     data <- data.frame(x = seq_along(y), y = y)
     ggplot(data = data, aes(x / 120, y)) + geom_line()+
          scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
          scale_y_continuous('Difference in activity', breaks = seq(min(y), max(y), 1.5))+
          geom_vline(xintercept = range(do.call('rbind', i$food)$t)/120, linetype = 2, size = 1)
})
sto_IiT <- lapply(sto, plot_IiT)
sto_spatial <- lapply(sto, plot_interactions)

ggpubr::ggarrange(plotlist = det_trails, common.legend = T)
ggpubr::ggarrange(plotlist = det_AiT, common.legend = T, labels = gsub('json', '', names(det)), 
                  font.label = list(size = 15, face = 'plain'), hjust = -0.75)
ggpubr::ggarrange(plotlist = det_diff, common.legend = T)
ggpubr::ggarrange(plotlist = det_diffXcap, common.legend = T)
ggpubr::ggarrange(plotlist = det_Idiff, common.legend = T)
ggpubr::ggarrange(plotlist = det_IiT, common.legend = T, labels = gsub('json', '', names(det)), 
                  font.label = list(size = 15, face = 'plain'), hjust = -0.75)
ggpubr::ggarrange(plotlist = det_spatial)

ggpubr::ggarrange(plotlist = sto_trails, common.legend = T)
ggpubr::ggarrange(plotlist = sto_AiT, common.legend = T, labels = gsub('json', '', names(sto)), 
                  font.label = list(size = 15, face = 'plain'), hjust = -0.75)
ggpubr::ggarrange(plotlist = sto_diff, common.legend = T)
ggpubr::ggarrange(plotlist = sto_diffXcap, common.legend = T)
ggpubr::ggarrange(plotlist = sto_Idiff, common.legend = T)
ggpubr::ggarrange(plotlist = sto_IiT, common.legend = T, labels = gsub('json', '', names(det)), 
                  font.label = list(size = 15, face = 'plain'), hjust = -0.75)
ggpubr::ggarrange(plotlist = sto_spatial)