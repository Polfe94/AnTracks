exp_condition <- c('')

source('~/research/gAnts/code/config.R')
source('~/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('~/research/2022/ANTS/AnTracks/src/generic.R')
source('~/research/2022/ANTS/AnTracks/src/coords.R')

load('/home/polfer/research/2022/ANTS/AnTracks/results/sto_coords.RData')
load('/home/polfer/research/2022/ANTS/AnTracks/results/det_coords.RData')

rm(list = ls()[grepl('json', ls())], visual, colony, current_dir, dirL, exp_condition, exps, conditionL)

########## DETS +++ ##########

DET_zp1 <- vector('list', length(det))
DET_zp2 <- vector('list', length(det))
DET_zp3 <- vector('list', length(det))
det_phase <- t(vapply(det, function(i){
        x <- food_detection(i)
        c(min(x, na.rm = T), max(x, na.rm = T))
}, numeric(2)))

# times javi 0-986 secs, 986-2696 secs, 2696-end
for(i in seq_along(det)){
     t <- 1:det_phase[i, 1]
     DET_zp1[[i]] <- pairwise_cov(det[[i]], t)
     t <- det_phase[i, 1]:det_phase[i, 2]
     DET_zp2[[i]] <- pairwise_cov(det[[i]], t)
     t <- det_phase[i, 2]:max(det[[i]]$data$Frame)
     DET_zp3[[i]] <- pairwise_cov(det[[i]], t)
}

dZp1 <- apply(do.call('rbind', DET_zp1), 2, mean)
dZp1[dZp1 < 0] <- norm_range(dZp1[dZp1 < 0], a = -1, b = -.Machine$double.eps)
dZp1[dZp1 > 0] <- norm_range(dZp1[dZp1 > 0], a = .Machine$double.eps, b = 1)
dZp2 <- apply(do.call('rbind', DET_zp2), 2, mean)
dZp2[dZp2 < 0] <- norm_range(dZp2[dZp2 < 0], a = -1, b = -.Machine$double.eps)
dZp2[dZp2 > 0] <- norm_range(dZp2[dZp2 > 0], a = .Machine$double.eps, b = 1)
dZp3 <- apply(do.call('rbind', DET_zp3), 2, mean)
dZp3[dZp3 < 0] <- norm_range(dZp3[dZp3 < 0], a = -1, b = -.Machine$double.eps)
dZp3[dZp3 > 0] <- norm_range(dZp3[dZp3 > 0], a = .Machine$double.eps, b = 1)

# draw_hexagons(det[[1]], z = z, add = draw_hexagons(det[[1]], size = 2.3, color = 'black'), size = 2)
dp1 <- draw_hexagons(det[[1]], z = dZp1, add = draw_hexagons(det[[1]], size = 1.8, color = 'black',
                                                             add = draw_FoodPatches(det[[1]])),
                     size = 1.6, show.legend = F) + 
        geom_point(data = hex[634, ], aes(x, y), size = 2, color = 'grey40', shape = 17)+
        theme_void() + theme(aspect.ratio = 0.5,
                             axis.title.y = element_text(angle = 90, 
                                                         margin = unit(c(0, 0, 0, 0), 'pt'),
                                                         size = 9))+
        ylab('Exploration') + ggtitle('                    DET')
        
dp2 <- draw_hexagons(det[[1]], z = dZp2, add = draw_hexagons(det[[1]], size = 1.8, color = 'black',
                                                             add = draw_FoodPatches(det[[1]])),
                     size = 1.6, show.legend = F)+ theme_void() + 
        geom_point(data = hex[634, ], aes(x, y), size = 2, color = 'grey40', shape = 17)+
        theme(aspect.ratio = 0.5,
                axis.title.y = element_text(angle = 90, 
                                            margin = unit(c(0, 0, 0, 0), 'pt'),
                                            size = 9))+
        ylab('Recruitment') + ggtitle('')
dp3 <- draw_hexagons(det[[1]], z = dZp3, add = draw_hexagons(det[[1]], size = 1.8, color = 'black',
                                                             add = draw_FoodPatches(det[[1]])),
                     size = 1.6, show.legend = F)+theme_void() + 
        geom_point(data = hex[634, ], aes(x, y), size = 2, color = 'grey40', shape = 17)+
        theme(aspect.ratio = 0.5, axis.title.y = element_text(angle = 90, 
                                                              margin = unit(c(0, 0, 0, 0), 'pt'),
                                                              size = 9))+
        ylab('Post recruitment') + ggtitle('')
grid.arrange(dp1,dp2,dp3,nrow = 3, ncol = 1)

########## STOS +++ ##########
sto_phase <- t(vapply(sto, function(i){
        x <- food_detection(i)
        c(min(x, na.rm = T), max(x, na.rm = T))
}, numeric(2)))


STO_zp1 <- vector('list', length(sto))
STO_zp2 <- vector('list', length(sto))
STO_zp3 <- vector('list', length(sto))

# times javi 0-986 secs, 986-2696 secs, 2696-end
for(i in seq_along(sto)){
     t <- 1:sto_phase[i, 1]
     STO_zp1[[i]] <- pairwise_cov(sto[[i]], t)
     t <- sto_phase[i, 1]:sto_phase[i, 2]
     STO_zp2[[i]] <- pairwise_cov(sto[[i]], t)
     t <- sto_phase[i, 2]:max(sto[[i]]$data$Frame)
     STO_zp3[[i]] <- pairwise_cov(sto[[i]], t)
}

sZp1 <- apply(do.call('rbind', STO_zp1[-2]), 2, mean)
sZp1[sZp1 < 0] <- norm_range(sZp1[sZp1 < 0], a = -1, b = -.Machine$double.eps)
sZp1[sZp1 > 0] <- norm_range(sZp1[sZp1 > 0], a = .Machine$double.eps, b = 1)
sZp2 <- apply(do.call('rbind', STO_zp2[-2]), 2, mean)
sZp2[sZp2 < 0] <- norm_range(sZp2[sZp2 < 0], a = -1, b = -.Machine$double.eps)
sZp2[sZp2 > 0] <- norm_range(sZp2[sZp2 > 0], a = .Machine$double.eps, b = 1)
sZp3 <- apply(do.call('rbind', STO_zp3[-2]), 2, mean)
sZp3[sZp3 < 0] <- norm_range(sZp3[sZp3 < 0], a = -1, b = -.Machine$double.eps)
sZp3[sZp3 > 0] <- norm_range(sZp3[sZp3 > 0], a = .Machine$double.eps, b = 1)

# draw_hexagons(sto[[1]], z = z, add = draw_hexagons(sto[[1]], size = 2.3, color = 'black'), size = 2)

sp1 <- draw_hexagons(sto[[1]], z = sZp1, add = draw_hexagons(sto[[1]], size = 1.8, color = 'black',
                                                             add = draw_FoodPatches(sto[-2])),
                     size = 1.6, show.legend = F) + theme_void()+ 
        geom_point(data = hex[634, ], aes(x, y), size = 2, color = 'grey40', shape = 17)+
        theme(aspect.ratio = 0.5) + ggtitle('                      STO')
        
sp2 <- draw_hexagons(sto[[1]], z = sZp2, add = draw_hexagons(sto[[1]], size = 1.8, color = 'black',
                                                             add = draw_FoodPatches(sto[-2])),
                     size = 1.6, show.legend = F) + theme_void()+ 
        geom_point(data = hex[634, ], aes(x, y), size = 2, color = 'grey40', shape = 17)+
        theme(aspect.ratio = 0.5) + ggtitle('')
sp3 <- draw_hexagons(sto[[1]], z = sZp3, add = draw_hexagons(sto[[1]], size = 1.8, color = 'black',
                                                             add = draw_FoodPatches(sto[-2])),
                     size = 1.6, show.legend = F) + theme_void()+ 
        geom_point(data = hex[634, ], aes(x, y), size = 2, color = 'grey40', shape = 17)+
        theme(aspect.ratio = 0.5) + ggtitle('')
grid.arrange(Sdp1,Sdp2,Sdp3,nrow = 3, ncol = 1)



####### BOTH PLOTS TOGETHER #######

png(filename = '~/research/2022/ANTS/Figs_spinGlasses/spatial_correlations.png', 8000, 8000, res = 1600)
ggarrange(dp1, sp1, dp2, sp2, dp3, sp3, nrow = 3, ncol = 2)
dev.off()
grid.arrange(dp1, sp1, dp2, sp2, dp3, sp3, nrow = 3, ncol = 2)
