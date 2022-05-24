exp_condition <- c('')
visual <- F
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


m <- coords2matrix(det[[1]])

########## DETS +++ ##########

zlist_p1 <- vector('list', length(det))
zlist_p2 <- vector('list', length(det))
zlist_p3 <- vector('list', length(det))
f <- round(apply(do.call('rbind', lapply(seq_along(det),
                             function(i) range(do.call('rbind', det[[i]]$food)$t))),
                   2, mean))

# times javi 0-986 secs, 986-2696 secs, 2696-end
for(i in seq_along(det)){
     t <- 1:f[1]
     zlist_p1[[i]] <- pairwise_cov(det[[i]], t)
     t <- f[1]:f[2]
     zlist_p2[[i]] <- pairwise_cov(det[[i]], t)
     t <- f[2]:max(det[[i]]$data$Frame)
     zlist_p3[[i]] <- pairwise_cov(det[[i]], t)
}

Zp1 <- apply(do.call('rbind', zlist_p1), 2, mean)
Zp1[Zp1 < 0] <- norm_range(Zp1[Zp1 < 0], a = -1, b = -.Machine$double.eps)
Zp1[Zp1 > 0] <- norm_range(Zp1[Zp1 > 0], a = .Machine$double.eps, b = 1)
Zp2 <- apply(do.call('rbind', zlist_p2), 2, mean)
Zp2[Zp2 < 0] <- norm_range(Zp2[Zp2 < 0], a = -1, b = -.Machine$double.eps)
Zp2[Zp2 > 0] <- norm_range(Zp2[Zp2 > 0], a = .Machine$double.eps, b = 1)
Zp3 <- apply(do.call('rbind', zlist_p3), 2, mean)
Zp3[Zp3 < 0] <- norm_range(Zp3[Zp3 < 0], a = -1, b = -.Machine$double.eps)
Zp3[Zp3 > 0] <- norm_range(Zp3[Zp3 > 0], a = .Machine$double.eps, b = 1)

# draw_hexagons(det[[1]], z = z, add = draw_hexagons(det[[1]], size = 2.3, color = 'black'), size = 2)
library(gridExtra)
dp1 <- draw_hexagons(det[[1]], z = Zp1, add = draw_hexagons(det[[1]], size = 3, color = 'black'),
                     size = 2.7, show.legend = F) + theme(aspect.ratio = 0.5)
dp2 <- draw_hexagons(det[[1]], z = Zp2, add = draw_hexagons(det[[1]], size = 3, color = 'black'),
                     size = 2.7, show.legend = F) + theme(aspect.ratio = 0.5)
dp3 <- draw_hexagons(det[[1]], z = Zp3, add = draw_hexagons(det[[1]], size = 3, color = 'black'),
                     size = 2.7, show.legend = F) + theme(aspect.ratio = 0.5)
grid.arrange(dp1,dp2,dp3,nrow = 3, ncol = 1)



# z[z < 0.01 | z > -0.01] <- -max(z)

########## STOS +++ ##########
# sto <- sto[-2]
for(i in seq_along(sto)){
     t <- food_detection(sto[[i]])
     sto[[i]]$food$GP1$t <- t[1:6]
     sto[[i]]$food$GP2$t <- t[7:12]
}
Szlist_p1 <- vector('list', length(sto))
Szlist_p2 <- vector('list', length(sto))
Szlist_p3 <- vector('list', length(sto))
Sf <- round(apply(do.call('rbind', lapply(seq_along(sto),
                                         function(i) range(do.call('rbind', sto[[i]]$food)$t, na.rm = T))),
                 2, mean))

# times javi 0-986 secs, 986-2696 secs, 2696-end
for(i in seq_along(sto)){
     t <- 1:Sf[1]
     Szlist_p1[[i]] <- pairwise_cov(sto[[i]], t)
     t <- Sf[1]:Sf[2]
     Szlist_p2[[i]] <- pairwise_cov(sto[[i]], t)
     t <- Sf[2]:max(sto[[i]]$data$Frame)
     Szlist_p3[[i]] <- pairwise_cov(sto[[i]], t)
}

sZp1 <- apply(do.call('rbind', Szlist_p1), 2, mean)
# sZp1[sZp1 == 0] <- NA
sZp1[sZp1 < 0] <- norm_range(sZp1[sZp1 < 0], a = -1, b = -.Machine$double.eps)
sZp1[sZp1 > 0] <- norm_range(sZp1[sZp1 > 0], a = .Machine$double.eps, b = 1)
sZp2 <- apply(do.call('rbind', Szlist_p2), 2, mean)
sZp2[sZp2 < 0] <- norm_range(sZp2[sZp2 < 0], a = -1, b = -.Machine$double.eps)
sZp2[sZp2 > 0] <- norm_range(sZp2[sZp2 > 0], a = .Machine$double.eps, b = 1)
sZp3 <- apply(do.call('rbind', Szlist_p3), 2, mean)
sZp3[sZp3 < 0] <- norm_range(sZp3[sZp3 < 0], a = -1, b = -.Machine$double.eps)
sZp3[sZp3 > 0] <- norm_range(sZp3[sZp3 > 0], a = .Machine$double.eps, b = 1)

# draw_hexagons(sto[[1]], z = z, add = draw_hexagons(sto[[1]], size = 2.3, color = 'black'), size = 2)
library(gridExtra)
Sdp1 <- draw_hexagons(sto[[1]], z = sZp1, add = draw_hexagons(sto[[1]], size = 3, color = 'black'),
                     size = 2.7, show.legend = F) + theme(aspect.ratio = 0.5)
Sdp2 <- draw_hexagons(sto[[1]], z = sZp2, add = draw_hexagons(sto[[1]], size = 3, color = 'black'),
                     size = 2.7, show.legend = F) + theme(aspect.ratio = 0.5)
Sdp3 <- draw_hexagons(sto[[1]], z = sZp3, add = draw_hexagons(sto[[1]], size = 3, color = 'black'),
                     size = 2.7, show.legend = F) + theme(aspect.ratio = 0.5)
grid.arrange(Sdp1,Sdp2,Sdp3,nrow = 3, ncol = 1)



####### BOTH PLOTS TOGETHER ###########

grid.arrange(dp1, Sdp1, dp2, Sdp2, dp3, Sdp3, nrow = 3, ncol = 2)
