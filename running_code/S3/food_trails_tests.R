source('~/research/2022/ANTS/AnTracks/src/generic.R')
source('~/research/2022/ANTS/AnTracks/src/nodes.R')
source('~/research/2022/ANTS/AnTracks/src/coords.R')

load('/home/polfer/research/2022/ANTS/AnTracks/results/det_coords.RData')
load('/home/polfer/research/2022/ANTS/AnTracks/results/sto_coords.RData')


# activity_heatmap.coords <- function(obj, t){
#         m <- coords2matrix(obj)
#         m <- m[t, ]
#         
#         z <- numeric(length(obj$segments$o))
#         s <- obj$segments[, c('o', 'd')]
#         
#         existing_segments <- s[s[, 'o'] %in% colnames(m) & s[, 'd'] %in% colnames(m), ]
# 
#         indices <- c()  
#         for(i in 1:nrow(existing_segments)){
#                 
#                 idx <- which(s$o == existing_segments[i, 1] & s$d == existing_segments[i, 2])
#                 for(i in idx){
#                         hex$nod
#                 }
#                 
#                 z[idx] <- infotheo::mutinformation(m[, existing_segments[i, 1]], m[, existing_segments[i, 2]])
#                 indices <- c(indices, idx)               
#         }
#         z
# }


# 
# food_trails.coords <- function(obj, method, prob){
# 
#      m <- coords2matrix(obj)
#      m[m == -1] <- 0
#      m[obj$intervals[1]:obj$intervals[2], ]
#      
# }
# 
# food_trails.nodes <- function(obj, method, prob){
#         
#      m <- obj$data
#      m[obj$intervals[1]:obj$intervals[2], ]
#         
# }


for(i in seq_along(det)){
        det[[i]]$trails <- food_trails(det[[i]], method = 'post-exploration', prob = 0.9) 
        det[[i]]$AiT <- activity_in_trail(det[[i]])
        sto[[i]]$trails <- food_trails(sto[[i]], method = 'post-exploration', prob = 0.9) 
        sto[[i]]$AiT <- activity_in_trail(sto[[i]])
}

# for(i in seq_along(sto)){
#         sto[[i]]$trails <- food_trails(sto[[i]], method = 'mid-recruitment', prob = 0.9) 
#         sto[[i]]$AiT <- activity_in_trail(sto[[i]], plot = F)
# }
# 
for(i in seq_along(det)){
        det[[i]]$AiT <- activity_in_trail(det[[i]], plot = F)
        sto[[i]]$AiT <- activity_in_trail(sto[[i]], plot = F)
}
# 
for(i in seq_along(det)){
        det[[i]]$AiT <- det[[i]]$AiT$data
        sto[[i]]$AiT <- sto[[i]]$AiT$data

}

for(i in seq_along(det)){
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































xlist <- lapply(det, food_trails)
ylist <- lapply(det, food_trails, method = 'in-recruitment')
# ylist <- lapply(xlist, function(i) transform_nodes(as.integer(names(i))))


plotL <- lapply(seq_along(xlist), function(i){
     draw_hexagons(det[[i]], add = draw_FoodPatches(det[[i]]) +
                        geom_point(data = data.frame(x = hex$x[xlist[[i]]],
                                                     y = hex$y[xlist[[i]]]),
                                   aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})

ggpubr::ggarrange(plotlist =plotL)

 plotY <- lapply(seq_along(ylist), function(i){
     draw_hexagons(det[[i]], add = draw_FoodPatches(det[[i]]) +
                        geom_point(data = data.frame(x = hex$x[ylist[[i]]],
                                                     y = hex$y[ylist[[i]]]),
                                   aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})

ggpubr::ggarrange(plotlist =plotY)





stolist <- lapply(sto, food_trails)
stoYlist <- lapply(sto, food_trails, method = 'in-recruitment')
# ylist <- lapply(xlist, function(i) transform_nodes(as.integer(names(i))))


plotSTO <- lapply(seq_along(stolist), function(i){
     draw_hexagons(sto[[i]], add = draw_FoodPatches(sto[[i]]) +
                        geom_point(data = data.frame(x = hex$x[stolist[[i]]],
                                                     y = hex$y[stolist[[i]]]),
                                   aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})

ggpubr::ggarrange(plotlist =plotSTO)




## TEST TONTO
activity_in_trail <- function(obj, ...){
        UseMethod('activity_in_trail')
}



nteractions_in_trail <- function(obj){
        UseMethod('interactions_in_trail')
}

ACTIVITY_DET <- list()
for(exp_idx in seq(10)){
        
        lerele <- coords2matrix(det[[exp_idx]])
        lerele[lerele == -1] <- 0
        df_test <- data.frame(t = as.integer(rownames(lerele)), 
                              N = rowSums(lerele),
                              inTrail = rowSums(lerele[, as.character(xlist[[exp_idx]])]),
                              outTrail = rowSums(lerele[, !colnames(lerele) %in% as.character(xlist[[exp_idx]])]))
        mlt_df_test <- reshape2::melt(df_test, id.vars = 't')
        ACTIVITY_DET[[exp_idx]] <- ggplot(data = mlt_df_test, aes(t/120, value, color = variable))+
                geom_line() + ylab('') +
                scale_color_viridis_d('', labels = c('Total activity', 'Activity in trail', 'Activity out trail'))+
                scale_x_continuous('Time (min)', breaks = seq(0, 180, 15))+
                geom_polygon(data = data.frame(x = range(do.call('rbind', det[[exp_idx]]$food)$t)[c(1:2, 2:1)]/120,
                                               y = c(0, 0, max(df_test$N), max(df_test$N))), aes(x, y),
                             fill = muted('green'), color = 'green', alpha = 0.15, linetype = 2, size = 1)+
                guides(color = guide_legend(override.aes = list(size = 2)))
}

png('~/Desktop/activity_trails_determinist.png', 4096, 4096/2, res = 2*100)
ggpubr::ggarrange(plotlist = ACTIVITY_DET, common.legend = T)
dev.off()


png('~/Desktop/activity_trails_stochastic.png', 4096, 4096/2, res = 2*100)
ggpubr::ggarrange(plotlist = ACTIVITY, common.legend = T)
dev.off()
        

# save(sto, file = '~/research/2022/ANTS/AnTracks/results/sto_coords.RData')

plotSTOy <- lapply(seq_along(stoYlist), function(i){
     draw_hexagons(sto[[i]], add = draw_FoodPatches(sto[[i]]) +
                        geom_point(data = data.frame(x = hex$x[stoYlist[[i]]],
                                                     y = hex$y[stoYlist[[i]]]),
                                   aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})

ggpubr::ggarrange(plotlist =plotSTOy)




draw_hexagons(det[[i]], add = draw_FoodPatches(det[[i]]) +
                   geom_point(data = data.frame(x = hex$x[stolist[[i]]],
                                                y = hex$y[stolist[[i]]]),
                              aes(x, y), fill = muted('blue'), size = 5, shape = 21))

draw_hexagons(sto[[10]], add = draw_FoodPatches(sto[[10]])+
     geom_point(data = cbind.data.frame(hex[colnames(doublecheck), c('x', 'y')],
                                        z = colSums(doublecheck))[colSums(doublecheck) > 500, ],
                aes(x, y, fill = z), shape = 21, size = 5))





# ml <- rbind(range(det[[1]]$food$GP1), range(det[[1]]$food$GP2))
# ml2 <- range(do.call('rbind', det[[1]]$food)$t)
# 
# d1 <- det[[1]]
# d1$refcoords <- hex
# # d1$data <- d1$data[d1$data$Frame[ml2[1]:ml2[2]], ]
# # d1$refcoords <- hex
# d1 <- closest_node(d1)
# d1_m <- coords2matrix(d1)
# d1_m[d1_m == -1] <- 0
# 
# sbst1 <- d1_m[ml2[1]:ml2[2], ]
# sbst2 <- d1_m[as.integer(mean(ml2)):(ml2[2] + 1200), ]
# sbst3 <- d1_m[ml2[2]:nrow(d1_m), ]
# 
# 
# q <- 0.85
# 
# x1 <- colSums(sbst1)
# x1[x1 < quantile(x1, prob = q)] <- 0
# x2 <- colSums(sbst2)
# x2[x2 < quantile(x2, prob = q)] <- 0
# x3 <- colSums(sbst3)
# x3[x3 < quantile(x3, prob = q)] <- 0
# 
# 
# 
# draw_hexagons(det[[1]], add = draw_FoodPatches(det[[1]]) +
#                    geom_point(data = data.frame(x = hex$x[as.integer(colnames(d1_m))],
#                                                 y = hex$y[as.integer(colnames(d1_m))],
#                                                 color = x1)[x1 > 0, ],
#                               aes(x, y, fill = color), size = 5, shape = 21))
# 
# draw_hexagons(det[[1]], add = draw_FoodPatches(det[[1]]) +
#                    geom_point(data = data.frame(x = hex$x[as.integer(colnames(d1_m))],
#                                                 y = hex$y[as.integer(colnames(d1_m))],
#                                                 color = x2)[x2 > 0, ],
#                               aes(x, y, fill = color), size = 5, shape = 21))
# 
# draw_hexagons(det[[1]], add = draw_FoodPatches(det[[1]]) +
#                    geom_point(data = data.frame(x = hex$x[as.integer(colnames(d1_m))],
#                                                 y = hex$y[as.integer(colnames(d1_m))],
#                                                 color = x3)[x3 > 0, ],
#                               aes(x, y, fill = color), size = 5, shape = 21))



test75 <- lapply(det,food_trails, prob = .75)
plot75 <- lapply(seq_along(test75), function(i){
        draw_hexagons(det[[i]], add = draw_FoodPatches(det[[i]]) +
                              geom_point(data = data.frame(x = hex$x[test75[[i]]],
                                                           y = hex$y[test75[[i]]]),
                                         aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})
ggpubr::ggarrange(plotlist = plot75)


test90 <- lapply(det,food_trails, prob = .90)
plot90 <- lapply(seq_along(test90), function(i){
        draw_hexagons(det[[i]], add = draw_FoodPatches(det[[i]]) +
                              geom_point(data = data.frame(x = hex$x[test90[[i]]],
                                                           y = hex$y[test90[[i]]]),
                                         aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})
ggpubr::ggarrange(plotlist = plot90)





stotest75 <- lapply(sto,food_trails, prob = .75)
plot75 <- lapply(seq_along(stotest75), function(i){
        draw_hexagons(sto[[i]], add = draw_FoodPatches(sto[[i]]) +
                              geom_point(data = data.frame(x = hex$x[stotest75[[i]]],
                                                           y = hex$y[stotest75[[i]]]),
                                         aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})
ggpubr::ggarrange(plotlist = plot75)

stotest75. <- lapply(sto,food_trails, method = 'in-recruitment', prob = .75)
plot75 <- lapply(seq_along(stotest75.), function(i){
        draw_hexagons(sto[[i]], add = draw_FoodPatches(sto[[i]]) +
                              geom_point(data = data.frame(x = hex$x[stotest75.[[i]]],
                                                           y = hex$y[stotest75.[[i]]]),
                                         aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})
ggpubr::ggarrange(plotlist = plot75)


stotest90 <- lapply(sto,food_trails, prob = .90)
plot90 <- lapply(seq_along(stotest90), function(i){
        draw_hexagons(sto[[i]], add = draw_FoodPatches(sto[[i]]) +
                              geom_point(data = data.frame(x = hex$x[stotest90[[i]]],
                                                           y = hex$y[stotest90[[i]]]),
                                         aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})
ggpubr::ggarrange(plotlist = plot90)





stotest90. <- lapply(sto,food_trails, method = 'post-exploration', prob = .9)
plot90 <- lapply(seq_along(stotest90.), function(i){
        draw_hexagons(sto[[i]], add = draw_FoodPatches(sto[[i]]) +
                              geom_point(data = data.frame(x = hex$x[stotest90.[[i]]],
                                                           y = hex$y[stotest90.[[i]]]),
                                         aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})
ggpubr::ggarrange(plotlist = plot90)


hmap_sto4 <- coords2matrix(sto[[4]])
hmap_sto4[hmap_sto4 == -1] = 0
z <- colSums(hmap_sto4)
z[z < 350] <- 0
z[z > 5000] <- 5000
z <- log(z)
xy <- do.call('rbind', lapply(as.integer(colnames(hmap_sto4)), function(i) hex[hex$old_node == i, c('x','y')]))
xy$z <- z
draw_hexagons(sto[[4]], add = draw_FoodPatches(sto[[4]], add = ggplot(data = xy)+
                      geom_point(aes(x, y, fill = z), shape = 21, size = 5)))
 


dettest90. <- lapply(det,food_trails, method = 'post-exploration', prob = .9)
plot90 <- lapply(seq_along(dettest90.), function(i){
        draw_hexagons(det[[i]], add = draw_FoodPatches(det[[i]]) +
                              geom_point(data = data.frame(x = hex$x[dettest90.[[i]]],
                                                           y = hex$y[dettest90.[[i]]]),
                                         aes(x, y), fill = muted('blue'), size = 5, shape = 21))
})
ggpubr::ggarrange(plotlist = plot90)


DET_ACT_IN_TRAIL <- activity_in_trail(det[[1]], plot = F)
