exp_condition <- c('')
rm(list = ls()[grepl('json', ls())], colony, current_dir, exp_condition, exps, visual)
source('G:/research/MoveAnts/code/config.R')
source('C:/Users/POL/Desktop/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')
load('G:/research/2022/AnTracks/results/det_coords.RData')
load('G:/research/2022/AnTracks/results/sto_coords.RData')

#### ++++ DETERMINISTIC ++++ ####
d <- exp_spreadsheet[exp_spreadsheet$CONDITION == 'DET', c('Date', 'MATI.TARDA')][1:10, ]
d <- paste(format(strptime(d[, 1], '%d/%m/%Y'), format = '%Y%m%d'), d[, 2], sep='')

for(i in seq_along(det)){
        det[[i]]$date <- d[i]
}

food <- get_foodPatches(det[[1]])

for(i in seq_along(det)){
        det[[i]]$food <- food[1:2]
}
draw_FoodPatches(det[[1]], add = draw_hexagons(det[[1]]))

for(i in seq_along(det)){
        x <- food_detection(det[[i]])
        det[[i]]$food[[1]]$t <- x[1:6]
        det[[i]]$food[[2]]$t <- x[7:12]
}

det_phase <- t(vapply(det, function(i){
        x <- food_detection(i)
        c(min(x, na.rm = T), max(x, na.rm = T))
}, numeric(2)))

# Mean cov
mcov_det <- vapply(1:10, function(i) mean(det[[i]]$local_cov), numeric(1))
mcov_t <- lapply(1:10, function(i) seq_along(det[[i]]$local_cov)* 5)
mcov_phase_det <- lapply(1:3, function(i){
        vapply(seq_along(det), function(x){
                if(i == 1){
                        idx <- which(mcov_t[[x]] <= det_phase[x, 1])
                        
                } else if (i == 2){
                        idx <- which(mcov_t[[x]] > det_phase[x, 1] & mcov_t[[x]] <= det_phase[x, 2])
                } else {
                        idx <- which(mcov_t[[x]] > det_phase[x, 2])
                }
                mean(det[[x]]$local_cov[idx])
        }, numeric(1))
})

# Population and interactions
for(i in seq_along(det)){
        det[[i]]$N <- get_N(det[[i]])
        det[[i]]$I <- get_I(det[[i]])
}


p1N <- numeric()
p2N <- numeric()
p3N <- numeric()
p1I <- numeric()
p2I <- numeric()
p3I <- numeric()
mN <- numeric()
mI <- numeric()
for(i in seq_along(det)){
        N <- det[[i]]$N
        I <- det[[i]]$I
        
        # N
        p1N <- c(p1N, mean(N[1:(det_phase[i, 1])]))
        p2N <- c(p2N, mean(N[(det_phase[i, 1] + 1):det_phase[i, 2]]))
        p3N <- c(p3N, mean(N[(det_phase[i, 2] + 1):length(N)]))
        
        # I
        p1I <- c(p1I, mean(I[1:(det_phase[i, 1])]))
        p2I <- c(p2I, mean(I[(det_phase[i, 1] + 1):det_phase[i, 2]]))
        p3I <- c(p3I, mean(I[(det_phase[i, 2] + 1):length(I)]))
        
        if(i > 1){
                if(length(N) < ncol(mN)){
                        N <- c(N, numeric(ncol(mN) - length(N)))
                        I <- c(I, numeric(ncol(mI) - length(I)))
                } else if (length(N) > ncol(mN)){
                        m <- matrix(0, nrow = nrow(mN), ncol = length(N) - ncol(mN))
                        mN <- cbind(mN, m)
                        mI <- cbind(mI, m)
                } else {
                        F
                }
        }
        
        mN <- rbind(mN, N)
        mI <- rbind(mI, I)
}

det_phases <- list(N = list(p1N = p1N, p2N = p2N, p3N = p3N),
                   I = list(p1I = p1I, p2I = p2I, p3I = p3I))

## BOXPLOTS per phase
# +++ N +++

ggplot(data = data.frame(x = c(rep('p1', 10), rep('p2',10), rep('p3', 10)),
                         y = c(p1N, p2N, p3N),
                         z = rep(seq(0, 50, length.out = 10), 3)), aes(x, y))+
        # geom_violin(trim = FALSE)+
        geom_violin(draw_quantiles = 0.5, trim = FALSE)+
        # geom_boxplot(width = 0.05)+
        geom_dotplot(binaxis = 'y', stackdir = 'center', aes(x,y, fill = z))+
        scale_fill_viridis()

## +++ I +++

# plot
ggplot(data = data.frame(x = c(rep('p1', 10), rep('p2',10), rep('p3', 10)),
                         y = c(p1I, p2I, p3I)), aes(x, y))+
        geom_violin(trim = FALSE)+
        geom_boxplot(width = 0.05)+
        geom_dotplot(binaxis = 'y', stackdir = 'center')


avgN_det <- colMeans(mN)
# semN <- apply(mN, 2, function(i){
#         sd(i) / sqrt(length(i))
# })
# df <- data.frame(x = 1:length(avgN), y = avgN, sem = semN)[1:21600, ]
# df$y <- as.numeric(smooth(df$y, twiceit = T))
# 
# ggplot(data = df)+
#         geom_errorbar(aes(x = x, ymin = y - sem, ymax = y + sem), alpha = 0.1, color = 'grey80')+
#         geom_line(aes(x, y))

avgI_det <- colMeans(mI)

df <- data.frame(x = 1:21600 /120, y = cumsum(avgI_det[1:21600])/ max(cumsum(avgI_det[1:21600])),
                 z = as.numeric(smooth(avgN_det[1:21600]), twiceit = T))

tiff('G:/research/2022/AnTracks/plots/det_.tiff', width = 2400, height = 1200, res = 180)
ggplot(data = df) + 
        geom_bracket(xmin = 0, xmax = mean(det_phase[, 1])/120,
                     y.position = 0.94, color = 'black',
                     label = 'Exploration', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(det_phase[, 1])/120, xmax = mean(det_phase[, 2])/120,
                     y.position = 1, color = 'black',
                     label = 'Food collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(det_phase[, 2])/120, xmax = 180,
                     y.position = 1.05, color = 'black',
                     label = 'Post-collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_line(aes(x, z / max(z), color = 'N'), size = 1.2)+
        geom_line(aes(x, y, color = 'IC'), size = 2)+
        
        scale_y_continuous('Cumulated Interactions (IC)', limits = c(0, 30/max(avgN[1:21600])), 
                           breaks = seq(0,30/max(avgN[1:21600]), length.out = 9),
                           labels = seq(0, 2000, 250),
                           sec.axis = sec_axis(~ ., 
                                               breaks = seq(0, 30/max(avgN[1:21600]), length.out = 7),
                                               labels = seq(0, 30, 5), name = 'Population (N)'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20))+
        scale_color_manual('', values = c('N' = viridis(1, begin= 0.25, end = 0.25),
                                          'IC' = viridis(1, begin= 0.7, end = 0.7)))+
        guides(color = guide_legend(override.aes = list(size = 3)))+
        theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
              axis.title.y.right =  element_text(margin = margin(0,0, 0, 10)))
dev.off()


#### ++++ STOCHASTIC ++++ ####

d <- exp_spreadsheet[exp_spreadsheet$CONDITION == 'STOCH', c('Date', 'MATI.TARDA')][1:10, ]
d <- paste(format(strptime(d[, 1], '%d/%m/%Y'), format = '%Y%m%d'), d[, 2], sep='')

for(i in seq_along(sto)){
        sto[[i]]$date <- d[i]
}

for(i in seq_along(sto)){
        food <- get_foodPatches(sto[[i]])
        sto[[i]]$food <- food[1:2]
}

draw_FoodPatches(sto[[1]], add = draw_hexagons(sto[[1]]))



sto_phase <- t(vapply(sto, function(i){
        x <- food_detection(i)
        c(min(x, na.rm = T), max(x, na.rm = T))
}, numeric(2)))

# Mean cov
mcov_sto <- vapply(1:10, function(i) mean(sto[[i]]$local_cov), numeric(1))
mcov_sto <- mcov_sto[-2]
mcov_t <- lapply(seq_along(sto), function(i) seq_along(sto[[i]]$local_cov)* 5)
mcov_phase_sto <- lapply(1:3, function(i){
        vapply(c(1, 3:10), function(x){
                if(i == 1){
                        idx <- which(mcov_t[[x]] <= sto_phase[x, 1])
                        
                } else if (i == 2){
                        idx <- which(mcov_t[[x]] > sto_phase[x, 1] & mcov_t[[x]] <= sto_phase[x, 2])
                } else {
                        idx <- which(mcov_t[[x]] > sto_phase[x, 2])
                }
                mean(sto[[x]]$local_cov[idx])
        }, numeric(1))
})

# Population and interactions
for(i in seq_along(sto)){
        sto[[i]]$N <- get_N(sto[[i]])
        sto[[i]]$I <- get_I(sto[[i]])
}

p1N <- numeric()
p2N <- numeric()
p3N <- numeric()
p1I <- numeric()
p2I <- numeric()
p3I <- numeric()
mN <- numeric()
mI <- numeric()
for(i in seq_along(sto)){
        N <- sto[[i]]$N
        I <- sto[[i]]$I
        
        # N
        p1N <- c(p1N, mean(N[1:(sto_phase[i, 1])]))
        p2N <- c(p2N, mean(N[(sto_phase[i, 1] + 1):sto_phase[i, 2]]))
        p3N <- c(p3N, mean(N[(sto_phase[i, 2] + 1):length(N)]))
        
        # I
        p1I <- c(p1I, mean(I[1:(sto_phase[i, 1])]))
        p2I <- c(p2I, mean(I[(sto_phase[i, 1] + 1):sto_phase[i, 2]]))
        p3I <- c(p3I, mean(I[(sto_phase[i, 2] + 1):length(I)]))
        
        if(i > 1){
                if(length(N) < ncol(mN)){
                        N <- c(N, numeric(ncol(mN) - length(N)))
                        I <- c(I, numeric(ncol(mI) - length(I)))
                } else if (length(N) > ncol(mN)){
                        m <- matrix(0, nrow = nrow(mN), ncol = length(N) - ncol(mN))
                        mN <- cbind(mN, m)
                        mI <- cbind(mI, m)
                } else {
                        F
                }
        }
        
        mN <- rbind(mN, N)
        mI <- rbind(mI, I)
}

sto_phases <- list(N = list(p1N = p1N, p2N = p2N, p3N = p3N),
                   I = list(p1I = p1I, p2I = p2I, p3I = p3I))

## BOXPLOTS per phase
## +++ N +++
ggplot(data = data.frame(x = c(rep('p1', 9), rep('p2',9), rep('p3', 9)),
                         y = c(p1N[-2], p2N[-2], p3N[-2]),
                         z = rep(seq(0, 40, length.out = 9), 3)), aes(x, y))+
        geom_boxplot()+
        geom_dotplot(binaxis = 'y', stackdir = 'center', aes(x,y, fill = z))+
        scale_fill_viridis()

## +++ I +++
ggplot(data = data.frame(x = c(rep('p1', 9), rep('p2',9), rep('p3', 9)),
                         y = c(p1I[-2], p2I[-2], p3I[-2]),
                         z = rep(seq(0, 40, length.out = 9), 3)), aes(x, y))+
        geom_boxplot()+
        geom_dotplot(binaxis = 'y', stackdir = 'center', aes(x,y, fill = z))+
        scale_fill_viridis()



avgN_sto <- colMeans(mN[-2, ])
avgI_sto <- colMeans(mI[-2, ])

df <- data.frame(x = 1:21600 /120, y = cumsum(avgI[1:21600])/ max(cumsum(avgI[1:21600])),
                 z = as.numeric(smooth(avgN[1:21600]), twiceit = T))

tiff('G:/research/2022/AnTracks/plots/sto_.tiff', width = 2400, height = 1200, res = 180)
ggplot(data = df) + 
        geom_bracket(xmin = 0, xmax = mean(sto_phase[, 1])/120,
                     y.position = 0.94, color = 'black',
                     label = 'Exploration', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(sto_phase[, 1])/120, xmax = mean(sto_phase[, 2])/120,
                     y.position = 1, color = 'black',
                     label = 'Food collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(sto_phase[, 2])/120, xmax = 180,
                     y.position = 1.05, color = 'black',
                     label = 'Post-collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_line(aes(x, z / max(z), color = 'N'), size = 1.2)+
        geom_line(aes(x, y, color = 'IC'), size = 2)+
        
        scale_y_continuous('Cumulated Interactions (IC)', limits = c(0, 37/max(avgN[1:21600])), 
                           breaks = seq(0,35/max(avgN[1:21600]), length.out = 9),
                           labels = seq(0, 1000, 125),
                           sec.axis = sec_axis(~ ., 
                                               breaks = seq(0, 35/max(avgN[1:21600]), length.out = 8),
                                               labels = seq(0, 35, 5), name = 'Population (N)'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20))+
        scale_color_manual('', values = c('N' = viridis(1, begin= 0.25, end = 0.25),
                                          'IC' = viridis(1, begin= 0.7, end = 0.7)))+
        guides(color = guide_legend(override.aes = list(size = 3)))+
        theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
              axis.title.y.right =  element_text(margin = margin(0,0, 0, 10)))
dev.off()


####### DET AND STO PLOTS COMBINED #######
library(ggnewscale)

## N
a <- ggplot(data = data.frame(x = substr(names(c(do.call('c', det_phases$N),do.call('c', sto_phases$N))), 2,2),
                         y = c(do.call('c', det_phases$N),do.call('c', sto_phases$N)),
                         z = c(rep('Determinist', 30), rep('Stochastic', 27))), aes(x, y))+

        geom_boxplot(aes(fill = factor(z)), width = 0.9, show.legend = F)+ 
        scale_fill_manual('', values = c('white', 'grey70'))+
        new_scale('fill')+ 
        scale_x_discrete('', breaks = c(1, 2, 3), labels = c('Exploration', 'Food collection', 'Post-collection')) + 
        scale_y_continuous('Population size', breaks = seq(0, 30, length.out = 7))+
        geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.05, position = 'dodge',
                     dotsize = 7.5, aes(fill = factor(z)))+
        scale_fill_manual('Experiment condition', values = c('black', 'black'))+
        guides(fill = guide_legend(override.aes = list(fill = c('white', 'grey70'))))

## I
b <- ggplot(data = data.frame(x = substr(names(c(do.call('c', det_phases$I),do.call('c', sto_phases$I))), 2,2),
                         y = c(do.call('c', det_phases$I),do.call('c', sto_phases$I)),
                         z = c(rep('Determinist', 30), rep('Stochastic', 27))), aes(x, y))+
        geom_boxplot(aes(fill = factor(z)), width = 0.9, show.legend = F)+ 
        scale_fill_manual('', values = c('white', 'grey70'))+
        new_scale('fill')+
        scale_y_continuous('Interactions', breaks = seq(0, 0.25, length.out = 6), limits = c(0, 0.25)) +
        scale_x_discrete('', breaks = c(1, 2, 3), labels = c('Exploration', 'Food collection', 'Post-collection')) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = 0.005, position = 'dodge',
                     dotsize = 0.75, aes(fill = factor(z)))+
        scale_fill_manual('Experiment condition', values = c('black', 'black'))+
        guides(fill = guide_legend(override.aes = list(fill = c('white', 'grey70'))))

c <- ggplot(data = data.frame(x = 1:21600 /120, y = norm_range(cumsum(avgI_det[1:21600]), a = 0, b = 1),
                              z = norm_range(as.numeric(smooth(avgN_det[1:21600]), twiceit = T), a = 0, b = 1))) + 
        geom_bracket(xmin = 0, xmax = mean(det_phase[, 1])/120,
                     y.position = 0.94, color = 'black',
                     label = 'Exploration', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(det_phase[, 1])/120, xmax = mean(det_phase[, 2])/120,
                     y.position = 1.02, color = 'black',
                     label = 'Food collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(det_phase[, 2])/120, xmax = 180,
                     y.position = 1.1, color = 'black',
                     label = 'Post-collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_line(aes(x, z / max(z), color = 'N'), size = 1.2)+
        geom_line(aes(x, y, color = 'IC'), size = 2)+
        
        scale_y_continuous('Cumulated Interactions (IC)', limits = c(0, 1.15), 
                           breaks = seq(0, 2000/max(cumsum(avgI_det[1:21600])), length.out = 9),
                           labels = seq(0, 2000, 250),
                           sec.axis = sec_axis(~ ., 
                                               breaks = seq(0, 2000/max(cumsum(avgI_det[1:21600])), length.out = 7),
                                               labels = seq(0, 30, 5), name = 'Population (N)'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20))+
        scale_color_manual('', values = c('N' = viridis(1, begin= 0.25, end = 0.25),
                                          'IC' = viridis(1, begin= 0.7, end = 0.7)))+
        guides(color = guide_legend(override.aes = list(size = 3)))+
        theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
              axis.title.y.right =  element_text(margin = margin(0,0, 0, 10)))


d <- ggplot(data = data.frame(x = 1:21600 /120, y = norm_range(cumsum(avgI_sto[1:21600]), a = 0, b = 1),
                              z = norm_range(as.numeric(smooth(avgN_sto[1:21600]), twiceit = T), a = 0, b = 1))) + 
        geom_bracket(xmin = 0, xmax = mean(sto_phase[, 1])/120,
                     y.position = 0.94, color = 'black',
                     label = 'Exploration', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(sto_phase[, 1])/120, xmax = mean(sto_phase[, 2])/120,
                     y.position = 1.02, color = 'black',
                     label = 'Food collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_bracket(xmin = mean(sto_phase[, 2])/120, xmax = 180,
                     y.position = 1.1, color = 'black',
                     label = 'Post-collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        geom_line(aes(x, z / max(z), color = 'N'), size = 1.2)+
        geom_line(aes(x, y, color = 'IC'), size = 2)+
        
        scale_y_continuous('Cumulated Interactions (IC)', limits = c(0, 1.15), 
                           breaks = seq(0, 1000/ max(cumsum(avgI_sto)), length.out = 9),
                           labels = seq(0, 1000, 125),
                           sec.axis = sec_axis(~ ., 
                                               breaks = seq(0, 1, length.out = 8),
                                               labels = seq(0, 35, 5), name = 'Population (N)'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20))+
        scale_color_manual('', values = c('N' = viridis(1, begin= 0.25, end = 0.25),
                                          'IC' = viridis(1, begin= 0.7, end = 0.7)))+
        guides(color = guide_legend(override.aes = list(size = 3)))+
        theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
              axis.title.y.right =  element_text(margin = margin(0,0, 0, 10)))


grid.arrange(a, b, c, d)
# grid.arrange(ggpubr::ggarrange(a, b, common.legend = T, legend = 'right'),
#              ggpubr::ggarrange(c, d, common.legend = T, legend = 'right'))
