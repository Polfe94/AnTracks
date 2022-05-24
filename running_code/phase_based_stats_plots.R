exp_condition <- c('det')
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

det_N <- lapply(det, get_N)
det_I <- lapply(det, get_I)

# Population and interactions
for(i in seq_along(det)){
        det[[i]]$N <- get_N(det[[i]])
        det[[i]]$I <- get_I(det[[i]])
}

mN <- numeric(0)
mI <- numeric(0)
for(i in seq_along(det)){
        N <- det[[i]]$N
        I <- det[[i]]$I
        
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


## BOXPLOTS per phase
## +++ N +++
p1N <- vapply(seq_along(det), function(i){
        mean(N[1:(det_phase[i, 1])])
}, numeric(1))
p2N <- vapply(seq_along(det), function(i){
        mean(N[(det_phase[i, 1] + 1):det_phase[i, 2]])
}, numeric(1))
p3N <- vapply(seq_along(det), function(i){
        mean(N[(det_phase[i, 2] + 1):length(N)])
}, numeric(1))

# plot
ggplot(data = data.frame(x = c(rep('p1', 10), rep('p2',10), rep('p3', 10)),
                         y = c(p1N, p2N, p3N),
                         z = rep(seq(0, 50, length.out = 10), 3)), aes(x, y))+
        # geom_violin(trim = FALSE)+
        geom_violin(draw_quantiles = 0.5, trim = FALSE)+
        # geom_boxplot(width = 0.05)+
        geom_dotplot(binaxis = 'y', stackdir = 'center', aes(x,y, fill = z))+
        scale_fill_viridis()

## +++ I +++
p1I <- vapply(seq_along(det), function(i){
        mean(I[1:(det_phase[i, 1])])
}, numeric(1))
p2I <- vapply(seq_along(det), function(i){
        mean(I[(det_phase[i, 1] + 1):det_phase[i, 2]])
}, numeric(1))
p3I <- vapply(seq_along(det), function(i){
        mean(I[(det_phase[i, 2] + 1):length(I)])
}, numeric(1))

# plot
ggplot(data = data.frame(x = c(rep('p1', 10), rep('p2',10), rep('p3', 10)),
                         y = c(p1I, p2I, p3I)), aes(x, y))+
        geom_violin(trim = FALSE)+
        geom_boxplot(width = 0.05)+
        geom_dotplot(binaxis = 'y', stackdir = 'center')


avgN <- colMeans(mN)
# semN <- apply(mN, 2, function(i){
#         sd(i) / sqrt(length(i))
# })
# df <- data.frame(x = 1:length(avgN), y = avgN, sem = semN)[1:21600, ]
# df$y <- as.numeric(smooth(df$y, twiceit = T))
# 
# ggplot(data = df)+
#         geom_errorbar(aes(x = x, ymin = y - sem, ymax = y + sem), alpha = 0.1, color = 'grey80')+
#         geom_line(aes(x, y))

avgI <- colMeans(mI)
 
df <- data.frame(x = 1:21600 /120, y = cumsum(avgI[1:21600])/ max(cumsum(avgI[1:21600])),
           z = as.numeric(smooth(avgN[1:21600]), twiceit = T))

tiff('G:/research/2022/AnTracks/plots/det_.tiff', width = 2400, height = 1200, res = 180)
ggplot(data = df) + 
        # geom_bracket(xmin = 0, xmax = mean(det_phase[, 1])/120,
        #              y.position = 0.94, color = 'black',
        #              label = 'Exploration', label.size = 5, size = 0.6, tip.length = 0.02) +
        # geom_bracket(xmin = mean(det_phase[, 1])/120, xmax = mean(det_phase[, 2])/120,
        #              y.position = 1, color = 'black',
        #              label = 'Food collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        # geom_bracket(xmin = mean(det_phase[, 2])/120, xmax = 180,
        #              y.position = 1.05, color = 'black',
        #              label = 'Post-collection', label.size = 5, size = 0.6, tip.length = 0.02) +
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
        food <- get_foodPatches(sto[[1]])
        sto[[i]]$food <- food[1:2]
}

draw_FoodPatches(sto[[1]], add = draw_hexagons(sto[[1]]))



sto_phase <- t(vapply(sto, function(i){
        x <- food_detection(i)
        c(min(x, na.rm = T), max(x, na.rm = T))
}, numeric(2)))

# Population and interactions
for(i in seq_along(sto)){
        sto[[i]]$N <- get_N(sto[[i]])
        sto[[i]]$I <- get_I(sto[[i]])
}

mN <- numeric(0)
mI <- numeric(0)
for(i in seq_along(sto)){
        N <- sto[[i]]$N
        I <- sto[[i]]$I
        
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


## BOXPLOTS per phase
## +++ N +++
mN <- mN[-2, ]
mI <- mI[-2, ]

p1N <- vapply(seq_along(sto), function(i){
        mean(N[1:(sto_phase[i, 1])])
}, numeric(1))
p2N <- vapply(seq_along(sto), function(i){
        mean(N[(sto_phase[i, 1] + 1):sto_phase[i, 2]])
}, numeric(1))
p3N <- vapply(seq_along(sto), function(i){
        mean(N[(sto_phase[i, 2] + 1):length(N)])
}, numeric(1))

# plot
ggplot(data = data.frame(x = c(rep('p1', 9), rep('p2',9), rep('p3', 9)),
                         y = c(p1N, p2N, p3N),
                         z = rep(seq(0, 40, length.out = 9), 3)), aes(x, y))+
        # geom_violin(trim = FALSE)+
        geom_violin(draw_quantiles = 0.5, trim = FALSE)+
        # geom_boxplot(width = 0.05)+
        geom_dotplot(binaxis = 'y', stackdir = 'center', aes(x,y))+
        scale_fill_viridis()

## +++ I +++


avgN <- colMeans(mN)
avgI <- colMeans(mI)

df <- data.frame(x = 1:21600 /120, y = cumsum(avgI[1:21600])/ max(cumsum(avgI[1:21600])),
                 z = as.numeric(smooth(avgN[1:21600]), twiceit = T))

tiff('G:/research/2022/AnTracks/plots/sto_.tiff', width = 2400, height = 1200, res = 180)
ggplot(data = df) + 
        # geom_bracket(xmin = 0, xmax = mean(sto_phase[, 1])/120,
        #              y.position = 0.94, color = 'black',
        #              label = 'Exploration', label.size = 5, size = 0.6, tip.length = 0.02) +
        # geom_bracket(xmin = mean(sto_phase[, 1])/120, xmax = mean(sto_phase[, 2])/120,
        #              y.position = 1, color = 'black',
        #              label = 'Food collection', label.size = 5, size = 0.6, tip.length = 0.02) +
        # geom_bracket(xmin = mean(sto_phase[, 2])/120, xmax = 180,
        #              y.position = 1.05, color = 'black',
        #              label = 'Post-collection', label.size = 5, size = 0.6, tip.length = 0.02) +
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
