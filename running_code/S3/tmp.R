source('~/research/gAnts/code/config.R')
source('~/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('~/research/2022/ANTS/AnTracks/src/generic.R')
source('~/research/2022/ANTS/AnTracks/src/coords.R')
det <- lapply(expL$det, function(i) init(i, refcoords = h,
                                         class = 'coords'))

d <- exp_spreadsheet[exp_spreadsheet$CONDITION == 'DET', c('Date', 'MATI.TARDA')][1:10, ]
d <- paste(format(strptime(d[, 1], '%d/%m/%Y'), format = '%Y%m%d'), d[, 2], sep='')

for(i in seq_along(det)){
     det[[i]]$date <- d[i]
}

s <- compute_segments(det[[1]])
for(i in seq_along(det)){
     det[[i]]$segments <- s
}

# visualization
draw_hexagons(det[[1]])
food <- get_foodPatches(det[[1]])
food <- lapply(food, function(i){
     hull <- chull(i)
     i[hull, ]
})

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
                         y = c(p1N, p2N, p3N)), aes(x, y))+
     # geom_violin(trim = FALSE)+
     geom_violin(draw_quantiles = 0.5, trim = FALSE)+
     # geom_boxplot(width = 0.05)+
     geom_dotplot(binaxis = 'y', stackdir = 'center')

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
semN <- apply(mN, 2, function(i){
     sd(i) / sqrt(length(i))
})
df <- data.frame(x = 1:length(avgN), y = avgN, sem = semN)[1:21600, ]
df$y <- as.numeric(smooth(df$y, twiceit = T))

ggplot(data = df)+
     geom_errorbar(aes(x = x, ymin = y - sem, ymax = y + sem), alpha = 0.1, color = 'grey80')+
     geom_line(aes(x, y))

avgI <- colMeans(mI)

normDF = cbind(df, i = avgI[1:21600])
normDF$y <- normDF$y / max(normDF$y)
normDF$ic <- cumsum(normDF$i) / max(cumsum(normDF$i))

idx <- vapply(seq(0, 2000, length.out = 5),
              function(i){
                   which.min(abs(cumsum(avgI[1:21600]) - i))
              }, integer(1))

ggplot(data = normDF) + geom_point(aes(x, y))+
     geom_point(color = 'red', aes(x, ic)) +
     scale_y_continuous(limits = c(0, 32 / max(df$y)),
                        breaks = seq(0, 30/max(df$y), length.out = 5),
                        labels = seq(0, 30/max(df$y), length.out = 5)*max(df$y),
                        sec.axis = ~ ic)+
     geom_vline(xintercept = idx)

# sec.axis = sec_axis(trans = ~.,
#                     labels = round(seq(0, 30/max(df$y)+0.25, 0.25)*max(cumsum(normDF$i)))))
                        

ggplot(data = cbind(df, i = avgI[1:21600])) + geom_point(aes(x, y))+
     geom_point(color = 'red', aes(x = x, y = cumsum(i) * max(y) / (max(cumsum(i)))))+
     scale_y_continuous(sec.axis = sec_axis(~.))+
     



