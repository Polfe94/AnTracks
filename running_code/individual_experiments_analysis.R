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

circleFun <- function(center = c(0,0),diameter = 1, npoints = 50){
        r = diameter / 2
        tt <- seq(0,2*pi,length.out = npoints)
        xx <- center[1] + r * cos(tt)
        yy <- center[2] + r * sin(tt)
        return(data.frame(x = xx, y = yy))
}


## BOLDNESS EXPERIMENTS ## 
setwd('C:/Users/POL/Desktop/CSV_Ind_experiments_test/Boldness/')

files <- list.files()
df <- list()
for(i in seq_along(files)){
        tmp <- read.csv(files[i])
        idx <- which(colnames(tmp) %in% c('fx', 'fy'))
        tmp$ID <- i
        if(length(idx)){
                df[[i]] <- tmp[, -idx]
        } else {
                df[[i]] <- tmp
        }
        
}
df <- do.call('rbind', df)

## speed by id
ggplot(data = df, aes(factor(ID), speed)) + geom_boxplot()

## cover distance by id
ggplot(data = df, aes(factor(ID), cdist)) + geom_boxplot()

## wall distance by id
ggplot(data = df, aes(factor(ID), wdist)) + geom_boxplot()

## time spent in mask by id
ggplot(data = df, aes(time, inmask)) + geom_line() + facet_wrap(~ factor(ID))

## displacement by id
ggplot(data = df, aes(time, cumdispl)) + geom_line() + facet_wrap(~ factor(ID))

# coordinates by id
ggplot(data = df, aes(cx, cy)) + geom_point(aes(color = time)) + facet_wrap(~ factor(ID), scales = 'free') +
 geom_polygon(aes(mx, my, fill = factor(ID))) + scale_fill_viridis_d()

## heatmap by id
ggarrange(plotlist = lapply(seq_len(max(df$ID)), function(i){
        ggplot(data = df[df$ID == i, ], aes(cx, cy)) + geom_density_2d_filled(bins = 30, show.legend = F) +
                scale_fill_manual(values = c('grey80', mixOmics::color.jet(30)))
}))


## SOCIABILITY EXPERIMENTS ##


setwd('C:/Users/POL/Desktop/CSV_Ind_experiments_test/Sociability/')

files <- list.files()
df <- list()
for(i in seq_along(files)){
        tmp <- read.csv(files[i])
        idx <- which(colnames(tmp) %in% c('fx', 'fy'))
        tmp$ID <- i
        if(length(idx)){
                df[[i]] <- tmp[, -idx]
        } else {
                df[[i]] <- tmp
        }
        
}
df <- do.call('rbind', df)

## speed by id
ggplot(data = df, aes(factor(ID), speed)) + geom_boxplot(outlier.shape = NA)+
        scale_y_continuous(limits = c(0, 7.5))

## cover distance by id
ggplot(data = df, aes(factor(ID), cdist)) + geom_boxplot()

## wall distance by id
ggplot(data = df, aes(factor(ID), wdist)) + geom_boxplot()

## time spent in mask by id
ggplot(data = df, aes(time, inmask)) + geom_line() + facet_wrap(~ factor(ID))

## displacement by id
ggplot(data = df, aes(time, cumdispl)) + geom_line() + facet_wrap(~ factor(ID))

## cordinates by id
ggplot(data = df, aes(cx, cy)) + geom_point(aes(color = time)) + facet_wrap(~ factor(ID))# +
        # geom_polygon(aes(mx, my, fill = factor(ID))) + scale_fill_viridis_d()

## with mask:
# ggarrange(plotlist = lapply(seq_len(max(df$ID)), function(i){
#         c <- circleFun(center = c(max(df[df$ID == i, 'mx'], na.rm = T) - (max(df[df$ID == i, 'mx'], na.rm = T) - min(df[df$ID == i, 'mx'], na.rm = T))/2,
#                                   max(df[df$ID == i, 'my'], na.rm = T) - (max(df[df$ID == i, 'my'], na.rm = T) - min(df[df$ID == i, 'my'], na.rm = T))/2),
#                        diameter = mean((max(df[df$ID == i, 'mx'], na.rm = T) - min(df[df$ID == i, 'mx'], na.rm = T)), 
#                                        (max(df[df$ID == i, 'my'], na.rm = T) - min(df[df$ID == i, 'my'], na.rm = T))),
#                        npoints = 50)
#         
#         ggplot(data = df[df$ID == i, ], aes(cx, cy)) + geom_point(aes(color = time)) + 
#                 geom_polygon(data = c, aes(x, y))
# }), common.legend = T)



## heatmap by id
# ggplot(data = df, aes(cx, cy)) + geom_density_2d_filled(bins = 30, show.legend = F) + facet_wrap(~ factor(ID), scales = 'free')+
#         scale_fill_manual(values = c('grey80', mixOmics::color.jet(30)))
ggarrange(plotlist = lapply(seq_len(max(df$ID)), function(i){
        ggplot(data = df[df$ID == i, ], aes(cx_c, cy_c)) + geom_density_2d_filled(bins = 30, show.legend = F) +
                scale_fill_manual(values = c('grey80', mixOmics::color.jet(30)))
}))

