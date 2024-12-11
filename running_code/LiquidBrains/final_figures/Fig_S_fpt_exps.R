source('~/research/gits/AnTracks/src/Experiment.R')

load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/nf.RData')


det_times <- unlist(lapply(det, function(i){
        f <- min(rbindlist(i@food)[['t']])
        f - setDT(i@data)[1, Frame]
}))
det_times <- c(range(det_times), mean(det_times))



R_exp <- seq(100, 1000, 100)

det_result <- rbindlist(lapply(det, function(p){
        data <- setDT(p@data)
        mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
        dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
        indices <- vapply(R_exp, function(i){
                which.max(dists > i)
        }, numeric(1))
        times <- mdata[indices, Frame]
        data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'det')
}))

nf_result <- rbindlist(lapply(nf[-c(1, 2)], function(p){
        data <- setDT(p@data)
        mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
        dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
        indices <- vapply(R_exp, function(i){
                which.max(dists > i)
        }, numeric(1))
        times <- mdata[indices, Frame]
        data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'nf')
}))



panel_a <- ggplot(data = rbindlist(list(det = det_result, nf = nf_result)),
                  aes(factor(R/10), t / 120, fill = factor(rho)))+
        geom_boxplot(outlier.shape = NA, show.legend = FALSE, alpha = 0.6, size = 1)+
        geom_dotplot(stackdir = 'center', binaxis = 'y', dotsize = 125, 
                     binwidth = 1/100, fill = 'lightgrey', show.legend = FALSE)+
        # geom_jitter(width = 0.15, alpha = 0.4, size = 3, show.legend = FALSE)+
        facet_wrap(~ factor(rho, labels = c('Experimental (Food)', 'Control (No-Food)')))+
        scale_y_continuous('First-passage time (min)', breaks = seq(0, 80, 10))+
        scale_x_discrete('Radius (cm)', breaks = R_exp / 10)+
        scale_fill_manual('', values = c('mediumpurple', 'gold3'))+ 
        theme(plot.title = element_text(size = 22))+
        theme(aspect.ratio = 0.75)
