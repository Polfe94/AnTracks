source('G:/research/gits/AnTracks/src/Experiment.R')
load('G:/research/gits/AnTracks/data/nf.RData')
load('G:/research/gits/AnTracks/data/det.RData')
source('G:/research/gits/AnTracks/src/fit_functions.R')

library(segmented)

bin_function <- function(v, res = 10, minv = 0){
        cut(v, breaks = seq(minv, ceiling(max(v)/res)*res, res))
}

DET_scouts <- rbindlist(lapply(det, function(p){
        t <- min(rbindlist(p@food)[['t']])
        data <- setDT(p@data)[Frame <= t]
        # data <- setDT(p@data)
        dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
        set(data, j = 'd', value = as.numeric(dmatrix))
        data[, .(maxd = max(d), meand = mean(d)), by = 'N_ind']
}))
DET_scouts[['tag']] <- 'DET'

NFD_scouts <- rbindlist(lapply(nf[-c(1, 2)], function(p){
        data <- setDT(p@data)
        dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
        set(data, j = 'd', value = as.numeric(dmatrix))
        data[, .(maxd = max(d), meand = mean(d)), by = 'N_ind']
}))
NFD_scouts[['tag']] <- 'NFD'


data_binned <- rbind(DET_scouts, NFD_scouts)
data_binned[, bin := bin_function(meand, res = 10, minv = 80), by = c('tag')]
data <- data_binned[, .N, by = c('bin', 'tag')]
fit_data <- data[, .(x = ((as.numeric(bin)+8)*10) / 10,
                     y = N), by = 'tag']
lmodels <- lapply(unique(fit_data[['tag']]),
                  function(m){
                          data <- fit_data[tag == m]
                          mdl <- lm(formula = log(y) ~ x, data = data)
                          mdl$call <- as.call(list(
                                  quote(lm),
                                  formula = quote(log(y) ~ x),
                                  data = substitute(fit_data[tag == m], list(m = m))
                          ))
                          mdl
                  })
smodels <- lapply(lmodels, segmented)
breakpoints <- unlist(sapply(smodels, function(i) i$psi[, 'Est.']))
st.err <- round(unlist(sapply(smodels, function(i) i$psi[, 'St.Err'])), 1)
x_values <- lapply(smodels, function(i){
        seq(min(i$model[, 2]), max(i$model[, 2]), length.out = 1000)
})

xv <- lapply(seq_along(x_values), function(i){
        x_values[[i]][c(1, which(x_values[[i]] > breakpoints[i])[1]-1,
                        which(x_values[[i]] > breakpoints[i])[1], length(x_values[[i]]))]
})
pred_ <- lapply(seq_along(xv), function(i) predict(smodels[[i]],
                                                   newdata = data.frame(x = xv[[i]])))

## elaborated graph
png('G:/research/gits/AnTracks/plots/liquid_brains/fig_dist_scales.png', width = 6000, height = 3000,
    res = 400)
ggplot(data = fit_data, aes(x, y, fill = tag)) + 
        scale_x_continuous('Average distance to nest (cm)', breaks = seq(0, 140, 10))+
        scale_y_log10('Frequency (log-scale)')+ 
        facet_wrap(~factor(tag, labels = c('Experimental', 'Control')), scales = 'free')+
        geom_line(data = data.frame(x = unlist(xv),
                                    y = exp(unlist(pred_)),
                                    tag = c(rep('DET', length(xv[[1]])),
                                            rep('NFD', length(xv[[1]])))),
                  aes(x, y), linewidth = 1)+
        annotation_logticks(sides = 'l')+
        geom_vline(data = data.frame(x = breakpoints,
                                     tag = c('DET', 'NFD')),
                   linetype = 2, aes(xintercept = x), linewidth = 1)+
        geom_point(size = 6, shape = 21, color = 'black', alpha = 0.5) +
        scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
        theme(axis.text = element_text(size = 18),
              axis.title = element_text(size = 20),
              strip.text.x = element_text(size = 25),
              legend.position = 'none')

dev.off()
