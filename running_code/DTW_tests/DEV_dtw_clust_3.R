#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')
load('~/research/gits/AnTracks/data/nf.RData')
ref_flat <- data.table(read.csv('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameter_space/bimodal_10/bimodal_10_65.csv'))
ndtw <- function(x, y, ...) {
	dtw::dtw(x, y, step.pattern = asymmetric,
		 distance.only = TRUE, ...)$normalizedDistance
}

Rcpp::cppFunction('
NumericVector fillVec(NumericVector x) {
    for (int i = 1; i < x.size(); i++) {
        if (NumericVector::is_na(x[i])) {
            x[i] = x[i - 1];
        }
    }
    return x;
}
')

process_data <- function(path, 
			 isFile = FALSE,
			 vars = c('Frame', 'N')){
	ref <- data.table(Frame = 1:21600)
	t <- movingAverage(ref[['Frame']], 60, 0)
	if(isFile){
		files <- tail(strsplit(path, '/')[[1]], 1)
		path <- gsub(files, '', path)
	} else {
		files <- list.files(path)
		files <- files[!grepl('food', files) & !grepl('position', files) & !grepl('data', files)]
	}
	lapply(files, function(i){
		y <- merge(ref, data.table(read.csv(paste0(path, i)))[, ..vars], all = TRUE, by = 'Frame')
		y[['N']] <- fillVec(y[['N']])
		n <- movingAverage(y[['N']], 60L, 0L)
		z <- zscore(n)
		data.table(Frame = t, N = n, Z = z)
	})
}

plot_clusters <- function(data_list, clusters){
	kprop <- 100* (table(clusters) / length(clusters))
	n <- names(kprop)
	xPlot <- rbindlist(lapply(seq_along(data_list), function(i){
		data.table(Frame = data_list[[i]][['Frame']], N = data_list[[i]][['N']], k = clusters[i], exp = i)
	}))
	
	ggplot(data = xPlot, aes(Frame, N)) + geom_path(aes(group = exp), color = 'grey70', alpha = 0.75)+
		geom_path(data = xPlot[, .(N = mean(N), exp = exp), by = c('Frame', 'k')], linewidth = 1)+
		facet_wrap(~factor(k, levels = c('norm', 'flat', 'late', 'low'),
				   labels = c(paste0('Experimental pattern (', 
				   		  ifelse('norm' %in% n, kprop[['norm']], 0),
				   		  '%)'), 
				   	   paste0('Plateau (', 
				   	          ifelse('flat' %in% n, kprop[['flat']], 0),
				   	          '%)'), 
				   	   paste0('Late start (', 
				   	          ifelse('late' %in% n, kprop[['late']], 0),
				   	          '%)'), 
				   	   paste0('Low activity (', 
				   	          ifelse('low' %in% n, kprop[['low']], 0),
				   	          '%)'))), 
			   nrow = 1) +
		ylab('Activity (number of ants in arena)') + 
		scale_x_continuous('Time (min)', breaks = seq(0, 21600, 30 * 120), labels = seq(0, 180, 30))+
		theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
		      aspect.ratio = 0.75)
}

foodEff <- function(data_path, label = NULL){
	files <- list.files(data_path)
	files <- files[grepl('food', files)]
	
	result <- rbindlist(lapply(seq_along(files), function(i){
		x <- read.csv(paste0(data_path, files[i]))[['t']]
		data.table(mint = min(x, na.rm = TRUE), 
			   maxt = max(x, na.rm = TRUE),
			   na = sum(is.na(x)),
			   inter_mean = mean(diff(sort(x))))
	}), idcol = TRUE)
	
	if(!is.null(label)){
		set(result, j = 'label', value = label)
	}
	result
}

knames <- c('norm', 'late', 'flat', 'low')
ref <- data.table(Frame = seq_len(21600))

avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(ref, i@data[, .N, by = 'Frame'], by = 'Frame', all = TRUE)[1:21600]
	x[is.na(N), 'N'] <- 0
	x
}))[, .(N = mean(N)), by = 'Frame']


exp_Eff <- rbind(rbindlist(lapply(det, function(i){rbindlist(i@food)})),
		 rbindlist(lapply(sto, function(i){rbindlist(i@food)})))
exp_Eff[['.id']] <- unlist(lapply(1:20, function(i) rep(i, 12)))

exp_Eff_mlt <- reshape2::melt(exp_Eff[, .(eff_1 = 1/ min(t/2), 
					  eff_2 = 1/(max(t/2)-min(t/2)),
					  label = 'experimental'), by = .id], id.vars = c('.id', 'label'))

exp_Eff_mlt <- reshape2::melt(exp_Eff[, .(eff_1 = 1/ min(t/2), 
					  eff_2 = 1/(max(t/2)-min(t/2)),
					  label = 'experimental'), by = .id], id.vars = c('.id', 'label'))


#### REFERENCE PATTERNS ####
zSeries <- list(norm = zscore(movingAverage(avg[['N']], t = 60L, overlap = 0L)))

sto_ref <- merge(ref, data.table(Frame = sto[[2]]@data[, 'Frame']), by = 'Frame', all = TRUE)[, .N, by = 'Frame'][1:21600]
sto_ref[is.na(N), 'N'] <- 0
low_ref <- merge(ref, data.table(Frame = nf[[6]]@data[, 'Frame']), by = 'Frame', all = TRUE)[, .N, by = 'Frame'][1:21600]
low_ref[is.na(N), 'N'] <- 0
flat_ref <- merge(ref, ref_flat[, c('Frame', 'N')], by = 'Frame', all = TRUE)[1:21600]
flat_ref[['N']] <- fillVec(flat_ref[['N']])


zSeries[['late']] <- zscore(movingAverage(sto_ref[['N']], t = 60L, overlap = 0L))
zSeries[['flat']] <- zscore(movingAverage(flat_ref[['N']], t = 60L, overlap = 0L))
zSeries[['low']] <- zscore(movingAverage(low_ref[['N']], t = 60L, overlap = 0L))

#### EXPERIMENTAL EFFICIENCY

exp_Eff_xCapita <- rbind(rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- rbindlist(i@food)
	mint <- min(x[['t']])
	maxt <- max(x[['t']])
	x[['N1']] <- i@data[Frame == mint, .N]
	x[['N2']] <- mean(i@data[Frame >= mint & Frame <= maxt, .N, by = 'Frame'][['N']])
	x
})),rbindlist(lapply(sto, function(i){
	setDT(i@data)
	x <- rbindlist(i@food)
	mint <- min(x[['t']])
	maxt <- max(x[['t']])
	x[['N1']] <- i@data[Frame == mint, .N]
	x[['N2']] <- mean(i@data[Frame >= mint & Frame <= maxt, .N, by = 'Frame'][['N']])
	x
})))
exp_Eff_xCapita[['.id']] <- unlist(lapply(1:20, function(i) rep(i, 12)))

exp_EffxCapita_mlt <- reshape2::melt(exp_Eff_xCapita[, .(eff_1 = 1/ (min(t/2)*mean(N1)), 
							 eff_2 = 1/((max(t/2)-min(t/2))*mean(N2)),
							 label = 'experimental'), by = .id], id.vars = c('.id', 'label'))


## clean environment
rm(det, sto, nf)
gc()

path_results <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/uniform/'
data <- process_data(path_results)
d <- rbindlist(lapply(data, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(d) <- knames
k <- knames[apply(d, 1, which.min)]
plot_clusters(data, k)


eff_data <- foodEff(path_results, label = 'uniform')
ggplot(data = reshape2::melt(eff_data[, .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint))]))



##################
files <- list.files(path_results)
dt <- data.table(read.csv(paste0(path_results, files[grepl('data', files)][7])))[-1, ]

dt2 <- dt[, .N, by = c('Frame', 'id')][, .(N = .N, ids = list(id)), by = 'Frame']
dt2 <- merge(data.table(Frame = 1:21600), dt2, all = TRUE, by = 'Frame')
dt2[sapply(ids, function(i)length(i) == 0), N := 0]
plot(dt2[, c('Frame',  'N')], type = 'l')
