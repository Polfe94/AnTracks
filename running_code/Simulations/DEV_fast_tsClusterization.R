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
		files <- files[!grepl('food', files) & !grepl('position', files)]
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

knames <- c('norm', 'late', 'flat', 'low')
ref <- data.table(Frame = seq_len(21600))

avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(ref, i@data[, .N, by = 'Frame'], by = 'Frame', all = TRUE)[1:21600]
	x[is.na(N), 'N'] <- 0
	x
}))[, .(N = mean(N)), by = 'Frame']


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

## clean environment
rm(det, sto, nf)
gc()

#### SOME PLOTS ####
t <- movingAverage(seq_len(21600), 60)
lims <- range(sapply(zSeries, range))

plot(t, zSeries[[1]], type = 'l', ylim = lims, ylab = 'z-Transformed N', xlab = 'Frame')
lines(t, zSeries[[2]], col = 'red')
lines(t, zSeries[[3]], col = 'blue')
lines(t, zSeries[[4]], col = 'orange')


#### FIRST PLOT, DISTRIBUTIONS
distribs <- list(bimodal_10 = c(rbeta(10000, 2, 6), rbeta(90000, 6, 2)),
		 bimodal_90 = c(rbeta(90000, 2, 6), rbeta(10000, 6, 2)),
		 bimodal_50 = c(rbeta(50000, 2, 6), rbeta(50000, 6, 2)),
		 bimodal = c(rbeta(100000, 0.5, 0.5)),
		 gaussian = c(rbeta(100000, 2, 2)),
		 uniform = runif(100000, 0, 1))
distribs <- reshape2::melt(data.table(do.call('cbind.data.frame', distribs)))
Aplot <- ggplot(data = distribs[distribs$variable != 'uniform', ], aes(value)) + 
	geom_histogram(breaks = seq(0, 1, 0.04), aes(y = after_stat(density)), color = 'black', fill = 'grey90')+
	geom_density(aes(y = after_stat(density)), color = 'blue', linewidth = 1, bw = 0.025)+
	facet_wrap(~factor(variable, 
			   levels = c('bimodal', 'gaussian', 'bimodal_50', 'bimodal_90', 'bimodal_10'),
			   labels = c('$\\beta(0.5, 0.5)$',
			   	   '$\\beta(2, 2)$',
			   	   '$0.5 \\times \\beta(2, 6)+0.5 \\times \\beta(6, 2)$',
			   	   '$0.9 \\times\\beta(2, 6)+0.1 \\times \\beta(6, 2)$',
			   	   '$0.1 \\times\\beta(2, 6)+0.9 \\times \\beta(6, 2)$')),
		   labeller = as_labeller(function(i){TeX(paste0(i))},
		   		       default = label_parsed), scales = 'free',
		   ncol = 1, nrow = 5)+
	theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
	      strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')))+
	ylab('') + xlab('')

#### PLOTS FOR EACH PATTERN ####

#### +++++++ WITH RECRUITMENT +++++++ ####
#### BIMODAL 10/90 ####
path_bimodal_10 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_10/'
bimodal_10 <- process_data(path_bimodal_10) # takes ~20 sec
dBimodal_10 <- rbindlist(lapply(bimodal_10, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dBimodal_10) <- knames
kBimodal_10 <- knames[apply(dBimodal_10, 1, which.min)]
plot_clusters(bimodal_10, kBimodal_10)

#### BIMODAL 90/10 ####
path_bimodal_90 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_90/'
bimodal_90 <- process_data(path_bimodal_90) # takes ~20 sec
dBimodal_90 <- rbindlist(lapply(bimodal_90, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dBimodal_90) <- knames
kBimodal_90 <- knames[apply(dBimodal_90, 1, which.min)]
plot_clusters(bimodal_90, kBimodal_90)

#### BIMODAL 50/50 ####
path_bimodal_50 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_50/'
bimodal_50 <- process_data(path_bimodal_50) # takes ~20 sec
dBimodal_50 <- rbindlist(lapply(bimodal_50, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dBimodal_50) <- knames
kBimodal_50 <- knames[apply(dBimodal_50, 1, which.min)]
plot_clusters(bimodal_50, kBimodal_50)

#### BIMODAL ####
path_bimodal <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal/'
bimodal <- process_data(path_bimodal) # takes ~20 sec
dBimodal <- rbindlist(lapply(bimodal, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dBimodal) <- knames
kBimodal <- knames[apply(dBimodal, 1, which.min)]
plot_clusters(bimodal, kBimodal)

#### GAUSSIAN-LIKE ####
path_gaussian <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gaussian/'
gaussian <- process_data(path_gaussian) # takes ~20 sec
dGaussian <- rbindlist(lapply(gaussian, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dGaussian) <- knames
kGaussian <- knames[apply(dGaussian, 1, which.min)]
plot_clusters(gaussian, kGaussian)


Bplot <- ggarrange(plot_clusters(bimodal, kBimodal),
		   plot_clusters(gaussian, kGaussian),
		   plot_clusters(bimodal_50, kBimodal_50),
		   plot_clusters(bimodal_90, kBimodal_90),
		   plot_clusters(bimodal_10, kBimodal_10), nrow = 5, ncol = 1)


png('~/research/gits/AnTracks/plots/PLOTS_PROPORTIONS.png', width = 8000, height = 8000, res = 450)
ggarrange(Aplot, Bplot, widths = c(1/4, 3/4), labels = c('A', 'B'), 
	  font.label = list(size = 18, face = 'plain'))
dev.off()



#### +++++++ WITHOUT RECRUITMENT +++++++ ####
#### GAUSSIAN ####
path_gaussian_norec <- '/home/polfer/research/gits/AutomatAnts/results/without_recruitment/parameters/gaussian/'
gaussian_norec <- process_data(path_gaussian_norec) # takes ~20 sec
dGaussian_norec <- rbindlist(lapply(gaussian_norec, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dGaussian_norec) <- knames
kGaussian_norec <- knames[apply(dGaussian_norec, 1, which.min)]

plot_clusters(gaussian_norec, kGaussian_norec)

#### UNIFORM ####
path_uniform_norec <- '/home/polfer/research/gits/AutomatAnts/results/without_recruitment/parameters/uniform/'
uniform_norec <- process_data(path_uniform_norec) # takes ~20 sec
dUniform_norec <- rbindlist(lapply(uniform_norec, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dUniform_norec) <- knames
kUniform_norec <- knames[apply(dUniform_norec, 1, which.min)]
plot_clusters(uniform_norec, kUniform_norec)


#### PLOT NO RECRUITMENT ####
Aplot_norec <- ggplot(data = distribs[distribs$variable %in% c('gaussian', 'uniform'), ], aes(value)) + 
	geom_histogram(breaks = seq(0, 1, 0.04), 
		       aes(y = after_stat(density)), color = 'black', fill = 'grey90')+
	geom_density(aes(y = after_stat(density)), color = 'blue', linewidth = 1, bw = 0.05)+
	facet_wrap(~factor(variable, 
			   levels = c('gaussian', 'uniform'),
			   labels = c('$\\beta(2, 2)$',
			   	   '$U(0,1)')),
		   labeller = as_labeller(function(i){TeX(paste0(i))},
		   		       default = label_parsed), scales = 'free',
		   ncol = 1, nrow = 2)+
	theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
	      strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')))+
	ylab('') + scale_x_continuous('')


png('~/research/gits/AnTracks/plots/PLOTS_PROPORTIONS_NOREC.png', width = 4000, height = 1600, res = 180)
ggarrange(Aplot_norec, 
	  ggarrange(plot_clusters(gaussian_norec, kGaussian_norec), 
	  	  plot_clusters(uniform_norec, kUniform_norec), nrow = 2, ncol = 1),
	  widths = c(1/4, 3/4), labels = c('A', 'B'), 
	  font.label = list(size = 18, face = 'plain'))
dev.off()


#### +++++++ PROPORTIONS TEST +++++++ ####
# B(2,2) vs 0.5*B(2,6)+0.5*B(6,2)
prop.test(x = c(sum(kGaussian == 'norm'), sum(kBimodal_50 == 'norm')), 
	  n = c(length(kGaussian), length(kBimodal_50)))
# B(2,2) with recruitment vs B(2,2) without recruitment
prop.test(x = c(sum(kGaussian == 'norm'), sum(kGaussian_norec == 'norm')), 
	  n = c(length(kGaussian), length(kGaussian_norec)))
# B(2,2) with recruitment vs U(0,1) without recruitment
prop.test(x = c(sum(kGaussian == 'norm'), sum(kUniform_norec == 'norm')), 
	  n = c(length(kGaussian), length(kUniform_norec)))
# B(2,2) without recruitment vs U(0,1) without recruitment
prop.test(x = c(sum(kGaussian_norec == 'norm'), sum(kUniform_norec == 'norm')), 
	  n = c(length(kGaussian_norec), length(kUniform_norec)))




#### UNIFORM + RECRUITMENT ####
path_uniform <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/uniform/'
uniform <- process_data(path_uniform) # takes ~20 sec
dUniform <- rbindlist(lapply(uniform, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dUniform) <- knames
kUniform <- knames[apply(dUniform, 1, which.min)]
plot_clusters(uniform, kUniform)

