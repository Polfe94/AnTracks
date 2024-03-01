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

n_data <- function(path, label){
	
	files <- list.files(path)
	files <- files[!grepl('food', files) & !grepl('position', files)]
	
	rbindlist(lapply(files, function(i){
		x <- data.table(read.csv(paste0(path, i)))[, c('Frame', 'N')]
		x[['label']] <- label
		x
	}), idcol = TRUE)
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
		files <- files[!grepl('food', files) & !grepl('position', files) & 
			       	!grepl('data', files)]
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


exp_Eff <- rbind(rbindlist(lapply(det, function(i){rbindlist(i@food)})),
		 rbindlist(lapply(sto, function(i){rbindlist(i@food)})))
exp_Eff[['.id']] <- unlist(lapply(1:20, function(i) rep(i, 12)))



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
	x[['N2']] <- mean(i@data[Frame > mint & Frame <= maxt, .N, by = 'Frame'][['N']])
	x
})),rbindlist(lapply(sto, function(i){
	setDT(i@data)
	x <- rbindlist(i@food)
	mint <- min(x[['t']])
	maxt <- max(x[['t']])
	x[['N1']] <- i@data[Frame == mint, .N]
	x[['N2']] <- mean(i@data[Frame > mint & Frame <= maxt, .N, by = 'Frame'][['N']])
	x
})))
exp_Eff_xCapita[['.id']] <- unlist(lapply(1:20, function(i) rep(i, 12)))

exp_EffxCapita_mlt <- reshape2::melt(exp_Eff_xCapita[, .(eff_1 = 1/ (min(t/2)*N1), 
							 eff_2 = 1/((max(t/2)-min(t/2))*N2),
							 label = 'experimental'), by = .id], id.vars = c('.id', 'label'))


#### PLOTS FOR EACH PATTERN ####

#### JIJ = 0 ####
path_Jij0 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/homogeneous_Jij/Jij_0/uniform/'
Jij0 <- process_data(path_Jij0) # takes ~20 sec
dJij0 <- rbindlist(lapply(Jij0, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij0) <- knames
kJij0 <- knames[apply(dJij0, 1, which.min)]
plot_clusters(Jij0, kJij0)

#### Jij = 0.4 ####
path_Jij0.4 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/homogeneous_Jij/Jij_0.4/uniform/'
Jij0.4 <- process_data(path_Jij0.4) # takes ~20 sec
dJij0.4 <- rbindlist(lapply(Jij0.4, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij0.4) <- knames
kJij0.4 <- knames[apply(dJij0.4, 1, which.min)]
plot_clusters(Jij0.4, kJij0.4)

#### JIJ = 1 ####
path_Jij1 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/homogeneous_Jij/Jij_1/uniform/'
Jij1 <- process_data(path_Jij1) # takes ~20 sec
dJij1 <- rbindlist(lapply(Jij1, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij1) <- knames
kJij1 <- knames[apply(dJij1, 1, which.min)]
plot_clusters(Jij1, kJij1)


#### JIJ = 1.5 ####
path_Jij1.5 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/homogeneous_Jij/Jij_1.5/uniform/'
Jij1.5 <- process_data(path_Jij1.5) # takes ~20 sec
dJij1.5 <- rbindlist(lapply(Jij1.5, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij1.5) <- knames
kJij1.5 <- knames[apply(dJij1.5, 1, which.min)]
plot_clusters(Jij1.5, kJij1.5)


klstrs <- ggarrange(plot_clusters(Jij0, kJij0) + ggtitle('Jij = 0'),
		    plot_clusters(Jij0.4, kJij0.4) + ggtitle('Jij = 0.4'),
		    plot_clusters(Jij1, kJij1) + ggtitle('Jij = 1'),
		    plot_clusters(Jij1.5, kJij1.5) + ggtitle('Jij = 1.5'),
		    nrow = 4, ncol = 1)


#### Jij = 0.4 + 1 ####
path_Jij0.4_1 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/heterogeneous_Jij/Jij_0.4_1/uniform/'
Jij0.4_1 <- process_data(path_Jij0.4_1) # takes ~20 sec
dJij0.4_1 <- rbindlist(lapply(Jij0.4_1, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij0.4_1) <- knames
kJij0.4_1 <- knames[apply(dJij0.4_1, 1, which.min)]
plot_clusters(Jij0.4_1, kJij0.4_1)

#### Jij = 0 + 1.5 ####
path_Jij0_1.5 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/heterogeneous_Jij/Jij_0_1.5/uniform/'
Jij0_1.5 <- process_data(path_Jij0_1.5) # takes ~20 sec
dJij0_1.5 <- rbindlist(lapply(Jij0_1.5, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij0_1.5) <- knames
kJij0_1.5 <- knames[apply(dJij0_1.5, 1, which.min)]
plot_clusters(Jij0_1.5, kJij0_1.5)

#### Jij = 0.4 + 1.5 ####
path_Jij0.4_1.5 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/heterogeneous_Jij/Jij_0.4_1.5/uniform/'
Jij0.4_1.5 <- process_data(path_Jij0.4_1.5) # takes ~20 sec
dJij0.4_1.5 <- rbindlist(lapply(Jij0.4_1.5, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij0.4_1.5) <- knames
kJij0.4_1.5 <- knames[apply(dJij0.4_1.5, 1, which.min)]
plot_clusters(Jij0.4_1.5, kJij0.4_1.5)

#### Jij = 0.4 + 3 ####
path_Jij0.4_3 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/heterogeneous_Jij/Jij_0.4_3/uniform/'
Jij0.4_3 <- process_data(path_Jij0.4_3) # takes ~20 sec
dJij0.4_3 <- rbindlist(lapply(Jij0.4_3, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij0.4_3) <- knames
kJij0.4_3 <- knames[apply(dJij0.4_3, 1, which.min)]
plot_clusters(Jij0.4_3, kJij0.4_3)

#### Jij = 0.4 + 5 ####
path_Jij0.4_5 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/heterogeneous_Jij/Jij_0.4_5/uniform/'
Jij0.4_5 <- process_data(path_Jij0.4_5) # takes ~20 sec
dJij0.4_5 <- rbindlist(lapply(Jij0.4_5, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dJij0.4_5) <- knames
kJij0.4_5 <- knames[apply(dJij0.4_5, 1, which.min)]
plot_clusters(Jij0.4_5, kJij0.4_5)



v <- rbind(foodEff(path_Jij0, 'Jij_0'),
	   foodEff(path_Jij0.4, 'Jij_0.4'),
	   foodEff(path_Jij1, 'Jij_1'),
	   foodEff(path_Jij1.5, 'Jij_1.5'),
	   foodEff(path_Jij0_1.5, 'Jij_0_1.5'),
	   foodEff(path_Jij0.4_1, 'Jij_0.4_1'),
	   foodEff(path_Jij0.4_1.5, 'Jij_0.4_1.5'),
	   foodEff(path_Jij0.4_3, 'Jij_0.4_3'),
	   foodEff(path_Jij0.4_5, 'Jij_0.4_5'))

dtw <- rbindlist(list(dJij0, dJij0.4, dJij1, dJij1.5, dJij0_1.5, 
		      dJij0.4_1, dJij0.4_1.5 , dJij0.4_3,  dJij0.4_5), idcol = T)

global_eff_mlt <- rbind(reshape2::melt(v[, .(eff_1 = 1/ mint, 
					     eff_2 = 1/(maxt-mint),
					     .id = mean(.id)), by = 'label'],
					     id.vars = c('label', '.id')), exp_Eff_mlt)

ggplot(data = global_eff_mlt, aes(factor(label, 
					 labels = c('Exp', '0',
					 	   '0 + 1.5',
					 	   '0.4', '0.4 + 1',
					 	   '0.4 + 1.5',
					 	   
					 	   '0.4 + 3',
					 	   '0.4 + 5',
					 	   '1',
					 	   '1.5'
					 	   
					 	   
					 	  )), 
				  value)) +
	facet_wrap(~ factor(variable, labels = c('Exploration Efficiency',
						 'Exploitation Efficiency')),
		   scales = 'free') + geom_boxplot() +
	scale_x_discrete('Information value', labels = TeX) + ylab('')+
	coord_cartesian(ylim = c(0, 0.01))
