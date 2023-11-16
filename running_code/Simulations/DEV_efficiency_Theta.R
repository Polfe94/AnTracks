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

## clean environment
rm(det, sto, nf)
gc()


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

#### THETA 10**-16 ####
path_theta16 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_16/'
theta16 <- process_data(path_theta16) # takes ~20 sec
dTheta16 <- rbindlist(lapply(theta16, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta16) <- knames
kTheta16 <- knames[apply(dTheta16, 1, which.min)]
plot_clusters(theta16, kTheta16)

#### THETA 10**-17 ####
path_theta17 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_17/'
theta17 <- process_data(path_theta17) # takes ~20 sec
dTheta17 <- rbindlist(lapply(theta17, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta17) <- knames
kTheta17 <- knames[apply(dTheta17, 1, which.min)]
plot_clusters(theta17, kTheta17)

#### THETA 10**-18 ####
path_theta18 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_18/'
theta18 <- process_data(path_theta18) # takes ~20 sec
dTheta18 <- rbindlist(lapply(theta18, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta18) <- knames
kTheta18 <- knames[apply(dTheta18, 1, which.min)]
plot_clusters(theta18, kTheta18)

#### THETA 10**-19 ####
path_theta19 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_19/'
theta19 <- process_data(path_theta19) # takes ~20 sec
dTheta19 <- rbindlist(lapply(theta19, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta19) <- knames
kTheta19 <- knames[apply(dTheta19, 1, which.min)]
plot_clusters(theta19, kTheta19)

#### THETA 10**-20 ####
path_theta20 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_20/'
theta20 <- process_data(path_theta20) # takes ~20 sec
dTheta20 <- rbindlist(lapply(theta20, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta20) <- knames
kTheta20 <- knames[apply(dTheta20, 1, which.min)]
plot_clusters(theta20, kTheta20)


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


global_Eff <- rbind(foodEff(path_theta16, 'theta_16'),
		    foodEff(path_theta16, 'theta_17'),
		    foodEff(path_theta16, 'theta_18'),
		    foodEff(path_theta16, 'theta_19'),
		    foodEff(path_theta16, 'theta_20'))

global_eff_mlt <- rbind(reshape2::melt(global_Eff[, .(eff_1 = 1/ mint, 
							eff_2 = 1/(maxt-mint),
							label = label), by = .id],
					id.vars = c('label', '.id')), exp_Eff_mlt)

# takes 1 min
sims_N <- rbind(n_data(path_theta16, label = 'theta_16'),
		n_data(path_theta17, label = 'theta_17'),
		n_data(path_theta18, label = 'theta_18'),
		n_data(path_theta19, label = 'theta_19'),
		n_data(path_theta20, label = 'theta_20'))

##### EFFICIENCY INCLUDING EXPERIMENTAL VALUES
global_Eff_eff[['Frame_min']] <- round(global_Eff_eff[['mint']]*2)
global_Eff_eff[['Frame_max']] <- round(global_Eff_eff[['maxt']]*2)

global_Eff_eff[['N1']] <- vapply(1:nrow(global_Eff_eff), function(i){
	mint = global_Eff_eff[i, Frame_min]
	maxt = global_Eff_eff[i, Frame_max]
	mean(sims_N[.id == global_Eff_eff[i, .id] &
		    	label == global_Eff_eff[i, label] &
		    	Frame == mint, N])
}, numeric(1))

global_Eff_eff[['N2']] <- vapply(1:nrow(global_Eff_eff), function(i){
	mint = global_Eff_eff[i, Frame_min]
	maxt = global_Eff_eff[i, Frame_max]
	mean(sims_N[.id == global_Eff_eff[i, .id] &
		    	label == global_Eff_eff[i, label] &
		    	Frame > mint & Frame <= maxt, N])
}, numeric(1))


global_EffxCapita_mlt <- rbind(reshape2::melt(global_Eff[, .(eff_1 = 1/(mint*N1), 
						       eff_2 = 1/((maxt-mint)*N2),
						       label = label, .id = .id)], id.vars = c('label', '.id')),
				exp_EffxCapita_mlt)

ggplot(data = global_EffxCapita_mlt, aes(label, value)) +
	facet_wrap(~ factor(variable, labels = c('Exploration Efficiency (per capita)',
						 'Exploitation Efficiency (per capita)')),
		   scales = 'free') + geom_boxplot()


wilcox.test(global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_1' & 
				  	global_EffxCapita_mlt$label == 'theta_16', 'value'],
	    global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_1' & 
	    	    	global_EffxCapita_mlt$label == 'experimental', 'value'])

wilcox.test(global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_2' & 
				   	global_EffxCapita_mlt$label == 'theta_16', 'value'],
	    global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_2' & 
	    		       	global_EffxCapita_mlt$label == 'experimental', 'value'])
