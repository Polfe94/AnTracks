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

#### PLOTS FOR EACH PATTERN ####

#### THETA 10**-5 ####
path_theta05 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_05/'
theta05 <- process_data(path_theta05) # takes ~20 sec
dTheta05 <- rbindlist(lapply(theta05, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta05) <- knames
kTheta05 <- knames[apply(dTheta05, 1, which.min)]
plot_clusters(theta05, kTheta05)

#### THETA 10**-8 ####
# path_theta08 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_08/'
# theta08 <- process_data(path_theta08) # takes ~20 sec
# dTheta08 <- rbindlist(lapply(theta08, function(i){
# 	data.table(t(sapply(zSeries, function(x){
# 		ndtw(i[['Z']], x)
# 	})))
# })) # 5sec
# colnames(dTheta08) <- knames
# kTheta08 <- knames[apply(dTheta08, 1, which.min)]
# plot_clusters(theta08, kTheta08)

#### THETA 10**-10 ####
path_theta10 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_10/'
theta10 <- process_data(path_theta10) # takes ~20 sec
dTheta10 <- rbindlist(lapply(theta10, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta10) <- knames
kTheta10 <- knames[apply(dTheta10, 1, which.min)]
plot_clusters(theta10, kTheta10)

#### THETA 10**-13 ####
# path_theta13 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_13/'
# theta13 <- process_data(path_theta13) # takes ~20 sec
# dTheta13 <- rbindlist(lapply(theta13, function(i){
# 	data.table(t(sapply(zSeries, function(x){
# 		ndtw(i[['Z']], x)
# 	})))
# })) # 5sec
# colnames(dTheta13) <- knames
# kTheta13 <- knames[apply(dTheta13, 1, which.min)]
# plot_clusters(theta13, kTheta13)

#### THETA 10**-15 ####
path_theta15 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/uniform/'
theta15 <- process_data(path_theta15) # takes ~20 sec
dTheta15 <- rbindlist(lapply(theta15, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta15) <- knames
kTheta15 <- knames[apply(dTheta15, 1, which.min)]
plot_clusters(theta15, kTheta15)

#### THETA 10**-16 ####
# path_theta16 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_16/'
# theta16 <- process_data(path_theta16) # takes ~20 sec
# dTheta16 <- rbindlist(lapply(theta16, function(i){
# 	data.table(t(sapply(zSeries, function(x){
# 		ndtw(i[['Z']], x)
# 	})))
# })) # 5sec
# colnames(dTheta16) <- knames
# kTheta16 <- knames[apply(dTheta16, 1, which.min)]
# plot_clusters(theta16, kTheta16)

#### THETA 10**-17 ####
# path_theta17 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_17/'
# theta17 <- process_data(path_theta17) # takes ~20 sec
# dTheta17 <- rbindlist(lapply(theta17, function(i){
# 	data.table(t(sapply(zSeries, function(x){
# 		ndtw(i[['Z']], x)
# 	})))
# })) # 5sec
# colnames(dTheta17) <- knames
# kTheta17 <- knames[apply(dTheta17, 1, which.min)]
# plot_clusters(theta17, kTheta17)

#### THETA 10**-18 ####
# path_theta18 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_18/'
# theta18 <- process_data(path_theta18) # takes ~20 sec
# dTheta18 <- rbindlist(lapply(theta18, function(i){
# 	data.table(t(sapply(zSeries, function(x){
# 		ndtw(i[['Z']], x)
# 	})))
# })) # 5sec
# colnames(dTheta18) <- knames
# kTheta18 <- knames[apply(dTheta18, 1, which.min)]
# plot_clusters(theta18, kTheta18)

#### THETA 10**-19 ####
# path_theta19 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_19/'
# theta19 <- process_data(path_theta19) # takes ~20 sec
# dTheta19 <- rbindlist(lapply(theta19, function(i){
# 	data.table(t(sapply(zSeries, function(x){
# 		ndtw(i[['Z']], x)
# 	})))
# })) # 5sec
# colnames(dTheta19) <- knames
# kTheta19 <- knames[apply(dTheta19, 1, which.min)]
# plot_clusters(theta19, kTheta19)

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

#### THETA 10**-25 ####
path_theta25 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_25/'
theta25 <- process_data(path_theta25) # takes ~25 sec
dTheta25 <- rbindlist(lapply(theta25, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta25) <- knames
kTheta25 <- knames[apply(dTheta25, 1, which.min)]
plot_clusters(theta25, kTheta25)

#### THETA 10**-30 ####
path_theta30 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_30/'
theta30 <- process_data(path_theta30) # takes ~30 sec
dTheta30 <- rbindlist(lapply(theta30, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta30) <- knames
kTheta30 <- knames[apply(dTheta30, 1, which.min)]
plot_clusters(theta30, kTheta30)

#### THETA 10**-35 ####
path_theta35 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_35/'
theta35 <- process_data(path_theta35) # takes ~35 sec
dTheta35 <- rbindlist(lapply(theta35, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta35) <- knames
kTheta35 <- knames[apply(dTheta35, 1, which.min)]
plot_clusters(theta35, kTheta35)

#### THETA 10**-40 ####
path_theta40 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_40/'
theta40 <- process_data(path_theta40) # takes ~40 sec
dTheta40 <- rbindlist(lapply(theta40, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta40) <- knames
kTheta40 <- knames[apply(dTheta40, 1, which.min)]
plot_clusters(theta40, kTheta40)

#### THETA 10**-45 ####
path_theta45 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_45/'
theta45 <- process_data(path_theta45) # takes ~45 sec
dTheta45 <- rbindlist(lapply(theta45, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta45) <- knames
kTheta45 <- knames[apply(dTheta45, 1, which.min)]
plot_clusters(theta45, kTheta45)

#### THETA 10**-50 ####
path_theta50 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_50/'
theta50 <- process_data(path_theta50) # takes ~50 sec
dTheta50 <- rbindlist(lapply(theta50, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta50) <- knames
kTheta50 <- knames[apply(dTheta50, 1, which.min)]
plot_clusters(theta50, kTheta50)

#### THETA 10**-60 ####
path_theta60 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_60/'
theta60 <- process_data(path_theta60) # takes ~60 sec
dTheta60 <- rbindlist(lapply(theta60, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta60) <- knames
kTheta60 <- knames[apply(dTheta60, 1, which.min)]
plot_clusters(theta60, kTheta60)

#### THETA 10**-70 ####
path_theta70 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_70/'
theta70 <- process_data(path_theta70) # takes ~70 sec
dTheta70 <- rbindlist(lapply(theta70, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta70) <- knames
kTheta70 <- knames[apply(dTheta70, 1, which.min)]
plot_clusters(theta70, kTheta70)

#### THETA 10**-80 ####
path_theta80 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_80/'
theta80 <- process_data(path_theta80) # takes ~80 sec
dTheta80 <- rbindlist(lapply(theta80, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta80) <- knames
kTheta80 <- knames[apply(dTheta80, 1, which.min)]
plot_clusters(theta80, kTheta80)

#### THETA 10**-90 ####
path_theta90 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_90/'
theta90 <- process_data(path_theta90) # takes ~90 sec
dTheta90 <- rbindlist(lapply(theta90, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 5sec
colnames(dTheta90) <- knames
kTheta90 <- knames[apply(dTheta90, 1, which.min)]
plot_clusters(theta90, kTheta90)

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

global_Eff <- rbind(foodEff(path_theta05, 'theta_05'),
		    #foodEff(path_theta08, 'theta_08'),
		    foodEff(path_theta10, 'theta_10'),
		    #foodEff(path_theta13, 'theta_13'),
		    foodEff(path_theta15, 'theta_15'),
		    #foodEff(path_theta16, 'theta_16'),
		    #foodEff(path_theta17, 'theta_17'),
		    #foodEff(path_theta18, 'theta_18'),
		    #foodEff(path_theta19, 'theta_19'),
		    foodEff(path_theta20, 'theta_20'),
		    foodEff(path_theta25, 'theta_25'),
		    foodEff(path_theta30, 'theta_30'),
		    foodEff(path_theta35, 'theta_35'),
		    foodEff(path_theta40, 'theta_40'),
		    foodEff(path_theta45, 'theta_45'),
		    foodEff(path_theta50, 'theta_50'),
		    foodEff(path_theta60, 'theta_60'),
		    foodEff(path_theta70, 'theta_70'),
		    foodEff(path_theta70, 'theta_80'),
		    foodEff(path_theta70, 'theta_90'))


global_eff_mlt <- rbind(reshape2::melt(global_Eff[, .(eff_1 = 1/ mint, 
							eff_2 = 1/(maxt-mint),
							.id = mean(.id)), by = 'label'],
					id.vars = c('label', '.id')), exp_Eff_mlt)

ggplot(data = global_Eff[, .(eff_1 = mean(1/mint), eff_2 = mean(1/(maxt-mint))), by = 'label'],
       aes(x = as.numeric(gsub('theta_', '',label)), y = eff_1)) + geom_point()+
	ylab('Average exploration efficiency')+ xlab('Theta')+
	geom_smooth(se = FALSE, method = 'gam', formula = y ~ log(as.numeric(gsub('theta_', '',x))))
ggplot(data = global_Eff[, .(eff_1 = median(1/mint), eff_2 = median(1/(maxt-mint))), by = 'label'],
       aes(x = as.numeric(gsub('theta_', '',label)), y = eff_1)) + geom_point()+
	ylab('Average exploitation efficiency')+ xlab('Theta')+
	geom_smooth(se = FALSE, method = 'gam', formula = y ~ log(as.numeric(gsub('theta_', '',x))))


######
png('/home/polfer/research/gits/AnTracks/plots/median_eff_thetas.png', width = 3840, height = 2160,
    res = 250)


data_plot <- data.table(reshape2::melt(global_Eff[label %in% c('theta_05','theta_10',
						    'theta_15','theta_20', 'theta_25',
						    'theta_30', 'theta_35', 'theta_40',
						    'theta_45', 'theta_50', 'theta_60', 'theta_70',
						    'theta_80', 'theta_90'), 
						    .(eff_1 = median(1/mint), eff_2 = median(1/(maxt-mint))),
						    by = 'label'], id.vars = 'label'))

# data_plot2 <- data.table(reshape2::melt(global_Eff[label %in% c('theta_05','theta_10',
# 							       'theta_15','theta_20', 'theta_25',
# 							       'theta_30', 'theta_35', 'theta_40',
# 							       'theta_45', 'theta_50', 'theta_60', 'theta_70') &
# 						   	is.finite(mint), 
# 						  .(eff_1 = median(1/mint), eff_2 = median(1/(maxt-mint))),
# 						  by = 'label'], id.vars = 'label'))
set(data_plot, j = 'q10', value = c(global_Eff[label %in% unique(data_plot$label), 
					       .(q10 = quantile(1/mint, p = 0.25)), by = 'label'][['q10']],
				    global_Eff[label %in% unique(data_plot$label), 
				    	   .(q10 = quantile(1/(maxt-mint), p = 0.25)), by = 'label'][['q10']]))
set(data_plot, j = 'q90', value = c(global_Eff[label %in% unique(data_plot$label), 
					       .(q90 = quantile(1/mint, p = 0.75)), by = 'label'][['q90']],
				    global_Eff[label %in% unique(data_plot$label), 
				    	   .(q90 = quantile(1/(maxt-mint), p = 0.75)), by = 'label'][['q90']]))

ggplot(data = data_plot,
       aes(x = as.numeric(gsub('theta_', '',label)), y = value)) + 
	geom_point()+
	geom_pointrange(aes(ymin = q10, ymax = q90))+
	ylab('')+ 
	scale_x_continuous(TeX('Risk aversion ($\\Theta$)'), breaks = seq(0, 100, 10), 
			   labels = function(i){
			   	TeX(paste0('$10^{-', i,'}$'))
			   })+
	geom_smooth(se = FALSE, method = 'gam', formula = y ~ log(as.numeric(gsub('theta_', '',x))))+
	theme(aspect.ratio = 0.75)+facet_wrap(~ factor(variable,
						       labels = c('Median exploration efficiency',
									 'Median exploitation efficiency')),
									 scales = 'free')

dev.off()


ggplot(data = global_Eff[label %in% c('theta_05','theta_10',
				      'theta_15','theta_20', 'theta_25',
				      'theta_30', 'theta_35', 'theta_40',
				      'theta_45', 'theta_50', 'theta_60', 'theta_70',
				      'theta_80', 'theta_90'), 
			 .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint)),
			 by = 'label'], aes(as.numeric(gsub('theta_', '', label)), eff_1)) +
	geom_jitter(alpha = 0.5) + 
	stat_summary(fun = 'median', geom = 'point', shape = 23, 
		     fill = 'blue', color = 'black', size = 3)+
	stat_summary(fun = 'mean', geom = 'point', shape = 23, 
		     fill = 'red', color = 'black', size = 3)+
	geom_smooth(se = FALSE, method = 'gam', formula = y ~ log(as.numeric(gsub('theta_', '',x))))


### PROBABILITY OF FINDING AT LEAST 1%, 50% AND 100% OF THE FOOD
global_Eff[label %in% c('theta_05','theta_10',
			'theta_15','theta_20', 'theta_25',
			'theta_30', 'theta_35', 'theta_40',
			'theta_45', 'theta_50', 'theta_60', 'theta_70',
			'theta_80', 'theta_90'), 
			.(p1 = sum(abs(1-as.integer((na)/12)))/100,
			  p50 = sum(abs(1-round((na)/12.1)))/100,
			  p100 = sum(abs(1-ceiling((na)/12.1)))/100),
			by = 'label']

label_data <- reshape2::melt(global_Eff[label %in% c('theta_05','theta_10',
			'theta_15','theta_20', 'theta_25',
			'theta_30', 'theta_35', 'theta_40',
			'theta_45', 'theta_50', 'theta_60', 'theta_70',
			'theta_80', 'theta_90'), 
			.(eff_1 = 1/mint, eff_2 = 1/(maxt-mint), na = abs(1-round((na)/12.1))),
			by = 'label'][, .(p = sum(na)/100,
					  eff_1 = runif(length(unique(label)),
					  	      min = 0.019, max = 0.021),
					  eff_2 = runif(length(unique(label)), 
					  	      min = 0.0019, max = 0.0021)),
					  by = 'label'],
			id.vars = c('label',  'p'))

ggplot(data = rbind(data.frame(value = 0.1, variable = c('eff_1', 'eff_2'), label = 'theta_0'),
		    reshape2::melt(global_Eff[label %in% c('theta_05','theta_10',
		    				       'theta_15','theta_20', 'theta_25',
		    				       'theta_30', 'theta_35', 'theta_40',
		    				       'theta_45', 'theta_50', 'theta_60', 'theta_70',
		    				       'theta_80', 'theta_90'), 
		    			  .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint)),
		    			  by = 'label'], id.vars = c('label'))),
       aes(x = as.numeric(gsub('theta_', '',label)),
           y = value, group = label))+
	geom_smooth(aes(group = 1), method = 'gam', formula = y ~ log(x)) +
	facet_wrap(~ variable, scales = 'free')


exp_Eff_mlt$label <- 'experimental'
exp_Eff_mlt$label <- 'theta_-5'
##### PLOT VALID (?)
png('~/research/gits/AnTracks/plots/eff_by_theta.png', 6000, 4000, res = 450)
ggplot(data = rbind(reshape2::melt(global_Eff[label %in% c('theta_05','theta_10',
			       'theta_15','theta_20', 'theta_25',
			       'theta_30', 'theta_35', 'theta_40',
			       'theta_45', 'theta_50', 'theta_60', 'theta_70',
			       'theta_80', 'theta_90'), 
		  .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint), .id = 1),
		  by = 'label'], id.vars = c('label', '.id')),
		  exp_Eff_mlt),
       aes(x = as.numeric(gsub('theta_', '',label)),
		  		   y = value, group = label))+ 
	geom_boxplot(aes(x = as.numeric(gsub('theta_', '',label))),
		     outlier.shape = NA, color = 'grey55')+
	geom_jitter(alpha = 0.25, width = 1)+
	geom_smooth(data =
				 reshape2::melt(global_Eff[label %in% c('theta_05','theta_10',
								  'theta_15','theta_20', 'theta_25',
								  'theta_30', 'theta_35', 'theta_40',
								  'theta_45', 'theta_50', 'theta_60', 'theta_70',
								  'theta_80', 'theta_90'), 
								  .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint)),
								  by = 'label'], id.vars = c('label')),
		    se = FALSE, method = 'gam', color = 'darkslateblue',
		    formula = y ~ log(x), aes(group = 1))+
	#geom_label(data = label_data, aes(label = p, x = as.numeric(gsub('theta_', '',label))))+
	stat_summary(fun = 'mean', geom = 'point', shape = 23, 
		     fill = 'dodgerblue3', color = 'black', size = 3.5)+
	scale_y_continuous(TeX('Efficiency ($s^{-1}$)'), breaks = function(i){
		seq(0, round(i[2], 3), length.out = 5)
	}) + theme(aspect.ratio = 1)+
	scale_x_continuous(TeX('Risk aversion ($\\Theta$)'), breaks = c(-5, seq(10, 100, 10)),
			   labels = function(i){
			   	c('Exp', TeX(paste0('$10^{-',i, '}$'))[-1])
			   	})+
	facet_wrap(~factor(variable, labels = c('Exploration', 'Exploitation')),
		   ncol = 2, nrow = 1, scales = 'free')+
	geom_vline(linetype = 'dashed', xintercept =0)

dev.off()
# ggplot(data = global_Eff[label %in% c('theta_05','theta_10',
# 				      'theta_15','theta_20', 'theta_25',
# 				      'theta_30', 'theta_35', 'theta_40',
# 				      'theta_45', 'theta_50', 'theta_60', 'theta_70'), 
# 			 .(eff_1 = median(1/mint), eff_2 = median(1/(maxt-mint))), by = 'label'],
#        aes(x = as.numeric(gsub('theta_', '',label)), y = eff_2)) + geom_point()+
# 	scale_x_continuous('Theta', breaks = seq(0, 50, 5), 
# 			   labels = function(i){
# 			   	TeX(paste0('$10^{-', i,'}$'))
# 			   })+ ylab('Median exploitation efficiency') +
# 	geom_smooth(se = FALSE, method = 'gam', formula = y ~ log(as.numeric(gsub('theta_', '',x))))+
# 	theme(aspect.ratio = 0.75)

dev.off()
ggplot(data = global_eff_mlt[global_eff_mlt$label %in% c('experimental', 'theta_05','theta_10',
							    'theta_15','theta_20', 'theta_25',
							    'theta_30', 'theta_35', 'theta_40',
							    'theta_45', 'theta_50', 'theta_60',
							 'theta_70'),],
       aes(factor(label), value)) +
	facet_wrap(~ factor(variable, labels = c('Exploration Efficiency',
						 'Exploitation Efficiency')),
		   scales = 'free') + geom_boxplot() +
	scale_x_discrete(TeX('Risk propensity ($\\Theta$)'), labels = function(i){

		c('Exp',
		TeX(paste0('$10^{-', gsub('theta_', '', i), '}$'))[-1])

		}) + ylab('Search efficiency')

png('~/research/gits/AnTracks/plots/theta_efficiency.png', width = 4000, height = 2000, res = 160)
ggplot(data = global_eff_mlt, aes(factor(label, labels = c('Exp', 
							   r'($\Theta = 10^{-05}$)',
							   r'($\Theta = 10^{-08}$)',
							   r'($\Theta = 10^{-10}$)',
							   r'($\Theta = 10^{-13}$)',
							   r'($\Theta = 10^{-16}$)',
							   r'($\Theta = 10^{-17}$)',
							   r'($\Theta = 10^{-18}$)',
							   r'($\Theta = 10^{-19}$)',
							   r'($\Theta = 10^{-20}$)',
							   r'($\Theta = 10^{-25}$)',
							   r'($\Theta = 10^{-30}$)',
							   r'($\Theta = 10^{-35}$)',
							   r'($\Theta = 10^{-40}$)',
							   r'($\Theta = 10^{-45}$)',
							   r'($\Theta = 10^{-50}$)')), 
				  value)) +
	facet_wrap(~ factor(variable, labels = c('Exploration Efficiency',
						 'Exploitation Efficiency')),
		   scales = 'free') + geom_boxplot() +
	scale_x_discrete('', labels = TeX) + ylab('')
dev.off()
# takes 1 min
sims_N <- rbind(n_data(path_theta05, label = 'theta_05'),
		n_data(path_theta08, label = 'theta_08'),
		n_data(path_theta10, label = 'theta_10'),
		n_data(path_theta13, label = 'theta_13'),
		n_data(path_theta16, label = 'theta_16'),
		n_data(path_theta17, label = 'theta_17'),
		n_data(path_theta18, label = 'theta_18'),
		n_data(path_theta19, label = 'theta_19'),
		n_data(path_theta20, label = 'theta_20'),
		n_data(path_theta25, label = 'theta_25'),
		n_data(path_theta30, label = 'theta_30'),
		n_data(path_theta35, label = 'theta_35'),
		n_data(path_theta40, label = 'theta_40'),
		n_data(path_theta45, label = 'theta_45'),
		n_data(path_theta50, label = 'theta_50'))

##### EFFICIENCY INCLUDING EXPERIMENTAL VALUES
global_Eff[['Frame_min']] <- round(global_Eff[['mint']]*2)
global_Eff[['Frame_max']] <- round(global_Eff[['maxt']]*2)

global_Eff[['N1']] <- vapply(1:nrow(global_Eff), function(i){
	mint = global_Eff[i, Frame_min]
	maxt = global_Eff[i, Frame_max]
	mean(sims_N[.id == global_Eff[i, .id] &
		    	label == global_Eff[i, label] &
		    	Frame >= (mint-1) & Frame <= (mint+1), N])
}, numeric(1))

global_Eff[['N2']] <- vapply(1:nrow(global_Eff), function(i){
	mint = global_Eff[i, Frame_min]
	maxt = global_Eff[i, Frame_max]
	mean(sims_N[.id == global_Eff[i, .id] &
		    	label == global_Eff[i, label] &
		    	Frame >= (mint-1) & Frame <= (maxt+1), N])
}, numeric(1))

for(i in 1:nrow(sbst)){
	idx <- which(global_Eff[['.id']] == sbst[i, .id] &
		global_Eff[['mint']] == sbst[i, mint] &
		global_Eff[['label']] == sbst[i, label])
	global_Eff[idx, 'N1'] <- 1
}

global_EffxCapita_mlt <- rbind(reshape2::melt(global_Eff[, .(eff_1 = 1/(mint*N1), 
						       eff_2 = 1/((maxt-mint)*N2),
						       label = label, .id = .id)],
					      id.vars = c('label', '.id')),
				exp_EffxCapita_mlt)

ggplot(data = global_EffxCapita_mlt, aes(factor(label, labels = c('Exp', 
								  r'($10^{-05}$)',
								  r'($10^{-08}$)',
								  r'($10^{-10}$)',
								  r'($10^{-13}$)',
								  r'($10^{-16}$)',
								  r'($10^{-17}$)',
								  r'($10^{-18}$)',
								  r'($10^{-19}$)',
								  r'($10^{-20}$)',
								  r'($10^{-25}$)',
								  r'($10^{-30}$)',
								  r'($10^{-35}$)',
								  r'($10^{-40}$)',
								  r'($10^{-45}$)',
								  r'($10^{-50}$)')),
					 value)) +
	facet_wrap(~ factor(variable, labels = c('Exploration Efficiency (per capita)',
						 'Exploitation Efficiency (per capita)')),
		   scales = 'free') + geom_boxplot()+
	ylab('') + scale_x_discrete('', labels = TeX)





######## WILCOXON TESTS ########
wilcox.test(global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_1' & 
				  	global_EffxCapita_mlt$label == 'theta_16', 'value'],
	    global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_1' & 
	    	    	global_EffxCapita_mlt$label == 'experimental', 'value'])

wilcox.test(global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_2' & 
				   	global_EffxCapita_mlt$label == 'theta_16', 'value'],
	    global_EffxCapita_mlt[global_EffxCapita_mlt$variable == 'eff_2' & 
	    		       	global_EffxCapita_mlt$label == 'experimental', 'value'])
