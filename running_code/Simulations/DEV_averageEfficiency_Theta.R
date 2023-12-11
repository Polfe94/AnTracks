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
		files <- files[!grepl('food', files) & !grepl('position', files) !grepl('data', files)]
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

#### THETA 10**-1 ####
path_theta01 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_01/'
theta01 <- process_data(path_theta01) # takes ~20 sec
dTheta01 <- rbindlist(lapply(theta01, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 1sec
colnames(dTheta01) <- knames
kTheta01 <- knames[apply(dTheta01, 1, which.min)]

#### THETA 10**-2 ####
path_theta02 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_02/'
theta02 <- process_data(path_theta02) # takes ~20 sec
dTheta02 <- rbindlist(lapply(theta02, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 2sec
colnames(dTheta02) <- knames
kTheta02 <- knames[apply(dTheta02, 1, which.min)]

#### THETA 10**-3 ####
path_theta03 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_03/'
theta03 <- process_data(path_theta03) # takes ~20 sec
dTheta03 <- rbindlist(lapply(theta03, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 3sec
colnames(dTheta03) <- knames
kTheta03 <- knames[apply(dTheta03, 1, which.min)]

#### THETA 10**-4 ####
path_theta04 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_04/'
theta04 <- process_data(path_theta04) # takes ~20 sec
dTheta04 <- rbindlist(lapply(theta04, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 4sec
colnames(dTheta04) <- knames
kTheta04 <- knames[apply(dTheta04, 1, which.min)]

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

#### THETA 10**-6 ####
path_theta06 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_06/'
theta06 <- process_data(path_theta06) # takes ~20 sec
dTheta06 <- rbindlist(lapply(theta06, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 6sec
colnames(dTheta06) <- knames
kTheta06 <- knames[apply(dTheta06, 1, which.min)]

#### THETA 10**-7 ####
path_theta07 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_07/'
theta07 <- process_data(path_theta07) # takes ~20 sec
dTheta07 <- rbindlist(lapply(theta07, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 7sec
colnames(dTheta07) <- knames
kTheta07 <- knames[apply(dTheta07, 1, which.min)]

#### THETA 10**-8 ####
path_theta08 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_08/'
theta08 <- process_data(path_theta08) # takes ~20 sec
dTheta08 <- rbindlist(lapply(theta08, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 8sec
colnames(dTheta08) <- knames
kTheta08 <- knames[apply(dTheta08, 1, which.min)]

#### THETA 10**-9 ####
path_theta09 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_09/'
theta09 <- process_data(path_theta09) # takes ~20 sec
dTheta09 <- rbindlist(lapply(theta09, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 9sec
colnames(dTheta09) <- knames
kTheta09 <- knames[apply(dTheta09, 1, which.min)]

#### THETA 10**-10 ####
path_theta10 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_10/'
theta10 <- process_data(path_theta10) # takes ~20 sec
dTheta10 <- rbindlist(lapply(theta10, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 10sec
colnames(dTheta10) <- knames
kTheta10 <- knames[apply(dTheta10, 1, which.min)]

#### THETA 10**-15 ####
path_theta15 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/uniform/'
theta15 <- process_data(path_theta15) # takes ~20 sec
dTheta15 <- rbindlist(lapply(theta15, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 15sec
colnames(dTheta15) <- knames
kTheta15 <- knames[apply(dTheta15, 1, which.min)]

#### THETA 10**-20 ####
path_theta20 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_20/'
theta20 <- process_data(path_theta20) # takes ~20 sec
dTheta20 <- rbindlist(lapply(theta20, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 20sec
colnames(dTheta20) <- knames
kTheta20 <- knames[apply(dTheta20, 1, which.min)]

#### THETA 10**-25 ####
path_theta25 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_25/'
theta25 <- process_data(path_theta25) # takes ~25 sec
dTheta25 <- rbindlist(lapply(theta25, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 25sec
colnames(dTheta25) <- knames
kTheta25 <- knames[apply(dTheta25, 1, which.min)]

#### THETA 10**-30 ####
path_theta30 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/theta/theta_30/'
theta30 <- process_data(path_theta30) # takes ~30 sec
dTheta30 <- rbindlist(lapply(theta30, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
})) # 30sec
colnames(dTheta30) <- knames
kTheta30 <- knames[apply(dTheta30, 1, which.min)]


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


global_Eff <- rbind(foodEff(path_theta01, 'theta_01'),
		    foodEff(path_theta02, 'theta_02'),
		    foodEff(path_theta03, 'theta_03'),
		    foodEff(path_theta04, 'theta_04'),
		    foodEff(path_theta05, 'theta_05'),
		    foodEff(path_theta06, 'theta_06'),
		    foodEff(path_theta07, 'theta_07'),
		    foodEff(path_theta08, 'theta_08'),
		    foodEff(path_theta09, 'theta_09'),
		    foodEff(path_theta10, 'theta_10'),
		    foodEff(path_theta15, 'theta_15'),
		    foodEff(path_theta20, 'theta_20'),
		    foodEff(path_theta25, 'theta_25'))


# global_Eff[!is.finite(mint), mint := 10800]
# global_Eff[!is.finite(maxt), maxt := 10800]
data_plot <- global_Eff[label %in% paste0('theta_',
			     sprintf('%02d', c(1:10, 15,20,25,30))),
	   .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint), .id = 1),
	   by = 'label']
# data_plot[!is.finite(eff_2), eff_2 := 0]
exp_Eff_mlt[['label']] <- 'theta_-05'



ggplot(data = rbind(reshape2::melt(data.table(reshape2::melt(data_plot[, c('label', 'eff_1', 'eff_2')], 
							     id.vars = 'label'))[, .(mean = mean(value), med = median(value)), by = .(var = variable, label = label)],
				   id.vars = c('var', 'label')),
		    reshape2::melt(data.table(exp_Eff_mlt)[, .(mean = mean(value), med = median(value),
		    			   label = unique(label)), by = .(var = variable)],
		    			   id.vars = c('var', 'label'))),
       aes(x = -as.numeric(gsub('theta_', '',label)),
           y = value, group = label, fill = factor(variable, levels = c('med', 'mean'),
           					labels = c('Median', 'Mean'))))+ 
	geom_point(size = 3.5, shape = 23)+
	

	scale_y_continuous(TeX('Efficiency ($s^{-1}$)'), breaks = function(i){
		c(0, seq(0, round(i[2], 3 ), length.out = 6))
	}, limits = c(0, NA)) + theme(aspect.ratio = 0.75)+
	scale_x_continuous(TeX('Risk aversion ($\\Theta$)'), breaks = c(5, seq(-3, -30, -3)),
			   labels = function(i){
			   	c('Exp', TeX(paste0('$10^{',i, '}$'))[-1])
			   })+
	facet_wrap(~factor(var, labels = c('Exploration', 'Exploitation')),
		   ncol = 2, nrow = 1, scales = 'free')+
	geom_vline(linetype = 'dashed', xintercept =2.5)+
	scale_fill_manual('', values = c('darkred', 'dodgerblue3'))















# ggplot(data = rbind(reshape2::melt(data_plot, id.vars = c('label', '.id')),
# 		    exp_Eff_mlt),
#        aes(x = as.numeric(gsub('theta_', '',label)),
#            y = value, group = label))+
# 	# geom_boxplot(aes(x = as.numeric(gsub('theta_', '',label))),
# 	# 	     outlier.shape = NA, color = 'grey55')+
# 	# geom_jitter(alpha = 0.25, width = 1)+
# 	# geom_smooth(data =
# 	# 	    	reshape2::melt(global_Eff[label %in% c('theta_05','theta_10',
# 	# 	    					       'theta_15','theta_20', 'theta_25',
# 	# 	    					       'theta_30', 'theta_35', 'theta_40',
# 	# 	    					       'theta_45', 'theta_50', 'theta_60', 'theta_70',
# 	# 	    					       'theta_80', 'theta_90'),
# 	# 	    				  .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint)),
# 	# 	    				  by = 'label'], id.vars = c('label')),
# 	# 	    se = FALSE, method = 'gam', color = 'darkslateblue',
# 	# 	    formula = y ~ log(x), aes(group = 1))+
# 
# 	geom_smooth(data =
# 		    	reshape2::melt(data_plot[label %in% paste0('theta_',
# 		    						    sprintf('%02d', c(3:10, 15,20,25)))],
# 		    		       id.vars = c('label', '.id')),
# 		    se = FALSE, method = 'gam', color = 'darkslateblue',
# 		    formula = y ~ x/1-x, aes(group = 1))+
# 	# geom_label(data = label_data, aes(label = p, x = as.numeric(gsub('theta_', '',label))))+
# 	stat_summary(fun = 'mean', geom = 'point', shape = 23,
# 		     fill = 'dodgerblue3', color = 'black', size = 3.5)+
# 	scale_y_continuous(TeX('Efficiency ($s^{-1}$)'), breaks = function(i){
# 		seq(0, round(i[2], 3), length.out = 5)
# 	}) + theme(aspect.ratio = 1)+
# 	scale_x_continuous(TeX('Risk aversion ($\\Theta$)'), breaks = c(-5, seq(10, 100, 10)),
# 			   labels = function(i){
# 			   	c('Exp', TeX(paste0('$10^{-',i, '}$'))[-1])
# 			   })+
# 	facet_wrap(~factor(variable, labels = c('Exploration', 'Exploitation')),
# 		   ncol = 2, nrow = 1, scales = 'free')+
# 	geom_vline(linetype = 'dashed', xintercept =0)

