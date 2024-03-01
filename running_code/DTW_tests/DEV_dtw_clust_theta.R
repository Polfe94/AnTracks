#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/Experiment.R')

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

prop.clusters <- function(distList, pattern = 'g_', label = NULL){
	a <- t(sapply(distList, function(i){
		table(factor(knames[apply(i, 1, which.min)],
			     levels = c('norm', 'late', 'low', 'flat')))
	}))
	a <- data.table(cbind(g = as.numeric(gsub(pattern, '', rownames(a))), a))
	if(!is.null(label)){
		a[['label']] <- label
	}
	a
}

meand.clusters <- function(distList, pattern = 'g_', label = NULL){
	a <- t(sapply(distList, function(i){
		table(factor(knames[apply(i, 1, which.min)],
			     levels = c('norm', 'late', 'low', 'flat')))
	}))
	a <- data.table(cbind(g = as.numeric(gsub(pattern, '', rownames(a))), a))
	if(!is.null(label)){
		a[['label']] <- label
	}
	a
}


knames <- c('norm', 'late', 'flat', 'low')

########## THETA 05 
path_theta05 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gains_thetas/theta_05/'
theta05 <- lapply(list.files(path_theta05), function(i) process_data(paste0(path_theta05,i,'/'))) # takes ~20 sec
dTheta05 <- lapply(theta05, function(k){
	a <- rbindlist(lapply(k, function(i){
		data.table(t(sapply(zSeries, function(x){
			ndtw(i[['Z']], x)
		})))
	}))
	colnames(a) <- knames
	a
})
names(dTheta05) <- paste0('g_0.', 1:9)
kTheta05 <- prop.clusters(dTheta05, pattern = 'g_', label = 'theta_05')

########## THETA 10 
path_theta10 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gains_thetas/theta_10/'
theta10 <- lapply(list.files(path_theta10), function(i) process_data(paste0(path_theta10,i,'/'))) # takes ~20 sec
dTheta10 <- lapply(theta10, function(k){
	a <- rbindlist(lapply(k, function(i){
		data.table(t(sapply(zSeries, function(x){
			ndtw(i[['Z']], x)
		})))
	}))
	colnames(a) <- knames
	a
})
names(dTheta10) <- paste0('g_0.', 1:9)
kTheta10 <- prop.clusters(dTheta10, pattern = 'g_', label = 'theta_10')

########## THETA 15 
path_theta15 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gains_thetas/theta_15/'
theta15 <- lapply(list.files(path_theta15), function(i) process_data(paste0(path_theta15,i,'/'))) # takes ~20 sec
dTheta15 <- lapply(theta15, function(k){
	a <- rbindlist(lapply(k, function(i){
		data.table(t(sapply(zSeries, function(x){
			ndtw(i[['Z']], x)
		})))
	}))
	colnames(a) <- knames
	a
})
names(dTheta15) <- paste0('g_0.', 1:9)
kTheta15 <- prop.clusters(dTheta15, pattern = 'g_', label = 'theta_15')

########## THETA 20 
path_theta20 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gains_thetas/theta_20/'
theta20 <- lapply(list.files(path_theta20), function(i) process_data(paste0(path_theta20,i,'/'))) # takes ~20 sec
dTheta20 <- lapply(theta20, function(k){
	a <- rbindlist(lapply(k, function(i){
		data.table(t(sapply(zSeries, function(x){
			ndtw(i[['Z']], x)
		})))
	}))
	colnames(a) <- knames
	a
})
names(dTheta20) <- paste0('g_0.', 1:9)
kTheta20 <- prop.clusters(dTheta20, pattern = 'g_', label = 'theta_20')

########## THETA 25 
path_theta25 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gains_thetas/theta_25/'
theta25 <- lapply(list.files(path_theta25), function(i) process_data(paste0(path_theta25,i,'/'))) # takes ~20 sec
dTheta25 <- lapply(theta25, function(k){
	a <- rbindlist(lapply(k, function(i){
		data.table(t(sapply(zSeries, function(x){
			ndtw(i[['Z']], x)
		})))
	}))
	colnames(a) <- knames
	a
})
names(dTheta25) <- paste0('g_0.', 1:9)
kTheta25 <- prop.clusters(dTheta25, pattern = 'g_', label = 'theta_25')

########## THETA 30 
path_theta30 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gains_thetas/theta_30/'
theta30 <- lapply(list.files(path_theta30), function(i) process_data(paste0(path_theta30,i,'/'))) # takes ~20 sec
dTheta30 <- lapply(theta30, function(k){
	a <- rbindlist(lapply(k, function(i){
		data.table(t(sapply(zSeries, function(x){
			ndtw(i[['Z']], x)
		})))
	}))
	colnames(a) <- knames
	a
})
names(dTheta30) <- paste0('g_0.', 1:9)
kTheta30 <- prop.clusters(dTheta30, pattern = 'g_', label = 'theta_30')


data_plot <- rbind(kTheta05, kTheta10, kTheta15, kTheta20, kTheta25, kTheta30)
data_plot$theta <- sapply((parse(text = paste0(gsub('theta_', '', data_plot$label)))), eval)

hmap_gTheta_Pexp <- data_plot
save(hmap_gTheta_Pexp, file = '/home/polfer/research/gits/AutomatAnts/data/hmap_gains_theta/hmap_gTheta_Pexp.RData')

ggplot(data = data_plot)+
	geom_tile(aes(g, theta,
		  fill = norm / 100))+
	scale_fill_viridis_c(TeX('$P_{exp}$'), limits = c(0, 1))+
	scale_y_continuous(TeX('Risk aversion ($\\Theta$)'),
			   breaks = seq(5, 30, 5), 
			   labels = c(10**-5,10**-10,10**-15,10**-20,10**-25,10**-30))+
	scale_x_continuous(TeX('Sensitivity to information ($g$)'), breaks = seq(0.1, 0.9, 0.2))


png('/home/polfer/research/gits/AnTracks/plots/g_theta_hmap.png', 4000, 4000, res = 400)
ggplot(data = reshape2::melt(data_plot, id.vars = c('g', 'label', 'theta')))+
	geom_tile(aes(g, theta,
		      fill = value / 100))+
	scale_fill_viridis_c(TeX('$Probability$'), limits = c(0, 1))+
	scale_y_continuous(TeX('Risk aversion ($\\Theta$)'),
			   breaks = seq(5, 30, 5), 
			   labels = function(i) TeX(paste0('$10^{-',i,'}')))+
	scale_x_continuous(TeX('Sensitivity to information ($g$)'), breaks = seq(0.1, 0.9, 0.2))+
	facet_wrap(~factor(variable, levels = c('norm', 'flat', 'late', 'low'),
			   labels = c('Experimental', 'Plateau', 'Late-start', 'Low activity')), 
		   nrow = 2, ncol = 2)+
	theme(aspect.ratio = 0.75)

dev.off()
# ggplot(data = data_plot)+
# 	geom_tile(aes(g, theta,
# 		      fill = flat / 100))+
# 	scale_fill_viridis_c(TeX('$P_{plat}$'), limits = c(0, 1))+
# 	scale_y_continuous(TeX('Risk aversion ($\\Theta$)'),
# 			   breaks = seq(5, 30, 5), 
# 			   labels = c(10**-5,10**-10,10**-15,10**-20,10**-25,10**-30))+
# 	scale_x_continuous(TeX('Sensitivity to information ($g$)'), breaks = seq(0.1, 0.9, 0.2))
# 
# ggplot(data = data_plot)+
# 	geom_tile(aes(g, theta,
# 		      fill = low / 100))+
# 	scale_fill_viridis_c(TeX('$P_{low}$'), limits = c(0, 1))+
# 	scale_y_continuous(TeX('Risk aversion ($\\Theta$)'),
# 			   breaks = seq(5, 30, 5), 
# 			   labels = c(10**-5,10**-10,10**-15,10**-20,10**-25,10**-30))+
# 	scale_x_continuous(TeX('Sensitivity to information ($g$)'), breaks = seq(0.1, 0.9, 0.2))
# 
# ggplot(data = data_plot)+
# 	geom_tile(aes(g, theta,
# 		      fill = late / 100))+
# 	scale_fill_viridis_c(TeX('$P_{late}$'), limits = c(0, 1))+
# 	scale_y_continuous(TeX('Risk aversion ($\\Theta$)'),
# 			   breaks = seq(5, 30, 5), 
# 			   labels = c(10**-5,10**-10,10**-15,10**-20,10**-25,10**-30))+
# 	scale_x_continuous(TeX('Sensitivity to information ($g$)'), breaks = seq(0.1, 0.9, 0.2))




tancurves <- data.table(g_0.1 = c(1, vapply(1:50, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x)};x}, numeric(1))),
			g_0.5 = c(1, vapply(1:50, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x)};x}, numeric(1))),
			g_0.9 = c(1, vapply(1:50, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x)};x}, numeric(1))),
			iter = 0:50, theta = 0)

tantheta_05 <- data.table(g_0.1 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-5))};x}, numeric(1))),
			  g_0.5 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-5))};x}, numeric(1))),
			  g_0.9 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-5))};x}, numeric(1))),
			  iter = 0:250, theta = 5)

tantheta_15 <- data.table(g_0.1 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-15))};x}, numeric(1))),
			  g_0.5 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-15))};x}, numeric(1))),
			  g_0.9 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-15))};x}, numeric(1))),
			  iter = 0:250, theta = 15)

tantheta_25 <- data.table(g_0.1 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-25))};x}, numeric(1))),
			  g_0.5 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-25))};x}, numeric(1))),
			  g_0.9 = c(1, vapply(1:250, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-25))};x}, numeric(1))),
			  iter = 0:250, theta = 25)

tancurves <- data.table(g_0.1 = vapply(1:100, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x)};x}, numeric(1)),
			g_0.5 = vapply(1:100, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x)};x}, numeric(1)),
			g_0.9 = vapply(1:100, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x)};x}, numeric(1)),
			iter = 1:100)




ggplot(data = reshape2::melt(rbind(
	data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x)};x}, numeric(1))),
		   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x)};x}, numeric(1))),
		   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x)};x}, numeric(1))),
		   g_1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x)};x}, numeric(1))),
		   iter = 0:500, theta = 0),
	data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-5))};x}, numeric(1))),
		   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-5))};x}, numeric(1))),
		   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-5))};x}, numeric(1))),
		   g_1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x-(10**-5))};x}, numeric(1))),
		   iter = 0:500, theta = 5),
	data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-15))};x}, numeric(1))),
		   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-15))};x}, numeric(1))),
		   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-15))};x}, numeric(1))),
		   g_1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x-(10**-15))};x}, numeric(1))),
		   iter = 0:500, theta = 15),
	data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-25))};x}, numeric(1))),
		   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-25))};x}, numeric(1))),
		   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-25))};x}, numeric(1))),
		   g_1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x-(10**-25))};x}, numeric(1))),
		   iter = 0:500, theta = 25)
), id.vars = c('iter', 'theta')), aes(iter, value, color = variable)) +
	geom_path(linewidth = 1) +
	scale_color_viridis_d('Sensitivity (g)', begin = 0.1, end = 0.8,
			      labels = c(0.1, 0.5, 0.9, 1))+
	scale_x_continuous('Iteration') +
	scale_y_continuous(TeX('$S_{i}(t)$'))+
	facet_wrap(~theta)+
	geom_vline(data = data.table())





bigdf <- rbind(
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x)};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x)};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x)};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x)};x}, numeric(1)),
		   iter = 1:1000, theta = 0),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-5))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-5))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-5))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x-(10**-5))};x}, numeric(1)),
		   iter = 1:1000, theta = 5),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-15))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-15))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-15))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x-(10**-15))};x}, numeric(1)),
		   iter = 1:1000, theta = 15),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-25))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-25))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-25))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x-(10**-25))};x}, numeric(1)),
		   iter = 1:1000, theta = 25))

bigdf_high_rate <- rbind(
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.1*(runif(1) + x)),tanh(0.1*x))
		};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.5*(runif(1)+x)),tanh(0.5*x))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.9*(runif(1)+x)),tanh(0.9*x))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(1*(runif(1)+x)),tanh(x))};x}, numeric(1)),
		   iter = 1:1000, theta = 0),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.1*(runif(1)+x-10**-5)),tanh(0.1*x-10**-5))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.5*(runif(1)+x-10**-5)),tanh(0.5*x-10**-5))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.9*(runif(1)+x-10**-5)),tanh(0.9*x-10**-5))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(1*(runif(1)+x-10**-5)),tanh(1*x-10**-5))};x}, numeric(1)),
		   iter = 1:1000, theta = 5),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.1*(runif(1)+x-10**-15)),tanh(0.1*x-10**-15))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.5*(runif(1)+x-10**-15)),tanh(0.5*x-10**-15))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.9*(runif(1)+x-10**-15)),tanh(0.9*x-10**-15))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(1*(runif(1)+x-10**-15)),tanh(1*x-10**-15))};x}, numeric(1)),
		   iter = 1:1000, theta = 15),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.1*(runif(1)+x-10**-25)),tanh(0.1*x-10**-25))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.5*(runif(1)+x-10**-25)),tanh(0.5*x-10**-25))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(0.9*(runif(1)+x-10**-25)),tanh(0.9*x-10**-25))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%10 == 0, tanh(1*(runif(1)+x-10**-25)),tanh(1*x-10**-25))};x}, numeric(1)),
		   iter = 1:1000, theta = 25))


bigdf_low_rate <- rbind(
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.1*(runif(1) + x)),tanh(0.1*x))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.5*(runif(1)+x)),tanh(0.5*x))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.9*(runif(1)+x)),tanh(0.9*x))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(1*(runif(1)+x)),tanh(1*x))};x}, numeric(1)),
		   iter = 1:1000, theta = 0),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.1*(runif(1)+x-10**-5)),tanh(0.1*x-10**-5))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.5*(runif(1)+x-10**-5)),tanh(0.5*x-10**-5))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.9*(runif(1)+x-10**-5)),tanh(0.9*x-10**-5))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(1*(runif(1)+x-10**-5)),tanh(1*x-10**-5))};x}, numeric(1)),
		   iter = 1:1000, theta = 5),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.1*(runif(1)+x-10**-15)),tanh(0.1*x-10**-15))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.5*(runif(1)+x-10**-15)),tanh(0.5*x-10**-15))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.9*(runif(1)+x-10**-15)),tanh(0.9*x-10**-15))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(1*(runif(1)+x-10**-15)),tanh(1*x-10**-15))};x}, numeric(1)),
		   iter = 1:1000, theta = 15),
	data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.1*(runif(1)+x-10**-25)),tanh(0.1*x-10**-25))};x}, numeric(1)),
		   g_0.5 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.5*(runif(1)+x-10**-25)),tanh(0.5*x-10**-25))};x}, numeric(1)),
		   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(0.9*(runif(1)+x-10**-25)),tanh(0.9*x-10**-25))};x}, numeric(1)),
		   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){if(is.na(x) | x < 0)x <- NA;x <- ifelse(i%%100 == 0, tanh(1*(runif(1)+x-10**-25)),tanh(1*x-10**-25))};x}, numeric(1)),
		   iter = 1:1000, theta = 25))


cutdf <- rbind(bigdf[, .(v = g_0.1[g_0.1 >= 0], g = 'g_0.1', iter = iter[g_0.1 >= 0]), by = 'theta'],
      bigdf[, .(v= g_0.5[g_0.5 >= 0], g = 'g_0.5', iter = iter[g_0.5 >= 0]), by = 'theta'],
      bigdf[, .(v = g_0.9[g_0.9 >= 0], g = 'g_0.9', iter = iter[g_0.9 >= 0]), by = 'theta'],
      bigdf[, .(v = g_1[g_1 >= 0], g = 'g_1', iter = iter[g_1 >= 0]), by = 'theta'])

cutdf_low_rate <- rbind(bigdf_low_rate[, .(v = g_0.1[g_0.1 >= 0], g = 'g_0.1', iter = iter[g_0.1 >= 0]), by = 'theta'],
			 bigdf_low_rate[, .(v= g_0.5[g_0.5 >= 0], g = 'g_0.5', iter = iter[g_0.5 >= 0]), by = 'theta'],
			 bigdf_low_rate[, .(v = g_0.9[g_0.9 >= 0], g = 'g_0.9', iter = iter[g_0.9 >= 0]), by = 'theta'],
			 bigdf_low_rate[, .(v = g_1[g_1 >= 0], g = 'g_1', iter = iter[g_1 >= 0]), by = 'theta'])


cutdf_high_rate <- rbind(bigdf_high_rate[, .(v = g_0.1[g_0.1 >= 0], g = 'g_0.1', iter = iter[g_0.1 >= 0]), by = 'theta'],
	       bigdf_high_rate[, .(v= g_0.5[g_0.5 >= 0], g = 'g_0.5', iter = iter[g_0.5 >= 0]), by = 'theta'],
	       bigdf_high_rate[, .(v = g_0.9[g_0.9 >= 0], g = 'g_0.9', iter = iter[g_0.9 >= 0]), by = 'theta'],
	       bigdf_high_rate[, .(v = g_1[g_1 >= 0], g = 'g_1', iter = iter[g_1 >= 0]), by = 'theta'])



png('/home/polfer/research/gits/AnTracks/plots/mega_plot_Si.png', height = 6000, width = 4000, res = 400)
ggarrange(
ggplot(data = cutdf) +
	geom_path(linewidth = 1, 
		  aes(iter, v, color = factor(g))) +
	scale_color_viridis_d('Sensitivity (g)', begin = 0, end = 0.95,
			      labels = c(0.1, 0.5, 0.9, 1))+
	scale_y_continuous(TeX('$S_{i}(t)$'), limits = c(NA, 1))+
	
scale_x_log10('Iteration') + facet_wrap(~factor(theta), labeller = as_labeller(function(i){
	c(TeX('$\\Theta = 0$'), TeX(paste0('$\\Theta = 10^{-',i,'}$')[-1]))
}, default = label_parsed)) +
	annotation_logticks(sides = 'b') + theme_bw() +
	geom_linerange(data = data.table::melt(bigdf[, .(g_0.1 = which.max(g_0.1 < 0)-1, 
					      g_0.5 = which.max(g_0.5 <0)-1,
					      g_0.9 = which.max(g_0.9 < 0)-1,
					      g_1 = which.max(g_1 < 0)-1), by = 'theta'], 
				    id.vars = 'theta', 
				    measure.vars = patterns('^g_'), value.name = 'iter')[iter == 0,
				    						     iter := Inf][!is.na(iter)],
		    linewidth = 1.5, aes(x = iter, ymin = -0.0175, ymax = 0.0175, color = variable),
		    show.legend = FALSE)+
	theme(axis.title = element_text(size  = 15, color = 'black'),
							    axis.text = element_text(size = 15, color = 'black'),
							    legend.text = element_text(size =15, color = 'black'),
							    legend.title = element_text(size = 15, color = 'black'),
							    strip.text = element_text(size = 15, margin = margin(t = 5, b = 5)),
	      aspect.ratio = 0.5),
ggplot(data = cutdf_high_rate) +
	geom_path(linewidth = 1, 
		  aes(iter, v, color = factor(g))) +
	scale_color_viridis_d('Sensitivity (g)', begin = 0, end = 0.95,
			      labels = c(0.1, 0.5, 0.9, 1))+
	scale_y_continuous(TeX('$S_{i}(t)$'))+
	
	scale_x_log10('Iteration') + facet_wrap(~factor(theta), labeller = as_labeller(function(i){
		c(TeX('$\\Theta = 0$'), TeX(paste0('$\\Theta = 10^{-',i,'}$')[-1]))
	}, default = label_parsed)) +
	annotation_logticks(sides = 'b') + theme_bw() +
	geom_linerange(data = data.table::melt(bigdf_high_rate[, .(g_0.1 = which.max(g_0.1 < 0)-1, 
							 g_0.5 = which.max(g_0.5 <0)-1,
							 g_0.9 = which.max(g_0.9 < 0)-1,
							 g_1 = which.max(g_1 < 0)-1), by = 'theta'], 
							 id.vars = 'theta', 
							 measure.vars = patterns('^g_'), value.name = 'iter')[iter == 0,
							 						     iter := Inf][!is.na(iter)],
							 linewidth = 1.5, aes(x = iter, ymin = -0.0175, ymax = 0.0175, color = variable),
							 show.legend = FALSE)+
	theme(axis.title = element_text(size  = 15, color = 'black'),
	      axis.text = element_text(size = 15, color = 'black'),
	      legend.text = element_text(size =15, color = 'black'),
	      legend.title = element_text(size = 15, color = 'black'),
	      strip.text = element_text(size = 15, margin = margin(t = 5, b = 5)),
	      aspect.ratio = 0.5),
ggplot(data = cutdf_low_rate) +
	geom_path(linewidth = 1, 
		  aes(iter, v, color = factor(g))) +
	scale_color_viridis_d('Sensitivity (g)', begin = 0, end = 0.95,
			      labels = c(0.1, 0.5, 0.9, 1))+
	scale_y_continuous(TeX('$S_{i}(t)$'))+
	
	scale_x_log10('Iteration') + facet_wrap(~factor(theta), labeller = as_labeller(function(i){
		c(TeX('$\\Theta = 0$'), TeX(paste0('$\\Theta = 10^{-',i,'}$')[-1]))
	}, default = label_parsed)) +
	annotation_logticks(sides = 'b') + theme_bw() +
	geom_linerange(data = data.table::melt(bigdf_low_rate[, .(g_0.1 = which.max(g_0.1 < 0)-1, 
							 g_0.5 = which.max(g_0.5 <0)-1,
							 g_0.9 = which.max(g_0.9 < 0)-1,
							 g_1 = which.max(g_1 < 0)-1), by = 'theta'], 
							 id.vars = 'theta', 
							 measure.vars = patterns('^g_'), value.name = 'iter')[iter == 0,
							 						     iter := Inf][!is.na(iter)],
							 linewidth = 1.5, aes(x = iter, ymin = -0.0175, ymax = 0.0175, color = variable),
							 show.legend = FALSE)+
	theme(axis.title = element_text(size  = 15, color = 'black'),
	      axis.text = element_text(size = 15, color = 'black'),
	      legend.text = element_text(size =15, color = 'black'),
	      legend.title = element_text(size = 15, color = 'black'),
	      strip.text = element_text(size = 15, margin = margin(t = 5, b = 5)),
	      aspect.ratio = 0.5)
, labels = c('a', 'b', 'c'), ncol = 1, font.label = list(size = 18, font = 'plain'))
dev.off()	



##############3
ggplot(data = 
reshape2::melt(data.table(g_0.1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x)};x}, numeric(1)),
	   g_0.5 = vapply(1:100, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x)};x}, numeric(1)),
	   g_0.9 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x)};x}, numeric(1)),
	   g_1 = vapply(1:1000, function(i){x <- 1;for(z in 1:i){x <- tanh(1*x)};x}, numeric(1)),
	   iter = 1:1000, theta = 0), id.vars = c('iter', 'theta'))
) + geom_path(aes(iter, value, color = factor(variable)))



rbind(
data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x)};x}, numeric(1))),
	   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x)};x}, numeric(1))),
	   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x)};x}, numeric(1))),
	   iter = 0:50, theta = 0),
data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-5))};x}, numeric(1))),
	   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-5))};x}, numeric(1))),
	   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-5))};x}, numeric(1))),
	   iter = 0:250, theta = 5),
data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-15))};x}, numeric(1))),
	   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-15))};x}, numeric(1))),
	   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-15))};x}, numeric(1))),
	   iter = 0:250, theta = 15),
data.table(g_0.1 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.1*x-(10**-25))};x}, numeric(1))),
	   g_0.5 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.5*x-(10**-25))};x}, numeric(1))),
	   g_0.9 = c(1, vapply(1:500, function(i){x <- 1;for(z in 1:i){x <- tanh(0.9*x-(10**-25))};x}, numeric(1))),
	   iter = 0:250, theta = 25)
)

ggplot(data = reshape2::melt(tancurves, id.vars = 'iter'), aes(iter, value, color = variable)) +
	geom_path(linewidth = 1) +
	scale_color_viridis_d('Sensitivity (g)', begin = 0.1, end = 0.8,
			      labels = c(0.1, 0.5, 0.9))+
	scale_x_continuous('Iteration') +
	scale_y_continuous(TeX('$S_{i}(t)$'))
	


ggplot(data = reshape2::melt(tantheta_05, id.vars = 'iter'), aes(iter, value, color = variable)) +
	geom_path(linewidth = 1) +
	scale_color_viridis_d('Sensitivity (g)', begin = 0.1, end = 0.8,
			      labels = c(0.1, 0.5, 0.9))+
	scale_x_continuous('Iteration') +
	scale_y_continuous(TeX('$S_{i}(t)$'))
