source('~/research/gits/AnTracks/src/Experiment.R')


##### CHECK DATES AND TIME STAMPS ####
library(lubridate)
load('~/research/gits/AnTracks/data/sto.RData')
l <- lapply(sto, function(i) lapply(range(rbindlist(i@food)[['t']])/2, seconds_to_period))
names(l) <- sapply(sto, function(i) i@date)


which.minN <- function(x, n = 2){
	len <- length(x)
	if(n > len) n <- len
	which(x == sort(x)[n]) 
}

data_dist <- rbindlist(lapply(1:9, function(i){
	f <- rbindlist(sto[[i]]@food)
	p1 <- ifelse(which.min(f[['t']]) %in% 1:6, 1, 2)
	p2 <- which(c(1,2) != p1)
	# data.table(d = as.numeric(pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(f[, c('x', 'y')]))),
	# 	   t = f[['t']], exp = i, patch = c(rep(p1, 6), rep(p2, 6)))
	data.table(d = as.numeric(pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(f[, c('x', 'y')]))),
		   t = f[['t']] - sto[[i]]@data[1, 'Frame'], exp = i, patch = c(rep(p1, 6), rep(p2, 6)))
}))

dist_between_patches <- unlist(lapply(1:9, function(i){
	f <- sto[[i]]@food
	p1 <- ifelse(which.min(rbindlist(f)[['t']]) %in% 1:6, 1, 2)
	p1_meanx <- mean(f[[p1]][['x']])
	p2_meanx <- mean(f[[-p1]][['x']])
	p1_meany <- mean(f[[p1]][['y']])
	p2_meany <- mean(f[[-p1]][['y']])
	as.numeric(pdist(matrix(c(p1_meanx, p1_meany), ncol = 2),
			 matrix(c(p2_meanx, p2_meany), ncol = 2)))
}))


######### EXPLORATION - EXPLOITATION DELAY ##########
ggplot(data = data_dist[, .(t = mean(t)[-which.min(t)] - t[which.min(t)]),
	  by = c('exp', 'patch')], aes(factor(patch), t/120)) + 
	geom_boxplot(outlier.shape = NA) + geom_jitter(size = 5, alpha = 0.8)+
	scale_y_continuous('<Patch time> - FPT (min)')+
	scale_x_discrete('', labels = c('First patch', 'Second patch'))


# R = 0.51, p = 2.2e-8
cor.test(x = data_dist[['t']], y = data_dist[['d']], method = 'pearson')
# R = 0.55, p = 5.7e-8
cor.test(x = data_dist[exp != 2, t], y = data_dist[exp != 2, d], method = 'pearson')
# rho = 0.66, p = 5.6e-15
cor.test(x = data_dist[, t], y = data_dist[, d], method = 'spearman')


######### PLOTS PATCH DISTANCE - FPT ##########
ggplot(data = data_dist, aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50))+
	scale_y_continuous('Distance to nest (cm)')+
	scale_fill_viridis_d('Exp')+coord_flip()


# R^2 around 0.25-0.3... not too relevant
# summary(lm(d ~ t, data = data_dist))
# summary(lm(d ~ t, data = data_dist[exp != 2]))
ggplot(data = data_dist[exp!=2], aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time (min)', breaks = seq(0, 60, 15))+
	scale_y_continuous('Distance to nest (cm)')+
	scale_fill_viridis_d('Exp')+coord_flip()



r2_TD_p1 <- round(summary(lm(t ~ d, data = data_dist[exp != 2 & patch == 1, .(d = d/10, t = t/120)]))$adj.r.squared, 3)
r2_TD_p2 <- round(summary(lm(t ~ d, data = data_dist[exp != 2 & patch == 2, .(d = d/10, t = t/120)]))$adj.r.squared, 3)

# summary(lm(d ~ t, data = data_dist[exp != 2 & patch == 1])) # R^2 ~ 0.4

data_plot_1 <- data_dist[exp!=2]
data_plot_1[['patch']] <- factor(data_plot_1[['patch']])
levels(data_plot_1[['patch']]) <- c(paste0('First patch (', '$R^{2}$',' = ',r2_TD_p1,')'),
				    paste0('Second patch (', '$R^{2}$',' = ',r2_TD_p2,')'))
ggplot(data =data_plot_1, aes(d/10, t/120)) + 
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black')+
	geom_smooth(method = 'lm', se = FALSE) +
	facet_wrap(~ patch, labeller = as_labeller(function(i){TeX(i)}, default = label_parsed))+
	scale_y_continuous('FPT (min)', breaks = seq(0, 60, 10))+
	scale_x_continuous('Distance to nest (cm)')+
	scale_fill_viridis_d('Experiment')+
	theme(legend.position = c(0.4, 0.8),
	      legend.background = element_rect(fill = NA, color = 'black'))

### WITHOUT THE OUTLIER 'EXP 2'
# ggplot(data = data_dist, aes(t/120, d / 10)) +
#         geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
#         scale_x_continuous('Time (min)', breaks = seq(0, 60, 15))+
#         scale_y_continuous('Distance to nest (cm)')+
#         scale_fill_viridis_d('Exp')+facet_wrap(~ patch)



######## CORRELATIONS ##########
# BOTH CORRELATIONS AND LMs GIVE NO SIGNIFICANT RESULTS IN EITHER ROW
# first row of plots
# cor.test(data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1])),
# 		   by = 'exp'][, d := dist_between_patches][['d']],
# 	 data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1])),
# 	 	  by = 'exp'][, d := dist_between_patches][['t']])
# 
# cor.test(data_dist[, .(t = min(t[patch == 2])),
# 		   by = 'exp'][, d := dist_between_patches][['d']],
# 	 data_dist[, .(t = min(t[patch == 2])),
# 	 	  by = 'exp'][, d := dist_between_patches][['t']])
# 
# cor.test(data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1])),
# 	  by = 'exp'][, d := dist_between_patches][['d']],
# data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1])),
# 	  by = 'exp'][, d := dist_between_patches][['t']])
# 
# 
# # second row of plots
# cor.test(data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
# 	      d = mean(d[patch == 2])),
# 	  by = 'exp'][['d']],
# data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
# 	      d = mean(d[patch == 2])),
# 	  by = 'exp'][['t']])
# 
# cor.test(data_dist[, .(t = min(t[patch == 2]),
# 	      d = mean(d[patch == 2])),
# 	  by = 'exp'][['d']],
# data_dist[, .(t = min(t[patch == 2]),
# 	      d = mean(d[patch == 2])),
# 	  by = 'exp'][['t']])
# 
# cor.test(data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1]),
# 	      d = mean(d[patch == 2])),
# 	  by = 'exp'][['d']],
# data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1]),
# 	      d = mean(d[patch == 2])),
# 	  by = 'exp'][['t']])

r2l_ <- round(summary(lm(t ~ d, data = data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
						     d = mean(d[patch == 2])),
						 by = 'exp'][, m := 1][exp != 3]))$adj.r.squared, 3)
r2m_ <- round(summary(lm(t ~ d, data = data_dist[, .(t = min(t[patch == 2]),
						     d = mean(d[patch == 2])),
						 by = 'exp'][, m := 2][exp != 2]))$adj.r.squared, 3)
r2r_ <- round(summary(lm(t ~ d, data = data_dist[, .(t = min(t[patch == 2]),
						     d = mean(d[patch == 2])),
						 by = 'exp'][, m := 2]))$adj.r.squared, 3)


ggarrange(
ggplot(data = rbind(
	data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1])),
		  by = 'exp'][, d := dist_between_patches][, m := 1],
	data_dist[, .(t = min(t[patch == 2])),
		  by = 'exp'][, d := dist_between_patches][, m := 2],
	data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1])),
		  by = 'exp'][, d := dist_between_patches][, m := 3]),
       aes(d / 10, t/120 )) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	geom_point(data = rbind(
		data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1])),
			  by = 'exp'][, d := dist_between_patches][, m := 1][exp == 3],
		data_dist[, .(t = min(t[patch == 2])),
			  by = 'exp'][, d := dist_between_patches][, m := 2][exp == 2],
		data_dist[, .(t = NA, d = NA),
			  by = 'exp'][, m := 3][exp == 1]),aes(fill = factor(exp)),
		shape = 4,size = 9, alpha = 0.9, color = 'red', show.legend = FALSE)+
	scale_y_continuous('Time (min)')+
	scale_x_continuous('Distance patch-patch (cm)', breaks = seq(0, 150, 25))+
	scale_fill_viridis_d('Experiment') +
	facet_wrap(~factor(m, labels = c('FPT2 - FPT1',
					 'FPT2',
					 'FPT2 - < FPT1 >')), ncol = 3, scales = 'free') + 
	# facet_wrap(~factor(m, labels = c('FPT2 - FPT1 (R^2 = 0.157)',
	# 				 'FPT2 (R^2 = 0.159)',
	# 				 'FPT2 - < FPT1 > (R^2 = -0.105)')), ncol = 3, scales = 'free') + 
	geom_smooth(data = rbind(
		data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1])),
			  by = 'exp'][, d := dist_between_patches][, m := 1][exp != 3],
		data_dist[, .(t = min(t[patch == 2])),
			  by = 'exp'][, d := dist_between_patches][, m := 2][exp != 2],
		data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1])),
			  by = 'exp'][, d := dist_between_patches][, m := 3]),
		method = 'lm', se = FALSE)
,

ggplot(data = rbind(
	data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
		      d = mean(d[patch == 2])),
		  by = 'exp'][, m := 1],
	data_dist[, .(t = min(t[patch == 2]),
		      d = mean(d[patch == 2])),
		  by = 'exp'][, m := 2],
	data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1]),
		      d = mean(d[patch == 2])),
		  by = 'exp'][, m := 3]),
       aes(d / 10, t/120 )) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	geom_point(data = rbind(
		data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
			      d = mean(d[patch == 2])),
			  by = 'exp'][, m := 1][exp == 3],
		data_dist[, .(t = min(t[patch == 2]),
			      d = mean(d[patch == 2])),
			  by = 'exp'][, m := 2][exp == 2],
		data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1]),
			      d = mean(d[patch == 2])),
			  by = 'exp'][, m := 3][exp == 3]),aes(fill = factor(exp)),
		shape = 4,size = 9, alpha = 0.9, color = 'red', show.legend = FALSE)+
	scale_y_continuous('Time (min)')+
	scale_x_continuous('Distance Nest-Patch 2 (cm)', breaks = seq(0, 150, 25))+
	scale_fill_viridis_d('Experiment') +
	facet_wrap(~factor(m, labels = c('FPT2 - FPT1 ',
					 'FPT2',
					 'FPT2 - < FPT1 >')), ncol = 3, scales = 'free') + 
	geom_smooth(data = rbind(
		data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
			      d = mean(d[patch == 2])),
			  by = 'exp'][, m := 1][exp != 3],
		data_dist[, .(t = min(t[patch == 2]),
			      d = mean(d[patch == 2])),
			  by = 'exp'][, m := 2][exp != 2],
		data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1]),
			      d = mean(d[patch == 2])),
			  by = 'exp'][, m := 3]),
		method = 'lm', se = FALSE)
, nrow = 2, ncol = 1, common.legend = T)



ggplot(data = data_dist[, .(fpt2 = min(t[patch == 2]), fpt1 = min(t[patch == 1]),
			    md = mean(d)), by = 'exp'][, dpatch := dist_between_patches],
       aes(fpt1/120, fpt2/120)) +
	geom_point(aes(size = dpatch/10, fill = md/10), shape = 21, alpha = 0.7) +
	scale_x_continuous('FPT patch 1 (min)', breaks = c(0, 5, 10, 15, 20, 115))+
	scale_y_continuous('FPT patch 2 (min)')+
	scale_size_continuous('Distance\nbetween\npatches (cm)', 
			      range = c(3, 15), breaks = c(25, 50, 75, 100))+
	scale_fill_viridis('Average\ndistance\npatch-nest (cm)')+
	theme(legend.position = c(0.7, 0.7),
	      legend.background = element_rect(fill = NA, color = 'black'))

# there seems not to be a correlation
ggplot(data = data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1])),
			    by = 'exp'][, d := dist_between_patches],
       aes(d / 10, t/120 )) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	# geom_smooth(data = data_dist[exp != 3, .(t = min(t[patch == 2]) - min(t[patch == 1])),
	# 		      by = 'exp'][, d := dist_between_patches[-3]],
	# 	    method = 'lm', se = FALSE, color = 'black')+
	scale_y_continuous('FPT2 - FPT1 (min)', breaks = seq(0, 60, 5))+
	scale_x_continuous('Distance patch-patch (cm)', breaks = seq(0, 150, 25))+
	scale_fill_viridis_d('Experiment') +
	theme(legend.position = c(0.8, 0.8),
	      legend.background = element_rect(fill = NA, color = 'black'))

ggplot(data = data_dist[, .(t = min(t[patch == 2])),
			by = 'exp'][, d := dist_between_patches],
       aes(d / 10, t/120)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_y_continuous('FPT2 (min)', breaks = c(seq(0, 60, 10), 120))+
	scale_x_continuous('Distance patch-patch (cm)')+
	scale_fill_viridis_d('Experiment')+
	theme(legend.position = c(0.8, 0.8),
	      legend.background = element_rect(fill = NA, color = 'black'))


ggplot(data = data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1])),
			by = 'exp'][, d := dist_between_patches],
       aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('FPT2 - < FPT1 > (min)', breaks = seq(0, 60, 5))+
	scale_y_continuous('Distance patch-patch (cm)', seq(0, 150, 25))+
	scale_fill_viridis_d('Experiment') +
coord_flip()+
	theme(legend.position = c(0.8, 0.8),
	      legend.background = element_rect(fill = NA, color = 'black'))


ggplot(data = data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
			    d = mean(d[patch == 2])),
			by = 'exp'],
       aes(d / 10, t/120)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	geom_point(data = data_dist[exp == 3, .(t = min(t[patch == 2]) - min(t[patch == 1]),
					d = mean(d[patch == 2])),
					by = 'exp'], shape = 4, size = 10,color = 'red')+
	geom_smooth(data = data_dist[exp != 3, .(t = min(t[patch == 2]) - min(t[patch == 1]),
						 d = mean(d[patch == 2])),
						 by = 'exp'],
		    method = 'lm', color = 'black', se = FALSE)+
	scale_y_continuous('Time difference between patches (min)', breaks = seq(0, 60, 10))+
	scale_x_continuous('Second patch-Nest distance (cm)')+
	scale_fill_viridis_d('Experiment') +
	theme(legend.position = c(0.8, 0.8),
	      legend.background = element_rect(fill = NA, color = 'black'))

ggplot(data = data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1]),
			    d = mean(d[patch == 2])),
			by = 'exp'],
       aes(d / 10,t/120)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_y_continuous('FPT2 - <Patch 1 times> (min)', breaks = seq(0, 60, 5))+
	scale_x_continuous('Second patch-Nest distance (cm)')+
	scale_fill_viridis_d('Experiment')+
	theme(legend.position = c(0.3, 0.8),
	      legend.background = element_rect(fill = NA, color = 'black'))

ggplot(data = data_dist[exp != 2, .(t = min(t[patch == 2]),
			    d = mean(d[patch == 2])),
			by = 'exp'],
       aes(d / 10, t/120)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_y_continuous('FPT second patch (min)', breaks = seq(0, 60, 10))+
	scale_x_continuous('Nest-Second patch distance (cm)')+
	scale_fill_viridis_d('Exp') # + coord_cartesian(xlim = c(0,20), ylim = c(0,20))


ggplot(data = melt(data_dist[exp != 2, .(tdiff = min(t[patch == 2]) - min(t[patch == 1]),
					 tp1 = min(t[patch == 1]),
					 d1 = d[which.min(t)], d2 = d[which.min(t[patch == 2])]), by = 'exp'][, c('exp', 'tdiff', 'tp1')],
		   id.vars = 'exp'),
       aes(variable, value / 120)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 5, alpha = 0.5)+
	ylab('Time (min)') + scale_x_discrete('', labels = c('Time difference', 'TP1'))

######### PROBS TIME ##########
dt <- rbind(
	data.frame(t = c('01:19','02:52','07:55', '08:04', '08:13', '08:35','09:16','09:30','12:03','12:19', '15:11',
			 '12:39', '15:43', '17:39', '19:40', '19:54', '20:40'),
			 p = c(0,0,1,0,1,1,0,1,0,1, 1,
			       1, 1, 1, 1, 1, 1),
			 exp = 1
	),
	data.frame(t = c('15:38', '25:13', '29:13', '32:45', '33:50', '34:57',
			 '03:26', '06:09', '09:17', '09:54', '13:53', '15:28'),
			 p = c(1,1 ,1,1, 1,1, 1,1, 1,1,1, 1),
			 exp = 2),
	data.frame(t = c('01:46', '32:55', '35:13', '36:12', '36:41', '37:12',
			 '39:37','41:16', '42:22', '43:03', '43:34', '45:04'),
			 p = c(1, 1, 1,1,1,1,1,1, 1, 1, 1,1),
			 exp = 3),
	data.frame(t = c('00:00', '02:29','4:29', '5:34', '05:37', '06:35', '07:10','07:45', '07:58', '08:05', '08:19',
			 '01:41','11:23', '11:48', '14:16', '15:26', '17:46'),
			 p = c(1,1, 0,0, 0,0,1, 1,0, 1, 1,1,1, 1,1,1, 1),
			 exp = 4),
	data.frame(t = c('02:55', '14:38', '15:01', '17:12', '15:20', '15:30', '16:01',
			 '12:07', '13:49','16:45', '17:07', '17:39', '18:10'),
			 p = c(1,1, 1, 1,0, 1, 1, 1, 1, 1, 1, 1, 1),
			 exp = 5),
	data.frame(t = c('13:58', '17:56', '19:54', '20:07', '20:24','25:52','31:20',
			 '01:29', '04:45', '04:49', '05:41', '08:17', '10:16', '10:46', '11:40'),
			 p = c(1, 1, 0, 1, 1,1,1, 1, 1, 1, 0, 1, 1, 0, 1), # ONE IS WRONG
			 exp = 6),
	data.frame(t = c('14:51', '16:32', '17:49', '18:57', '20:14', '21:28',
			 '01:21', '04:16', '04:45', '05:21', '05:50', '06:25', '07:17', '07:26', '07:39',
			 '08:35', '08:43', '09:07', '10:00'),
			 p = c(1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0,0, 0, 1),
			 exp = 7),
	data.frame(t = c('03:18', '10:20', '12:45', '15:17', '15:20', '15:30','18:10',
			 '18:33','18:42','20:08','20:45','22:26','22:50', '23:00'),
			 p = c(1,1,1,1,1,0,1,1,1,1,1, 0,1, 1),
			 exp = 8),
	data.frame(t = c('17:14', '23:06', '24:31','24:34', '27:08', '27:59',
			 '04:14', '08:13', '12:57', '14:26', '14:44', '14:54'),
			 p = c(1,1,1,1,1,1, 1, 1, 1, 1, 1, 1),
			 exp = 9))


mean(dt[['p']])
# sum(dt[['p']]) # CHECK; sums 12 * num of exps (108)