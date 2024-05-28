source('~/research/gits/AnTracks/src/Experiment.R')


### CHECK DATES AND TIME STAMPS
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
	data.table(d = as.numeric(pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(f[, c('x', 'y')]))),
		   t = f[['t']], exp = i, patch = c(rep(p1, 6), rep(p2, 6)))
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

# plot(data_dist[, .(t = t[which.minN(t, n= 2)] - t[which.min(t)],
# 	      d = mean(c(d[which.minN(t, n= 2)], d[which.min(t)]))), by = c('exp', 'patch')][, c('t', 'd')])
# 
# plot(data_dist[, .(t = mean(t)[-which.min(t)] - t[which.min(t)],
# 		   d = mean(c(d[which.minN(t, n= 2)], d[which.min(t)]))), by = c('exp', 'patch')][, c('t', 'd')])


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



ggplot(data = data_dist, aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50))+
	scale_y_continuous('Distance to nest (cm)')+
	scale_fill_viridis_d('Exp')


# R^2 around 0.25-0.3... not too relevant
# summary(lm(d ~ t, data = data_dist))
# summary(lm(d ~ t, data = data_dist[exp != 2]))
ggplot(data = data_dist[exp!=2], aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time (min)', breaks = seq(0, 60, 15))+
	scale_y_continuous('Distance to nest (cm)')+
	scale_fill_viridis_d('Exp')

# summary(lm(d ~ t, data = data_dist[exp != 2 & patch == 1])) # R^2 ~ 0.4
ggplot(data = data_dist[exp!=2], aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time (min)', breaks = seq(0, 60, 15))+
	scale_y_continuous('Distance to nest (cm)')+
	scale_fill_viridis_d('Exp')+facet_wrap(~ patch)


### WITHOUT THE OUTLIER 'EXP 2'
# ggplot(data = data_dist, aes(t/120, d / 10)) +
#         geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
#         scale_x_continuous('Time (min)', breaks = seq(0, 60, 15))+
#         scale_y_continuous('Distance to nest (cm)')+
#         scale_fill_viridis_d('Exp')+facet_wrap(~ patch)


ggplot(data = data_dist[, .(fpt2 = min(t[patch == 2]), fpt1 = min(t[patch == 1]),
			    md = mean(d)), by = 'exp'][, dpatch := dist_between_patches],
       aes(fpt1/120, fpt2/120)) +
	geom_point(aes(size = dpatch/10, fill = md/10), shape = 21, alpha = 0.7) +
	scale_x_continuous('FPT patch 1 (min)', breaks = c(5, 10, 15, 20, 25, 115))+
	scale_y_continuous('FPT patch 2 (min)')+
	scale_size_continuous('Mean\ndistance\nbetween\npatches (cm)', 
			      range = c(3, 15), breaks = c(25, 50, 75, 100))+
	scale_fill_continuous('Mean\ndistance\npatch-nest (cm)')
	      strip.text = element_blank())

# there seems not to be a correlation
ggplot(data = data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1])),
			    by = 'exp'][, d := dist_between_patches],
       aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time difference between patches (min)', breaks = seq(0, 60, 10))+
	scale_y_continuous('Mean distance patch-patch (cm)')+
	scale_fill_viridis_d('Exp') 

ggplot(data = data_dist[, .(t = min(t[patch == 2])),
			by = 'exp'][, d := dist_between_patches],
       aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time difference between patches (min)', breaks = seq(0, 60, 10))+
	scale_y_continuous('Mean distance patch-patch (cm)')+
	scale_fill_viridis_d('Exp')

ggplot(data = data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1])),
			by = 'exp'][, d := dist_between_patches],
       aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time difference between patches (min)', breaks = seq(0, 60, 10))+
	scale_y_continuous('Mean distance patch-patch (cm)')+
	scale_fill_viridis_d('Exp') 
coord_flip()

ggplot(data = data_dist[, .(t = min(t[patch == 2]) - min(t[patch == 1]),
			    d = mean(d[patch == 2])),
			by = 'exp'],
       aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time difference between patches (min)', breaks = seq(0, 60, 10))+
	scale_y_continuous('Mean distance nest-patch (cm)')+
	scale_fill_viridis_d('Exp') # + coord_cartesian(xlim = c(0,20), ylim = c(0,20))

ggplot(data = data_dist[, .(t = min(t[patch == 2]) - mean(t[patch == 1]),
			    d = mean(d[patch == 2])),
			by = 'exp'],
       aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time difference between patches (min)', breaks = seq(0, 60, 10))+
	scale_y_continuous('Mean distance nest-patch (cm)')+
	scale_fill_viridis_d('Exp')+coord_flip()

ggplot(data = data_dist[exp != 2, .(t = min(t[patch == 2]),
			    d = mean(d[patch == 2])),
			by = 'exp'],
       aes(t/120, d / 10)) +
	geom_point(aes(fill = factor(exp)), shape = 21,size = 6, alpha = 0.9, color = 'black') +
	scale_x_continuous('Time difference between patches (min)', breaks = seq(0, 60, 10))+
	scale_y_continuous('Mean distance nest-patch (cm)')+
	scale_fill_viridis_d('Exp') # + coord_cartesian(xlim = c(0,20), ylim = c(0,20))


ggplot(data = melt(data_dist[exp != 2, .(tdiff = min(t[patch == 2]) - min(t[patch == 1]),
					 tp1 = min(t[patch == 1]),
					 d1 = d[which.min(t)], d2 = d[which.min(t[patch == 2])]), by = 'exp'][, c('exp', 'tdiff', 'tp1')],
		   id.vars = 'exp'),
       aes(variable, value / 120)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 5, alpha = 0.5)+
	ylab('Time (min)') + scale_x_discrete('', labels = c('Time difference', 'TP1'))

# ggplot(data = melt(data_dist[exp != 2, .(tdiff = min(t[patch == 2]) - min(t[patch == 1]),
# 					 tp1 = min(t[patch == 1]),
# 					 d1 = d[which.min(t)],
# 					 d2 = d[which.min(t[patch == 2])]),
# 					 by = 'exp'][, .(tp1 = d1/tp1, d2/tdiff), by = 'exp'],
# 		   id.vars = 'exp'),
#        aes(variable, value)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 5, alpha = 0.5)

# ggplot(data = data_dist, aes(t/120, d / 10)) +
#         geom_point(aes(color = factor(exp), shape = factor(exp)), size = 6, alpha = 0.5) +
#         scale_x_continuous('Time (min)', breaks = seq(0, 150, 50))+
#         scale_y_continuous('Distance to nest (cm)')+
#         scale_color_viridis_d('Exp')+
#         scale_shape_manual('Exp', values = c(3,4,8,9,10,12,
#                                              15,17,19))


######### PROBS TIME ##########

# start time: '00:10:00'
# probs <- data.frame(t = c('01:19','02:52','07:55', '08:04', '08:13', '08:35','09:16','09:30','12:03','12:19', '15:11',
#  '12:39', '15:43', '17:39', '19:40', '19:54', '20:40'),
#  p = c(0,0,1,0,1,1,0,1,0,1, 1,
#        1, 1, 1, 1, 1, 1),
#  exp = 1
# )

# start time:
# probs <- data.frame(t = c('15:38', '25:13', '29:13', '32:45', '33:50', '34:57',
#  '03:26', '06:09', '09:17', '09:54', '13:53', '15:28'),
#  p = c(1,1 ,1,1, 1,1, 1,1, 1,1,1, 1),
#  exp = 2)


# start time:
# probs <- data.frame(t = c('01:46', '32:55', '35:13', '36:12', '36:41', '37:12',
#  '39:37','41:16', '42:22', '43:03', '43:34', '45:04'),
#  p = c(1, 1, 1,1,1,1,1,1, 1, 1, 1,1),
#  exp = 3)

# probs <- data.frame(t = c('00:00', '02:29','4:29', '5:34', '05:37', '06:35', '07:10','07:45', '07:58', '08:05', '08:19',
#  '01:41','11:23', '11:48', '14:16', '15:26', '17:46'),
#    p = c(1,1, 0,0, 0,0,1, 1,0, 1, 1,1,1, 1,1,1, 1),
#    exp = 4)

# probs <- data.frame(t = c('02:55', '14:38', '15:01', '17:12', '15:20', '15:30', '16:01',
#  '12:07', '13:49','16:45', '17:07', '17:39', '18:10'),
#  p = c(1,1, 1, 1,0, 1, 1, 1, 1, 1, 1, 1, 1),
#  exp = 5)


# one of the observations is wrong (the probability vector is right, but the timestamps do not match some of the values)
# probs <- data.frame(t = c('13:58', '17:56', '19:54', '20:07', '20:24','25:52','31:20',
#  '01:29', '04:45', '04:49', '05:41', '08:17', '10:16', '10:46', '11:40'),
#  p = c(1, 1, 0, 1, 1,1,1, 1, 1, 1, 0, 1, 1, 0, 1), # ONE IS WRONG
#  exp = 6)

# probs <- data.frame(t = c('14:51', '16:32', '17:49', '18:57', '20:14', '21:28',
#  '01:21', '04:16', '04:45', '05:21', '05:50', '06:25', '07:17', '07:26', '07:39',
#  '08:35', '08:43', '09:07', '10:00'),
#  p = c(1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0,0, 0, 1),
#  exp = 7)

# probs <- data.frame(t = c('03:18', '10:20', '12:45', '15:17', '15:20', '15:30','18:10',
#  '18:33','18:42','20:08','20:45','22:26','22:50', '23:00'),
#  p = c(1,1,1,1,1,0,1,1,1,1,1, 0,1, 1),
#  exp = 8)
#
# probs <- data.frame(t = c('17:14', '23:06', '24:31','24:34', '27:08', '27:59',
#  '04:14', '08:13', '12:57', '14:26', '14:44', '14:54'),
#  p = c(1,1,1,1,1,1, 1, 1, 1, 1, 1, 1),
#  exp = 9)


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