source('~/research/gits/AnTracks/src/Experiment.R')

load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')

rho_I <- function(exp){
	x <- exp@data[, c('Frame', 'Crossings', 'node')]
	setDT(x)
	
	a <- x[, lapply(.SD, min), .SDcols = 'Frame',by = 'node']
	a2 <- x[, .N, by = c('Frame', 'node')]
	b <- setDT(data.frame(start = c(1, a[['Frame']]),
			      end = c(a[['Frame']], 22000), n = c(0L, as.integer(rownames(a)))))
	b2 <- setDT(data.frame(start = c(1, a2[['Frame']]),
		    end = c(a2[['Frame']], 22000), n = c(0L, a2[['N']])))
	ns <- x[, .N, by = 'Frame']
	
	x <- x[b, on = .(Frame >= start, Frame < end), n := i.n][]
	x <- x[b2, on = .(Frame >= start, Frame < end), d := i.n][]
	x <- x[ns, on = .(Frame == Frame), N := i.N][]
	
	x <- x[, lapply(.SD, mean), .SDcols = c('Crossings', 'n', 'N', 'd'), by = 'Frame']
	x <- x[, c('dIn', 'dIN', 'densI') := .(log(Crossings / n), log((Crossings / N)/n),
					       log(Crossings / d))]
	x[!is.finite(dIn), 'dIn'] <- min(x[is.finite(dIn), 'dIn'])
	x[!is.finite(dIN), 'dIN'] <- min(x[is.finite(dIN), 'dIN'])
	x[!is.finite(densI), 'densI'] <- min(x[is.finite(densI), 'densI'])
	# set(x, j= 'avg_X', value = moving_average(x[['Crossings']],60,30))
	# set(x, j= 'avg_dIn', value = moving_average(-x[['dIn']],60,30))
	# set(x, j= 'avg_dIN', value = moving_average(-x[['dIN']],60,30))
	return(x)
}

test <- rho_I(det[[1]])

plotL_det <- lapply(seq_along(det), function(i){
	data <- rho_I(det[[i]])
	ggplot(data = data,aes(x = Frame)) + 
		geom_path(aes(y = norm_range(moving_average(densI, 2, 1), 0,1 )),
			  color = 'darkorange')+
		#geom_path(aes(y = norm_range(moving_average(Crossings, 2, 1), 0, 1))) + 
		geom_path(aes(y = norm_range(moving_average(dIn, 2, 1), 0, 1)), color = 'darkred')+

		geom_path(aes(y = norm_range(moving_average(N, 2, 1), 0,1 )), color = 'purple') +
		# geom_path(aes(y = norm_range(cumsum(norm_range(dIn, 0, 1)), 0, 1)),
		# 	  color = 'darkorange', linewidth = 1)+
		scale_x_continuous('Time (min)', breaks = seq(0, 21600, 15*120),
				   labels = seq(0, 180, 15))+
		
		scale_y_continuous('')
})

ggarrange(plotlist = plotL_det)


plotL_sto <- lapply(seq_along(sto), function(i){
	data <- rho_I(sto[[i]])
	ggplot(data = data,aes(x = Frame)) + 
		geom_path(aes(y = norm_range(moving_average(densI, 2, 1), 0,1 )),
			  color = 'darkorange')+
		#geom_path(aes(y = norm_range(moving_average(Crossings, 2, 1), 0, 1))) + 
		geom_path(aes(y = norm_range(moving_average(dIn, 2, 1), 0, 1)), color = 'darkred')+
		
		geom_path(aes(y = norm_range(moving_average(N, 2, 1), 0,1 )), color = 'purple') +
		# geom_path(aes(y = norm_range(cumsum(norm_range(dIn, 0, 1)), 0, 1)),
		# 	  color = 'darkorange', linewidth = 1)+
		scale_x_continuous('Time (min)', breaks = seq(0, 21600, 15*120),
				   labels = seq(0, 180, 15))+
		
		scale_y_continuous('')
})

ggarrange(plotlist = plotL_sto)








# mlt <- reshape2::melt(det2, id.vars = c('Frame', 'n', 'N'))
# 
# ggplot(data = mlt[mlt$variable %in% c('avg_dIn', 'avg_dIN, avg_X'), c('variable', 'value', 'Frame')],
#        aes(Frame, value, color = variable)) + geom_path()
# 
# ggplot(data = det2, aes(Frame, norm_range(-dIn, 0, 1))) + geom_path()
# ggplot(data = det2, aes(Frame, norm_range(moving_average(Crossings, 60, 30), 0, 1))) + geom_path()+
# 	geom_path(aes(y = norm_range(moving_average(N, 60, 30), 0, 1)), color = 'red')
# ggplot(data = det2, aes(Frame, norm_range(moving_average(-dIn, 60, 30), 0, 1))) + geom_path()+
# 	geom_path(aes(y = norm_range(moving_average(N, 60, 30), 0, 1)), color = 'red')
# ggplot(data = det2, aes(Frame, norm_range(moving_average(-dIN, 60, 30), 0, 1))) + geom_path()+
# 	geom_path(aes(y = norm_range(moving_average(N, 60, 30), 0, 1)), color = 'red')
# 
# ggplot(data = det2, aes(Frame, 
# 			norm_range(log((Crossings/N)/n), 0, 1))) +
# 	geom_line() 
# ggplot(data = det2, aes(Frame, log(Crossings/n))) + geom_line()
# 
# 
# 
# ggplot(data = l, aes(Frame, log((Crossings/N)/n))) + geom_line() 
# ggplot(data = test, aes(Frame, log((Crossings/N)/n))) + geom_line() 
# ggplot(data = test, aes(Frame, log(Crossings/n))) + geom_line()

