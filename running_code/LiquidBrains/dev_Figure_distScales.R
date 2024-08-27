source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')
source('~/research/gits/AnTracks/src/fit_functions.R')
# load('~/research/gits/AnTracks/results/mov_rho_R.RData')
library(segmented)
library(trajr)

bin_function <- function(v, res = 10, minv = 0){
	cut(v, breaks = seq(minv, ceiling(max(v)/res)*res, res))
}

DET_scouts <- rbindlist(lapply(det, function(p){
	t <- min(rbindlist(p@food)[['t']])
	data <- setDT(p@data)[Frame <= t]
	# data <- setDT(p@data)
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	data[, .(maxd = max(d), meand = mean(d)), by = 'N_ind']
}))
DET_scouts[['tag']] <- 'DET'

NFD_scouts <- rbindlist(lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	data[, .(maxd = max(d), meand = mean(d)), by = 'N_ind']
}))
NFD_scouts[['tag']] <- 'NFD'


data_binned <- rbind(DET_scouts, NFD_scouts)
data_binned[, bin := bin_function(meand, res = 10, minv = 80), by = c('tag')]
data <- data_binned[, .N, by = c('bin', 'tag')]
fit_data <- data[, .(x = ((as.numeric(bin)+8)*10) / 10,
		 y = N), by = 'tag']
lmodels <- lapply(unique(fit_data[['tag']]),
		  function(m){
		  	data <- fit_data[tag == m]
		  	mdl <- lm(formula = log(y) ~ x, data = data)
		  	mdl$call <- as.call(list(
		  		quote(lm),
		  		formula = quote(log(y) ~ x),
		  		data = substitute(fit_data[tag == m], list(m = m))
		  	))
		  	mdl
		  })
smodels <- lapply(lmodels, segmented)
breakpoints <- unlist(sapply(smodels, function(i) i$psi[, 'Est.']))
st.err <- round(unlist(sapply(smodels, function(i) i$psi[, 'St.Err'])), 1)
x_values <- lapply(smodels, function(i){
	seq(min(i$model[, 2]), max(i$model[, 2]), length.out = 1000)
})

xv <- lapply(seq_along(x_values), function(i){
	x_values[[i]][c(1, which(x_values[[i]] > breakpoints[i])[1]-1,
		      which(x_values[[i]] > breakpoints[i])[1], length(x_values[[i]]))]
})
pred_ <- lapply(seq_along(xv), function(i) predict(smodels[[i]],
						   newdata = data.frame(x = xv[[i]])))


ggplot(data = fit_data, aes(x, y)) + 

	scale_x_continuous('Average distance to nest (cm)', breaks = seq(0, 140, 10))+
	# scale_y_continuous() # seems single exponential
	scale_y_log10('Frequency')+ # reveals a double exponential+
	facet_wrap(~tag, scales = 'free')+
	geom_line(data = data.frame(x = unlist(xv),
				    y = exp(unlist(pred_)),
				    tag = c(rep('DET', length(xv[[1]])),
				    	rep('NFD', length(xv[[1]])))),
		  aes(x, y))+
	annotation_logticks(sides = 'l')+
	geom_vline(data = data.frame(x = breakpoints,
				     tag = c('DET', 'NFD')),
		   linetype = 2, aes(xintercept = x))+
geom_point(size = 4.5, shape = 1) 









foods <- rbindlist(lapply(det, function(i){
	d <- rbindlist(i@food)
	idx <- which.min(d[['t']])
	d[idx]
}))[, node := get_node(x, y), by = seq_along(det)]


ids <- lapply(seq_along(det), function(i){
	id <- setDT(det[[i]]@data)[Frame == foods[i, t] & node == foods[i, node]][['N_ind']]
	data <- det[[i]]@data[N_ind == id]
	Dt <- diff(data[, Frame])
	stops <- which(Dt != 1)
	stop_duration <- Dt[Dt != 1]
	list(data = data, stops = stops, stop_duration = stop_duration, interactions = sum(data[['Crossings']]))
})

result_LR <- rbindlist(lapply(ids, function(i){
	sbst <- i[['data']]
	dmatrix <- pdist(as.matrix(sbst[, c('Xmm', 'Ymm')]),
			 as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(sbst, j = 'd', value = as.numeric(dmatrix))
	xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
	data.frame(s = TrajStraightness(xy[xy$Frame <= sbst[which.max(d), Frame], ]),
		   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
		   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
		   v2 = Mod(TrajMeanVelocity(xy)),
		   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
		   l = TrajLength(xy),
		   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
		   t = TrajDuration(xy),
		   sin = TrajSinuosity2(xy))
	
}), idcol = TRUE)

result_SR <- rbindlist(lapply(seq_along(det), function(i){
	id <- setDT(det[[i]]@data)[Frame == foods[i, t] & node == foods[i, node]][['N_ind']]
	data <- det[[i]]@data[N_ind != id]
	rbindlist(lapply(unique(data[['N_ind']]), function(j){
		sbst <- data[N_ind == j]
		if(nrow(sbst) > 10){
			dmatrix <- pdist(as.matrix(sbst[, c('Xmm', 'Ymm')]),
					 as.matrix(hex[hex$node == 634, c('x', 'y')]))
			set(sbst, j = 'd', value = as.numeric(dmatrix))
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy))
		}

	}), idcol = 'id')


	
}), idcol = TRUE)

result_SR <- rbindlist(lapply(seq_along(det), function(i){
	id <- setDT(det[[i]]@data)[Frame == foods[i, t] & node == foods[i, node]][['N_ind']]
	data <- det[[i]]@data[N_ind != id]

}))

boxplot(result_LR[['v']])
boxplot(result_LR[['acc']])
boxplot(result_LR[['d']]/10)
boxplot(result_LR[['t']]/120)

boxplot(result_LR[['s']])


maxt <- 3000
maxd <- 500
DET_metrics <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('det_LR_',i), scout = 'LR', condition = 'det')
		}
	}))
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('det_SR_',i), scout = 'SR', condition = 'det')
		}
		
		
	}))
	rbind(result_LR, result_SR)
}), idcol = TRUE)

NFD_metrics <- rbindlist(lapply(nf, function(p){
	data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('nf_LR_',i), scout = 'LR', condition = 'nf')
		}
	}))
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('nf_SR_',i), scout = 'SR', condition = 'nf')
		}
		
		
	}))
	rbind(result_LR, result_SR)
}), idcol = TRUE)

mlt_det_metrics <- data.table(melt(DET_metrics, id.vars = c('tag', 'scout', 'condition', '.id')))
mlt_nf_metrics <- data.table(melt(NFD_metrics, id.vars = c('tag', 'scout', 'condition', '.id')))

## selected cols 
# sdcols <- c('v', 'acc', 'l', 'd', 't')
sdcols <- c('v', 'acc', 'd', 't')
general_metrics <- rbind(mlt_det_metrics[variable %in% sdcols], mlt_nf_metrics[variable %in% sdcols])
general_metrics[variable %in% sdcols[-4], value := value / 10]
general_metrics[variable %in% sdcols[4], value := value / 60]
# variable = c(v = TeX("Velocity ($cm*s^{-1}$)"),
# 	     acc = TeX("Acceleration ($cm*s^{-2}$)"),
# 	     l ="Cumulative distance (cm)",
# 	     d ='Max distance (cm)',
# 	     t = 'Time in arena (min)'), .default = label_parsed)

labeller <- function(string){
	TeX(paste0(string))
}

## CODE TO OBTAIN THE PAIRWISED WILCOXON TESTS
# tests <- general_metrics
# tests[['label']] <- apply(tests[, c('scout', 'condition')], 1, paste0, collapse = '', sep='')
# 
# for(v in sdcols){
# 	print(v)
# 	print(
# 	pairwise.wilcox.test(x = tests[variable == v, value], 
# 			     g= tests[variable == v, label], p.adjust.method = 'bonferroni'))
# }


pvalues <- data.table(scout = rep(c('LR', 'SR'), 8), 
		      condition = rep(c(rep('det', 2), rep('nf', 2)), 4),
		      label = c('a', 'a', 'b', 'b',
		      	  'a', 'bc', 'ab', 'c',
		      	  'a', 'a', 'b', 'c',
		      	  'a', 'a', 'b', 'b'),
		      variable = c(vapply(sdcols, rep,FUN.VALUE = character(4), 4), recursive = T),
		      value = NA)
pvalues[['y']] <- sapply(1:nrow(pvalues), function(i) general_metrics[scout == pvalues[i, scout] & 
								      	variable == pvalues[i, variable]&
								      	condition == pvalues[i, condition],
								      min(quantile(value, p = 0.75) + 1.5 * (quantile(value, p = 0.75)-quantile(value, p = 0.25)),
								          max(value))  + median(value)* 0.1])
pvalues[['y']] <- pvalues[['y']] + c(rep(0.25, 4), 0.5, 1, 1,1, 7.5,20, 10, 15, 0.5,1.5, 1.25, 1.25)
# pvalues[['y']] <- pvalues[['y']] + c(rep(0, 5), 0.5, 0.5,0.5, 0,10, 0, 10, 0,0.5, 0.5, 0.5)

# pvalues[variable == 't', 'y'] <- pvalues[variable == 't', y / 60]
# pvalues[variable == 'd', 'y'] <- pvalues[variable == 'd', y / 100]
pvalues[['variable']] <- factor(pvalues[['variable']], levels = sdcols)

panel_A <- ggplot(data = general_metrics, aes(scout, value, fill = condition)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = NA)+
	geom_jitter(size = 3, alpha = 0.1, 
		    position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.15),
		    show.legend = FALSE)+
	## for some reason labels are sorted wrong...
	geom_text(data = pvalues, aes(label = label, y = pvalues[order(variable, scout), y]),
		  size = 7, position = position_nudge(x = c(-0.19, -1+0.19, 1-0.19, 0.19)))+
	scale_fill_manual('', labels = c('DET', 'NFD'), values = c('mediumpurple','gold3'))+
	ylab('')+xlab('') +
	facet_wrap(~ factor(variable, 
			    labels = c("Velocity ($cm\\cdot s^{-1}$)",
			    	   "Acceleration ($cm\\cdot s^{-2}$)",
			    	   # "`Cumulative distance (cm)`",
			    	   '`Maximum distance (cm)`',
			    	   '`Time in arena (min)`')),
		   scales = 'free_y', 
		   labeller = as_labeller(labeller, default = label_parsed)) +
	theme(legend.position = c(0.35, 0.35), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))+
	geom_hline(data = melt(data.table(result_LR[, .(v = v / 10,
						   acc = acc /10,
						   d = d /10,
						   t = t / 120
	)])), aes(yintercept = value))+
	geom_hline(data = melt(data.table(setDT(result_SR)[, .(v = v / 10,
							acc = acc /10,
							d = d /10,
							t = t / 120
	)])), aes(yintercept = value), color = 'red')



# 
# 
# 
# 
# maxt <- 20*120
# 
# library(trajr)
# library(arrow)
# library(egg)
# library(latex2exp)
# library(segmented)
# 
# DET_scouts <- rbindlist(lapply(c(5, 10, 15, 20, 25) * 120,
# 			       function(t){
# 			       	
# 			       	rbindlist(lapply(det, function(p){
# 			       		data <- setDT(p@data)[Frame < t + min(Frame)]
# 			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
# 			       		set(data, j = 'd', value = as.numeric(dmatrix))
# 			       		data[, .(maxd = max(d),
# 			       			 maxt = t,
# 			       			 meand = mean(d)), by = 'N_ind']
# 			       	}))
# 			       }))
# DET_scouts[['tag']] <- 'DET'
# 
# NFD_scouts <- rbindlist(lapply(c(5, 10, 15, 20, 25) * 120,
# 			       function(t){
# 			       	
# 			       	rbindlist(lapply(nf[-c(1, 2)], function(p){
# 			       		data <- setDT(p@data)[Frame < maxt + min(Frame)]
# 			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
# 			       		set(data, j = 'd', value = as.numeric(dmatrix))
# 			       		data[, .(maxd = max(d),
# 			       			 maxt = t,
# 			       			 meand = mean(d)), by = 'N_ind']
# 			       	}))
# 			       }))
# NFD_scouts[['tag']] <- 'NFD'
# 
# # NFD_scouts <- rbindlist(lapply(c(180) * 120,
# # 			       function(t){
# # 			       	
# # 			       	rbindlist(lapply(nf[-c(1, 2)], function(p){
# # 			       		data <- setDT(p@data)[Frame < maxt + min(Frame)]
# # 			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
# # 			       		set(data, j = 'd', value = as.numeric(dmatrix))
# # 			       		data[, .(maxd = max(d),
# # 			       			 maxt = t,
# # 			       			 meand = mean(d)), by = 'N_ind']
# # 			       	}))
# # 			       }))
# # NFD_scouts[['tag']] <- 'NFD'
# 
# bin_resolution <- 10
# data_binned <- rbind(DET_scouts, NFD_scouts)
# ## MAX DIST
# data_binned[['bin']] <- cut(data_binned[['maxd']],
# 			    breaks = seq(0, ceiling(max(data_binned[['maxd']])/ 100)*100,
# 			    	     bin_resolution))
# data <- data_binned[, .N, by = c('bin', 'maxt', 'tag')][, logN := log(N)]
# 
# 
# ggplot(data = data[maxt == 3000 & N > 1],
#        aes(x = as.numeric(bin)*bin_resolution, y = N)) + 
# 	geom_point() +
# 	scale_x_continuous()+
# 	scale_y_log10()
# 
# ## MEAN DIST
# data_binned[['bin']] <- cut(data_binned[['meand']],
# 			    breaks = seq(0, ceiling(max(data_binned[['meand']])/ 100)*100,
# 			    	     bin_resolution))
# data <- data_binned[, .N, by = c('bin', 'maxt', 'tag')][, logN := log(N)]
# 
# 
# ggplot(data = data[maxt == 3000 & N > 1 & 
# 		   	as.numeric(bin)*bin_resolution > 100 &
# 		   	tag == 'DET'],
#        aes(x = as.numeric(bin)*bin_resolution, y = N)) + 
# 	geom_point() +
# 	scale_x_continuous('<Distance to nest> (cm)', breaks = seq(100, 1400, 100),
# 			   labels = seq(100, 1400, 100)/10)+
# 	# scale_y_continuous() # seems single exponential
# 	scale_y_log10('Freq') # reveals a double exponential
# 
# 
# fit_data <- data[maxt == 3000 & N > 1 & 
# 		 	as.numeric(bin)*bin_resolution > 100 &
# 		 	tag == 'DET', .(y = N,
# 		 			x = as.numeric(bin)*bin_resolution)]
# lmodel <- lm(formula = log(y) ~ x, data = fit_data)
# smodel <- segmented(lmodel)
# breakpoint <- smodel$psi[, 'Est.']
# st.err <- smodel$psi[, 'St.Err']
# x_values <- seq(min(smodel$model[, 2]), max(smodel$model[, 2]), length.out = 1000)
# xv <- x_values[c(1, which(x_values > breakpoint)[1]-1,
# 	   which(x_values > breakpoint)[1], length(x_values))]
# pred_ <- predict(smodel, newdata = data.frame(x = xv))
# 
# 
# ggplot(data = data[maxt == 3000 & N > 1 & 
# 		   	as.numeric(bin)*bin_resolution > 100 &
# 		   	tag == 'DET'],
#        aes(x = as.numeric(bin)*bin_resolution, y = N)) + 
# 	geom_point() +
# 	scale_x_continuous('<Distance to nest> (cm)', breaks = seq(0, 1400, 100),
# 			   labels = seq(0, 1400, 100)/10)+
# 	# scale_y_continuous() # seems single exponential
# 	scale_y_log10('Freq')+ # reveals a double exponential
# 	geom_line(data = data.frame(x = xv,
# 				    y = exp(pred_)),
# 		  aes(x, y))+
# 	annotation_logticks(sides = 'l')+
# 	coord_cartesian(ylim = c(1, 300))+
# 	geom_vline(linetype = 2, xintercept = breakpoint)
# 
# 
# 
# fit_data <- data[as.numeric(bin)*bin_resolution > 100 &
# 		 	tag == 'NFD', .(y = N,
# 		 			x = as.numeric(bin)*bin_resolution)][x < 900]
# lmodel <- lm(formula = log(y) ~ x, data = fit_data)
# smodel <- segmented(lmodel)
# breakpoint <- smodel$psi[, 'Est.']
# st.err <- smodel$psi[, 'St.Err']
# x_values <- seq(min(smodel$model[, 2]), max(smodel$model[, 2]), length.out = 1000)
# xv <- x_values[c(1, which(x_values > breakpoint)[1]-1,
# 		 which(x_values > breakpoint)[1], length(x_values))]
# pred_ <- predict(smodel, newdata = data.frame(x = xv))
# ggplot(data = data[as.numeric(bin)*bin_resolution > 100 &
# 		   	as.numeric(bin)*bin_resolution < 900 &
# 		   	tag == 'NFD'],
#        aes(x = as.numeric(bin)*bin_resolution, y = N)) + 
# 	geom_point() +
# 	scale_x_continuous('<Distance to nest> (cm)', breaks = seq(0, 1400, 100),
# 			   labels = seq(0, 1400, 100)/10)+
# 	# scale_y_continuous() # seems single exponential
# 	scale_y_log10('Freq')+ # reveals a double exponential
# 	geom_line(data = data.frame(x = xv,
# 				    y = exp(pred_)),
# 		  aes(x, y))+
# 	annotation_logticks(sides = 'l')+
# 	coord_cartesian(ylim = c(1, 300))+
# 	geom_vline(linetype = 2, xintercept = breakpoint)
# fit_ <- fit_doub_nexp(data[maxt == 3000 & N > 1 & 
# 			   	as.numeric(bin)*bin_resolution > 100 &
# 			   	tag == 'DET', .(y = N,
# 			   			x = as.numeric(bin)*bin_resolution)],
# 			   	a_0 = 0.5, k_0 = 10**-3, a2_0 = 0.1, k2_0 = 10**-4,
# 			   	alims = c(-10,10000), klims = c(-10, 10),
# 			   	a2lims =  c(-10,10000), k2lims = c(-10, 10))
# pred_ <- data.table(x = data[maxt == 3000 & N > 1 & 
# 			 	as.numeric(bin)*bin_resolution > 100 &
# 			 	tag == 'DET', .(y = N,
# 			 			x = as.numeric(bin)*bin_resolution)][['x']],
# 			 	y = doub_nexp(data[maxt == 3000 & N > 1 & 
# 			 			   	as.numeric(bin)*bin_resolution > 100 &
# 			 			   	tag == 'DET', .(y = N,
# 			 			   			x = as.numeric(bin)*bin_resolution)][['x']],
# 			 			   	k = fit_[['k']], a = fit_[['a']],
# 			 			   	k2 = fit_[['k2']], a2 = fit_[['a2']]))
# 
# 
# 
# ggplot(data = smodel$model, aes(x, exp(`log(y)`))) + geom_point()+scale_y_log10()
# 
# fits <- lapply(seq_along(data_list), function(i){
# 	fit_doub_nexp(data_list[[i]], a_0 = 0.5, k_0 = 10**-3, a2_0 = 0.1, k2_0 = 10**-4,
# 		      alims = c(-10, 10), klims = c(-10, 10), a2lims = c(-10, 10), k2lims = c(-10, 10))
# })
# names(fits) <- names(data_list)
# pred_data <- lapply(seq_along(data_list), function(i){
# 	data.table(x = data_list[[i]][['x']], 
# 		   y = doub_nexp(data_list[[i]][['x']], a = fits[[i]][['a']], 
# 		   	      k = fits[[i]][['k']], a2 = fits[[i]][['a2']], k2 = fits[[i]][['k2']]))
# })
# names(pred_data) <- names(data_list)
# 
# data_list <- rbindlist(data_list, idcol = TRUE)
# labels_unique <- do.call('rbind', strsplit(data_list[['.id']], '_'))
# data_list[['scout']] <- labels_unique[, 1]
# data_list[['condition']] <- labels_unique[, 2]
# 
# pred_data <- rbindlist(pred_data)
# pred_data[['scout']] <- labels_unique[, 1]
# pred_data[['condition']] <- labels_unique[, 2]
# 
# 
# png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_S1_2.png', 6000, 3000, res = 350)
# ggplot(data = rbind(DET_scouts, NFD_scouts), aes(maxd, fill = tag)) + 
# 	geom_density(aes(fill = tag, y = after_stat(density)),
# 		     alpha = 0.9,  color = 'black')+
# 	geom_histogram(alpha = 0.5, color = 'black', aes(y = after_stat(density)),
# 		       bins = 25) +
# 	facet_wrap(~ maxt, scales = 'free',
# 		   labeller = as_labeller(function(i){
# 		   	paste0(
# 		   		as.numeric(i)/120,
# 		   		' min')
# 		   }))+
# 	ylab('Density') + 
# 	scale_x_continuous('Radius (mm)', breaks = seq(0, 1500, 250))+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
# 	theme(legend.position = c(0.8, 0.25))+
# 	geom_vline(xintercept = 650, linetype = 2)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ############# COSES DE LA VIDA SENSE INTERACCIONS ################
# 
# 
# DET_scouts <- rbindlist(lapply(c(5, 10, 15, 20, 25) * 120,
# 			       function(t){
# 			       	
# 			       	rbindlist(lapply(det, function(p){
# 			       		data <- setDT(p@data)[Frame < t + min(Frame)]
# 			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
# 			       		set(data, j = 'd', value = as.numeric(dmatrix))
# 			       		ints <- data[, .(i = sum(Crossings)), by = 'N_ind']
# 			       		inds <- ints[i == 0, 'N_ind']
# 			       		data[N_ind %in% inds, .(maxd = max(d),
# 			       			 maxt = t,
# 			       			 meand = mean(d)), by = 'N_ind']
# 			       	}))
# 			       }))
# DET_scouts[['tag']] <- 'DET'
# 
# NFD_scouts <- rbindlist(lapply(c(5, 10, 15, 20, 25) * 120,
# 			       function(t){
# 			       	
# 			       	rbindlist(lapply(nf[-c(1, 2)], function(p){
# 			       		data <- setDT(p@data)[Frame < maxt + min(Frame)]
# 			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
# 			       		set(data, j = 'd', value = as.numeric(dmatrix))
# 			       		data[, .(maxd = max(d),
# 			       			 maxt = t,
# 			       			 meand = mean(d)), by = 'N_ind']
# 			       	}))
# 			       }))
# NFD_scouts[['tag']] <- 'NFD'
# 
# # NFD_scouts <- rbindlist(lapply(c(180) * 120,
# # 			       function(t){
# # 			       	
# # 			       	rbindlist(lapply(nf[-c(1, 2)], function(p){
# # 			       		data <- setDT(p@data)[Frame < maxt + min(Frame)]
# # 			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
# # 			       		set(data, j = 'd', value = as.numeric(dmatrix))
# # 			       		data[, .(maxd = max(d),
# # 			       			 maxt = t,
# # 			       			 meand = mean(d)), by = 'N_ind']
# # 			       	}))
# # 			       }))
# # NFD_scouts[['tag']] <- 'NFD'
# 
# 
# 
# 
# 
# 
# 
# 
# test <- det[[1]]@data[Frame < 25*120 +  min(Frame)]
# 
# dmatrix <- pdist(as.matrix(test[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
# set(test, j = 'd', value = as.numeric(dmatrix))
# ints <- test[, .(i = sum(Crossings)), by = 'N_ind']
# inds <- as.numeric(as.character(ints[i == 0, N_ind]))
# test[N_ind %in% inds, .(maxd = max(d),
# 			meand = mean(d)), by = 'N_ind']
# 
# hist(test[N_ind %in% inds, .(maxd = max(d),
# 			     meand = mean(d)), by = 'N_ind'][['meand']])
