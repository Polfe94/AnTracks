source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/fit_functions.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

det <- lapply(det, class_ids)
nf <- lapply(nf, class_ids)

det_stats <- rbindlist(lapply(det, move_stats))
nf_stats <- rbindlist(lapply(nf[-c(1, 2)], move_stats))

stats <- rbindlist(list(det = det_stats, nf = nf_stats), idcol = 'exp')
stats <- na.omit(stats)

## 
stats <- stats[d > 1]
stats <- stats[!(type == 'Recruit' & exp == 'det' & d > 500)]

stats_mlt <- data.table(melt(stats[, c('exp', 'v', 'acc', 'd', 't', 'type')],
		  id.vars = c('exp', 'type')))
stats_mlt[['interaction']] <- apply(stats_mlt[, c('exp', 'type')], 1,
				    paste0, collapse = '', sep='')

## transformation of units 
stats_mlt[variable == 'v', value := value / 10] # to cm/s
stats_mlt[variable == 'acc', value := value / 10] # to cm/s^2
stats_mlt[variable == 'd', value := value / 10] # to cm
stats_mlt[variable == 't', value := value / 60] # to min




pairwise.wilcox.test(stats_mlt[variable == 'v', value], 
		     stats_mlt[variable == 'v', interaction], 
		     p.adjust.method = 'bonferroni')

pairwise.wilcox.test(stats_mlt[variable == 'acc', value], 
		     stats_mlt[variable == 'acc', interaction], 
		     p.adjust.method = 'bonferroni')

pairwise.wilcox.test(stats_mlt[variable == 'd', value], 
		     stats_mlt[variable == 'd', interaction], 
		     p.adjust.method = 'bonferroni')

pairwise.wilcox.test(stats_mlt[variable == 't', value], 
		     stats_mlt[variable == 't', interaction], 
		     p.adjust.method = 'bonferroni')


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


plot_ <- ggplot(data = stats_mlt, aes(type, value, fill = exp)) + 
	# geom_violin(alpha = 0.6)+
	geom_boxplot(alpha = 0.6, outlier.shape = NA)+
	# geom_jitter(size = 3, alpha = 0.1, 
	# 	    position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.225),
	# 	    show.legend = FALSE)+
	## for some reason labels are sorted wrong...
	# geom_text(data = pvalues, aes(label = label, y = pvalues[order(variable, scout), y]),
	# 	  size = 7, position = position_nudge(x = c(-0.19, -1+0.19, 1-0.19, 0.19)))+
	
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	facet_wrap(~ factor(variable, 
			    labels = c("Velocity ($cm\\cdot s^{-1}$)",
			    	   "Acceleration ($cm\\cdot s^{-2}$)",
			    	   # "`Cumulative distance (cm)`",
			    	   '`Maximum distance (cm)`',
			    	   '`Time in arena (min)`')),
		   scales = 'free_y', 
		   labeller = as_labeller(labeller_TeX, default = label_parsed)) +
	theme(legend.position = c(0.35, 0.35), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))+
	scale_y_continuous(limits = function(x) ifelse (unique(stats_mlt$variable) == 't',
						   c(0, 20), range(x)))
















DET_scouts <- rbindlist(lapply(det, function(p){
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
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(maxL >= 100){
			idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
		} else {
			idxs <- nrow(xy)
		}
		straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
		data.frame(l = Ls[idxs], s = straightness, 
			   tag = paste0('det_LR_',i), scout = 'LR', condition = 'det')
		
	}), idcol = TRUE)
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 2){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness, 
				   tag = paste0('det_SR_',i), scout = 'SR', condition = 'det')
		}
		
	}), idcol = TRUE)
	rbind(result_LR, result_SR)
}), idcol = TRUE)





setMethod('get_straightness', 'Experiment', function(.Object, ...){
	require(trajr)
	.Object <- class_ids(.Object, ...)
	data <- .Object@data[type != 'UNKNOWN']
	ids <- data[, unique(N_ind)]
	
	rbindlist(lapply(seq_along(ids), function(i){
		sbst <- data[N_ind == ids[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		
		# slight filter
		if(nrow(xy) > 10){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		} else {
			df <- data.frame(l = NA, s = NA)
		}
		df[['type']] <- sbst[, unique(type)]
		df[d > 0]
	}), idcol = 'ind')
})



