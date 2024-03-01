tseq <- seq(5*120, 25*120, 5*120) # up to the maximum exploration time (25 mins)
LR_det <- vector('list', length(tseq))
names(LR_det) <- tseq
for(i in seq_along(tseq)){
	print(paste0('Iter: ', i))
	LR_det[[i]] <- lapply(det, move_stats, scouts = 'LR', maxt = tseq[i], maxd = 500, min_pos = 10)
}
LR_det <- rbindlist(lapply(LR_det, rbindlist), idcol = TRUE)[, scouts := 'LR']
mlted_LR_det <- melt(LR_det, id.vars = c('.id', 'scouts'))

SR_det <- vector('list', length(tseq))
names(SR_det) <- tseq
for(i in seq_along(tseq)){
	print(paste0('Iter: ', i))
	SR_det[[i]] <- lapply(det, move_stats, scouts = 'SR', maxt = tseq[i], maxd = 500, min_pos = 10)
}
SR_det <- rbindlist(lapply(SR_det, rbindlist), idcol = TRUE)[, scouts := 'SR']
mlted_SR_det <- melt(SR_det, id.vars = c('.id', 'scouts'))

ggplot(data = rbind(mlted_LR_det, mlted_SR_det),
       aes(scouts, value, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	# geom_jitter(aes(color = .id), 
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio') +
	facet_wrap(~ variable, scales = 'free')
