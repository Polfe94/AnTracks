source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

det <- lapply(det, class_ids)
nf <- lapply(nf, class_ids)

det_table <- sapply(det, function(i){
	rws <- i@data[type != 'UNKNOWN', .N, by = 'N_ind']
	setDT(i@data)[rowid(N_ind) == 1 & type != 'UNKNOWN' & N_ind %in% rws[N > 10, N_ind], ][, .N, by = 'type'][, .(N = N / sum(N), type = type)][type == 'Scout', N]

})

nf_table <- sapply(nf[-c(1, 2)], function(i){
	rws <- i@data[type != 'UNKNOWN', .N, by = 'N_ind'][['N']]
	setDT(i@data)[rowid(N_ind) == 1 & type != 'UNKNOWN', ][, .N, by = 'type'][, .(N = N / sum(N), type = type)][type == 'Scout', N]
})


# png()
ggplot(data = data.frame(y = det_table),
       aes(x = 'Proportion of scouts', y = y)) +
	geom_boxplot(outlier.shape = NA, fill = 'mediumpurple',
		     alpha = 0.5, color = 'mediumpurple')+
	geom_jitter(shape = 19, size = 3, color = 'mediumpurple', 
		    width = 0.2)+
	geom_point(shape = 23, size = 6, fill = 'blue', alpha = 0.8,
		   data = data.table(y = det_table)[, .(y = mean(y))])+
	theme(axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      axis.title.x = element_blank(),
	      aspect.ratio = 1.5)+
	scale_y_continuous('Proportion of scouts', limits = c(0, 1.01))


