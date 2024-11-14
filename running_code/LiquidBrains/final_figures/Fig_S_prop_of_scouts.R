source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

det <- lapply(det, class_ids) # old method
nf <- lapply(nf, class_ids)

det_table <- sapply(det, function(i){
	rws <- i@data[type != 'UNKNOWN', .N, by = 'N_ind']
	setDT(i@data)[rowid(N_ind) == 1 & type != 'UNKNOWN' & N_ind %in% rws[N > 10, N_ind], ][, .N, by = 'type'][, .(N = N / sum(N), type = type)][type == 'Scout', N]

})

# nf_table <- sapply(nf[-c(1, 2)], function(i){
# 	rws <- i@data[type != 'UNKNOWN', .N, by = 'N_ind'][['N']]
# 	setDT(i@data)[rowid(N_ind) == 1 & type != 'UNKNOWN', ][, .N, by = 'type'][, .(N = N / sum(N), type = type)][type == 'Scout', N]
# })


ggplot(data = data.frame(y = det_table),
       aes(x = 'Proportion of scouts', y = y)) +
	geom_boxplot(outlier.shape = NA,
		     alpha = 0.5, size = 1, width = 0.5)+
	geom_dotplot(binaxis = 'y', stackdir = 'center', fill = 'grey80',
		     binwidth = 1/100, dotsize = 3) +
	geom_hline(yintercept = round(mean(det_table), 1), linetype = 2, linewidth = 1.1)+
	theme(axis.text.y = element_blank(),
	      axis.ticks.y = element_blank(),
	      axis.title.y = element_blank(),
	      axis.line.y = element_blank())+
	scale_y_continuous('Proportion of scouts', limits = c(0, 1))+
	coord_flip()


