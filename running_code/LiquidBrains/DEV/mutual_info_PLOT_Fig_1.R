library(latex2exp)
source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/det.RData')

load("/home/polfer/research/gits/AnTracks/results/MI_det.RData")
load("/home/polfer/research/gits/AnTracks/results/MI_default.RData")
mi_sim <- rbindlist(lapply(seq_along(mutual_info_result), function(i){
	x <- mutual_info_result[[i]]
	n <- names(x)
	rbindlist(lapply(seq_along(x), function(ii){
		data.frame(t = as.numeric(n)[ii], x = as.numeric(x[[ii]]))
	}))
}), idcol = TRUE)

mi_exp <- rbindlist(lapply(seq_along(det_mi), function(i){

	x <- det_mi[[i]]
	n <- names(x)
	rbindlist(lapply(seq_along(x), function(ii){
		data.frame(t = as.numeric(n)[ii], x = as.numeric(x[[ii]]))
	}))
	

}), idcol = TRUE)


ggplot(data = mi_sim, aes(t, x, group = t)) + geom_boxplot()
ggplot(data = mi_sim[, .(x = mean(x)), by = 't'], aes(t, x))+
	geom_point(size = 3, shape = 21) +
	geom_path()+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120,
			   labels = seq(0, 150, 50))+
	scale_y_continuous('<MI>')
ggplot(data = mi_exp[, .(x = mean(x)), by = 't'], aes(t, x))+
	geom_point(size = 3, shape = 21) +
	geom_path()+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120,
			   labels = seq(0, 150, 50))+
	scale_y_continuous('<MI>')




# png('~/research/gits/AnTracks/plots/provisional_MI.png', 4000)
d <- ggplot(data = mi_sim[, .(x = mean(x)), by = 't'][!is.na(x), .(x = norm_range(x, 0, 1), t = t)], aes(t, x), alpha = 0.5)+
	geom_point(size = 3, shape = 21, aes(color = 'Simulations')) +
	geom_path(aes(color = 'Simulations'))+
	geom_point(size = 3, shape = 21, 
		   data = mi_exp[, .(x = mean(x)), by = 't'][!is.na(x), .(x = norm_range(x, 0, 1), t = t)], aes(color = 'Experiments')) +
	geom_path(data = mi_exp[, .(x = mean(x)), by = 't'][!is.na(x), .(x = norm_range(x, 0, 1), t = t)], aes(color = 'Experiments'))+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120,
			   labels = seq(0, 150, 50))+
	scale_y_continuous(TeX('$<MI_{sim}>$'), labels = seq(0, 1, length.out = 5)*
			   	round(max(mi_sim[, .(x = mean(x)), by = 't'][['x']], na.rm = TRUE), 4),
			   sec.axis = sec_axis(trans = function(i) i, 
			   		    labels = seq(0, 1, length.out = 5)*
			   		    	round(max(mi_exp[, .(x = mean(x)), by = 't'][['x']], na.rm = TRUE), 4),
			   		    name = TeX('$<MI_{exp}>$'))) +
	scale_color_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.position = c(0.825, 0.775), 
	      legend.background = element_rect(color = 'black', fill = NA),
	      legend.title = element_blank())



png('/home/polfer/research/gits/AnTracks/plots/summary_default_prov.png', 8000, 6000, res = 550)
grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5, plot.title = element_text(face = 'bold'))+
			  	ggtitle('A'), b+ggtitle('B')+theme(plot.title = element_text(face = 'bold')), 
			  ggarrange(c2, c, ncol = 1, common.legend = TRUE,
			  	  legend = 'bottom', labels = c('C'), hjust = -8.25, vjust = 0.5 ), 
			  d + ggtitle('D')+theme(plot.title = element_text(face = 'bold'))),
			  layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
			  		      c(3, 2), c(3, 4), c(3, 4)))
dev.off()
