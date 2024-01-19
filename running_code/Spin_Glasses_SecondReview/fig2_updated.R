source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')
load('~/research/gits/AnTracks/data/nf.RData')

library(scales)

matrix(data = c(det_food_times, sto_food_times)/120, byrow = TRUE, ncol = 2, 
       dimnames = list(c('DET', 'STO'), c('TP1', 'TP2')))


# EXPERIMENTS TO DISCARD !!
# DET: 20180730M
# STO: 20180809T, 20180813M, 20180823T
det <- det[lapply(det, function(i) i@date) != '20180730M']
sto <- sto[!lapply(sto, function(i) i@date) %in% c('20180809T', '20180813M', '20180823T')]

ref <- data.table(Frame = seq_len(21600))
cols <- viridis(3, begin = 0, end = 0.9, direction = -1)

#### ----- ACTIVITY AND INTERACTINONS ----- ####
det_ <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	I <- i@data[Crossings > 0, .(I = .N), by = 'Frame']
	x <- merge(i@data[, .N, by = 'Frame'], ref, all = TRUE, by = 'Frame')
	x <- merge(x, I, by = 'Frame', all = TRUE)[1:21600]
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x[is.na(I), 'I'] <- 0
	x

}), idcol = T)
det_[['cond']] <- 'DET'
det_food_times <- as.integer(apply(do.call('rbind', lapply(seq_along(det), function(i){

	x <- do.call('rbind', det[[i]]@food)
	range(x$t)	

})), 2, mean))

det_[['color']] <- 'a'

sto_ <- rbindlist(lapply(sto, function(i){

	setDT(i@data)
	I <- i@data[Crossings > 0, .(I = .N), by = 'Frame']
	x <- merge(i@data[, .N, by = 'Frame'], ref, all = TRUE, by = 'Frame')
	x <- merge(x, I, by = 'Frame', all = TRUE)[1:21600]
	x[is.na(N), 'N'] <- 0
	x[is.na(I), 'I'] <- 0
	x


}), idcol = T)
sto_[['cond']] <- 'STO'
sto_food <- lapply(sto, function(i) i@food)
sto_food_times <- as.integer(apply(do.call('rbind', lapply(seq_along(sto), function(i){

	x <- do.call('rbind', sto[[i]]@food)
	range(x$t)
})), 2, mean))

sto_[['color']] <- 'a'

nf_ <- rbindlist(lapply(nf, function(i){
	setDT(i@data)
	I <- i@data[Crossings > 0, .(I = .N), by = 'Frame']
	x <- merge(i@data[, .N, by = 'Frame'], ref, all = TRUE, by = 'Frame')
	x <- merge(x, I, by = 'Frame', all = TRUE)[1:21600]
	idx <- which.max(x[['N']] == 1)
	x[1:(idx-1), 'N'] <- NA
	x[is.na(N), 'N'] <- 0
	x[is.na(I), 'I'] <- 0
	x
}), idcol = T)
nf_[['cond']] <- 'NFD'

nf_[['color']] <- 'a'

fulldf <- rbind(det_, sto_, nf_)

A <- ggplot(data = fulldf, aes(Frame, N)) +
	geom_path(alpha = 0.5, color = 'grey80', aes(group = .id))+
	geom_path(data = fulldf[, .(N = mean(N, na.rm = T), col = color), 
				by = c('Frame', 'cond')],
		  aes(color = col), linewidth = 1.1)+
	geom_path(data = fulldf[, .(I = mean(I, na.rm = T)), 
				by = c('Frame', 'cond')][,
							 .(I = norm_range(cumsum(I),
							 		 a = 0,
							 		 b = 60 * max(cumsum(I))/2000), 
							   col = factor(5),
							   Frame = Frame),
							 by = 'cond'],
		  aes(color = col, y = I), linewidth = 1.1, show.legend = FALSE)+
	geom_vline(data = data.frame(x = c(det_food_times, sto_food_times, NA, NA), 
				     cond = c('DET', 'DET', 'STO', 'STO', 'NFD', 'NFD')),
		   aes(xintercept = x), linetype = 'dashed', linewidth = 1)+
	facet_wrap(~ factor(cond, levels = c('DET', 'STO', 'NFD'))) +
	scale_x_continuous('Time (min)', breaks = seq(0, 180*120, 50*120),
			   labels = seq(0, 180, 50))+
	scale_color_manual('', values = c('dodgerblue3', rep('black', 4)),
			   labels = c('Interactions', 'Occupancy'))+
	scale_y_continuous('Occupancy (number of ants)',
			   breaks = seq(0, 60, length.out = 7),
			   labels = seq(0, 60, length.out = 7),
			   
			   sec.axis = sec_axis(trans = ~.*1,
			   		    breaks = seq(0, 60, length.out = 9),
			   		    labels = seq(0, 2000, length.out = 9),
			   		    name = 'Interactions (cumulative sum)'))+
	theme(legend.position = c(0.9, 0.7),
	      legend.background = element_rect(fill = 'white', linewidth = 1,
	      				 linetype = 'solid', color = 'black'),
	      legend.title = element_blank(),
	      strip.text = element_text(size = 18, margin = margin(b = 5, t = 5)),
	      aspect.ratio = 0.75)

png('~/research/papers/SPIN_GLASSES/Figs_spinGlasses/Fig_Occupancy.png', width = 3600, height = 1000, res = 200)
A
dev.off()



#### ----- ACTIVITY AND INTERACTINONS ----- ####
det_phases <- list(p1 = vector('list', length(det)),
		   p2 = vector('list', length(det)),
		   p3 = vector('list', length(det)))

t0 <- Sys.time()
for(i in seq_along(det)){
	t <- range(do.call('rbind', det[[i]]@food)$t)
	det_phases$p1[[i]] <- mutual_info(det[[i]], t = 1:t[1], edges)
	det_phases$p2[[i]] <- mutual_info(det[[i]], t = (t[1]+1):t[2], edges)
	det_phases$p3[[i]] <- mutual_info(det[[i]], t = (t[2]+1):length(det[[i]]@N), edges)
}
Sys.time() - t0
beepr::beep(3)


zdet <- list(z1 = apply(do.call('rbind', det_phases$p1), 2, mean),
	     z2 = apply(do.call('rbind', det_phases$p2), 2, mean),
	     z3 = apply(do.call('rbind', det_phases$p3), 2, mean))


sto_phases <- list(p1 = vector('list', length(sto)),
		   p2 = vector('list', length(sto)),
		   p3 = vector('list', length(sto)))

t0 <- Sys.time()
for(i in seq_along(sto)){
	t <- range(do.call('rbind', sto[[i]]@food)$t)
	sto_phases$p1[[i]] <- mutual_info(sto[[i]], t = 1:t[1], edges)
	sto_phases$p2[[i]] <- mutual_info(sto[[i]], t = (t[1]+1):t[2], edges)
	sto_phases$p3[[i]] <- mutual_info(sto[[i]], t = (t[2]+1):length(sto[[i]]@N), edges)
}
Sys.time() - t0
beepr::beep(3)

t0 <- Sys.time()
for(i in seq_along(nf)){
	nf_phases <- mutual_info(nf[[i]], t = 1:21600, edges)
}
Sys.time() - t0
beepr::beep(3)

zdet <- list(z1 = apply(do.call('rbind', det_phases$p1), 2, mean),
	     z2 = apply(do.call('rbind', det_phases$p2), 2, mean),
	     z3 = apply(do.call('rbind', det_phases$p3), 2, mean))
zsto <- list(z1 = apply(do.call('rbind', sto_phases$p1), 2, mean),
	     z2 = apply(do.call('rbind', sto_phases$p2), 2, mean),
	     z3 = apply(do.call('rbind', sto_phases$p3), 2, mean))
znf <- list(z1 = nf_phases)

supervec <- rank(c(unlist(zdet), unlist(zsto), nf_phases, recursive = T))
ztotal <- list(d1 = supervec[1:(length(supervec)/7)],
	       d2 = supervec[(length(supervec)/7 + 1):(2*length(supervec)/7)],
	       d3 = supervec[(2*length(supervec)/7 + 1):(3*length(supervec)/7)],
	       s1 = supervec[(3*length(supervec)/7 + 1):(4*length(supervec)/7)],
	       s2 = supervec[(4*length(supervec)/7 + 1):(5*length(supervec)/7)],
	       s3 = supervec[(5*length(supervec)/7 + 1):(6*length(supervec)/7)],
	       nf = supervec[(6*length(supervec)/7 + 1):(length(supervec))])

ztotaldt <- rbindlist(list(data.table(mi = supervec[1:(length(supervec)/7)], exptype = 'DET', phase = 'Exploration'),
	       data.table(mi = supervec[(length(supervec)/7 + 1):(2*length(supervec)/7)],
	       	   exptype = 'DET', phase = 'Exploitation'),
	       data.table(mi = supervec[(2*length(supervec)/7 + 1):(3*length(supervec)/7)],
	       	   exptype = 'DET', phase = 'Post-exploitation'),
	       data.table(mi = supervec[(3*length(supervec)/7 + 1):(4*length(supervec)/7)],
	       	   exptype = 'STO', phase = 'Exploration'),
	       data.table(mi = supervec[(4*length(supervec)/7 + 1):(5*length(supervec)/7)],
	       	   exptype = 'STO', phase = 'Exploitation'),
	       data.table(mi = supervec[(5*length(supervec)/7 + 1):(6*length(supervec)/7)],
	       	   exptype = 'STO', phase = 'Post-exploitation'),
	       data.table(mi = supervec[(6*length(supervec)/7 + 1):(length(supervec))],
	       	   exptype = 'NFD', phase = 'Exploration')))


TOGETHER_ <- list(ggplot()+theme_void(), draw_hexagons(z = ztotal$nf,
						       size = 3, lineend = 'round',
						       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round',
						       		    add = geom_foodpatches(
						       		    	food = det[[1]]@food, fill = 'mediumpurple1', color = 'white')
						       )) +
		  	ylab('NFD')+ 
		  	theme(aspect.ratio = 0.5, legend.position = 'none')+
		  	scico::scale_color_scico(palette = 'bilbao', direction = 1)+
		  	geom_point(data = data.frame(x = hex$x[hex$node == 634],
		  				     y = hex$y[hex$node == 634] - 35),
		  		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
		  draw_hexagons(z = apply(rbind(do.call('rbind', zdet), do.call('rbind', zsto), do.call('rbind', znf)), 2, max),
		  	      size = 3, lineend = 'round',
		  	      add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', 
		  	      		    add = geom_foodpatches(
		  	      		    	food = det[[1]]@food, fill = 'mediumpurple1', color = 'white')
		  	      )) +
		  	theme(aspect.ratio = 0.5, legend.position = 'left',
		  	      legend.key.width = unit(1, 'cm'),
		  	      legend.direction = 'horizontal',
		  	      legend.title.align = 0.5,
		  	      legend.background = element_blank(),
		  	      legend.box.background = element_rect(fill = NA, color = 'black', linewidth = 1),
		  	      legend.box.margin = unit(c(5, 10, 5, 15), 'pt'),
		  	      legend.box.spacing = margin(r = 5, unit = 'cm'),
		  	      plot.margin = margin(l = 0, r = -600, unit = 'pt'))+
		  	scico::scale_color_scico('Mutual information', palette = 'bilbao',
		  				 direction = 1, limits = c(0, NA))+
		  	geom_point(data = data.frame(x = hex$x[hex$node == 634],
		  				     y = hex$y[hex$node == 634] - 35),
		  		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3)+
		  	guides(color = guide_colorbar(title.position = 'top', )),
		  draw_hexagons(z = ztotal$d1,
			   size = 3, lineend = 'round',
			   add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round',
			   		    add = geom_foodpatches(
			   		    	food = det[[1]]@food, fill = 'mediumpurple1', color = 'white')
			   )) +
	ylab('DET')+ scale_x_continuous('Exploration', position = 'top')+
	theme(aspect.ratio = 0.5, legend.position = 'none')+
	scico::scale_color_scico(palette = 'bilbao', direction = 1)+
	geom_point(data = data.frame(x = hex$x[hex$node == 634],
				     y = hex$y[hex$node == 634] - 35),
		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),

draw_hexagons(z = ztotal$d2,
	      size = 3, lineend = 'round',
	      add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add =
	      		    	geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1', color = 'white'))) +
	theme(aspect.ratio = 0.5, legend.position = 'none')+
	scale_x_continuous('Exploitation', position = 'top')+
	scico::scale_color_scico(palette = 'bilbao', direction = 1)+
	geom_point(data = data.frame(x = hex$x[hex$node == 634],
				     y = hex$y[hex$node == 634] - 35),
		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
draw_hexagons(z = ztotal$d3,
	      size = 3, lineend = 'round',
	      add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add =
	      		    	geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1',
	      		    			 color = 'white'))) +
	scale_x_continuous('Relaxation', position = 'top')+
	theme(aspect.ratio = 0.5, legend.position = 'none') +
	scico::scale_color_scico(palette = 'bilbao', direction = 1)+
	geom_point(data = data.frame(x = hex$x[hex$node == 634],
				     y = hex$y[hex$node == 634] - 35),
		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
draw_hexagons(z = ztotal$s1, size = 3, lineend = 'round',
	      add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add =
	      		    	geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
	theme(aspect.ratio = 0.5, legend.position = 'none') +
	scico::scale_color_scico(palette = 'bilbao')+
	geom_point(data = data.frame(x = hex$x[hex$node == 634],
				     y = hex$y[hex$node == 634] - 35),
		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3)+
	ylab('STO'),
draw_hexagons(z = ztotal$s2, size = 3, lineend = 'round',
	      add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add =
	      		    	geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
	theme(aspect.ratio = 0.5, legend.position = 'none') +
	scico::scale_color_scico(palette = 'bilbao')+
	geom_point(data = data.frame(x = hex$x[hex$node == 634],
				     y = hex$y[hex$node == 634] - 35),
		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
draw_hexagons(z = ztotal$s3, size = 3, lineend = 'round',
	      add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add =
	      		    	geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
	theme(aspect.ratio = 0.5, legend.position = 'none') +
	scico::scale_color_scico(palette = 'bilbao')+
	geom_point(data = data.frame(x = hex$x[hex$node == 634],
				     y = hex$y[hex$node == 634] - 35),
		   aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3)
)


png('~/research/papers/SPIN_GLASSES/Figs_spinGlasses/Fig_MI.png', width = 8000, height = 4000, res = 400)
ggarrange(plotlist = TOGETHER_, nrow =3, ncol = 3)
dev.off()