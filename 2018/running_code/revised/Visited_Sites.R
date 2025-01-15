ldet <- lapply(det, function(i){
	x <- setDT(i@data)[, .SD[1], by = 'node'][order(Frame), Frame]
	y <- seq_along(x)
	data.table(t = x - i@data[, min(Frame)], visited_sites = y)
})




lnf <- lapply(nf[-c(1, 2)], function(i){
	x <- setDT(i@data)[, .SD[1], by = 'node'][order(Frame), Frame]
	y <- seq_along(x)
	data.table(t = x - i@data[, min(Frame)], visited_sites = y)
})


ggplot(data = rbindlist(list(det = rbindlist(ldet)[visited_sites < 414,
		.(t = sum(t) / .N), by = 'visited_sites'],

nf = rbindlist(lnf)[visited_sites < 414, .(t = sum(t) / .N), by = 'visited_sites']),
idcol = 'exp'), aes(t / 120, visited_sites, color = exp))+
	geom_point() + geom_path()


ggplot(data = rbindlist(list(det = rbindlist(ldet),
			     
			     nf = rbindlist(lnf)),
			idcol = 'exp'), aes(t / 120, color = exp, fill = exp))+
	geom_density(linewidth = 1, alpha = 0.5, bw = 5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3'))+
	scale_x_continuous('Time (min)') + ylab('Density')+
	geom_vline(linetype = 2, linewidth = 1, 
		   xintercept = median(sapply(det, function(i){
		   	min(rbindlist(i@food)[['t']]) / 120
		   })))
