source('~/research/gits/AnTracks/src/Experiment.R')

load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/nf.RData')


## ANALYSIS 1 : FPT ####

foodnodes <- unlist(lapply(det[1], function(i){
	p <- vapply(i@food[1], function(x){
		l <- pdist(as.matrix(x[, 1:2]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		r <- range(l)
		avg <- mean(l)
		c(r, avg)
	}, numeric(3))
	data.table(mind = p[1], maxd = p[2], avgd = p[3])
}))

det_times <- unlist(lapply(det, function(i){
	f <- min(rbindlist(i@food)[['t']])
	f - i@data[1, Frame]
}))
det_times <- c(range(det_times), mean(det_times))



R_exp <- seq(2.5, 20, 1.25) * 50


det_result <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)
	mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_exp, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'det')
}))

nf_result <- rbindlist(lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)
	mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_exp, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'nf')
}))



panel_a <- ggplot(data = rbindlist(list(det = det_result, nf = nf_result))[R %in% seq(2.5*50, 20*50, 50*2.5)],
		  aes(factor(R/10), t / 120, fill = factor(rho)))+
	geom_boxplot(outlier.shape = NA, show.legend = FALSE, alpha = 0.6)+
	geom_jitter(width = 0.15, alpha = 0.4, size = 3, show.legend = FALSE)+
	facet_wrap(~ factor(rho, labels = c('DET', 'NFD')))+
	scale_y_continuous('First passage time (min)', breaks = seq(0, 80, 15))+
	scale_x_discrete('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+ ggtitle('A')+
	theme(plot.title = element_text(size = 22))



fit_DET <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
					data = det_result[, .(t = median(t)), by = 'R'],
					supp_errors = 'Y', iter = 3e3,
					start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
					start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))
fit_NFD <- lm(t ~ poly(R,2), data = nf_result[, .(t = median(t)), by = 'R'])

R_exp_hd <- seq(2.5, 20, length.out = 150)*50

data_fits <- rbind(data.table(R = R_exp_hd, t = SSlogis(R_exp_hd, coef(fit_DET)[['Asym']],
							coef(fit_DET)[['xmid']], coef(fit_DET)[['scal']]),
							rho = 'det'),
		   data.table(R = R_exp_hd, t = predict(fit_NFD,
		   				     newdata = data.frame(R = R_exp_hd)),
		   	   rho = 'nf'))


panel_b <- ggplot(data = rbindlist(list(det = det_result, nf = nf_result))[, .(t = median(t)),
									   by = c('R', 'rho')],
		  aes(R/10, t / 120, fill = factor(rho), color = factor(rho)))+
	geom_polygon(fill = 'grey90', color = NA,
		     data = data.frame(x = foodnodes[c(1, 1, 2, 2)]/10,
		     		  y = c(-Inf, Inf, Inf, -Inf)),
		     aes(x = x, y = y)) +
	geom_polygon(fill = 'grey90', color = NA,
		     data = data.frame(x = c(-Inf, Inf, Inf, -Inf),
		     		  y = det_times[c(1, 1, 2, 2)]/120),
		     aes(x = x, y = y)) +
	geom_vline(xintercept = foodnodes[3]/10, linewidth = 1, linetype = 2, color = 'grey35')+
	geom_hline(yintercept = det_times[3]/120, linewidth = 1, linetype = 2, color = 'grey35')+
	scale_y_continuous('Median first passage time (min)', breaks = seq(0, 80, 10))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'),
			  labels = c('DET', 'NFD'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3'),
			   labels = c('DET', 'NFD'))+
	geom_path(data = data_fits,
		  size = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, shape = 21, color = 'black', alpha = 0.6)+
	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22)) +
	ggtitle('B')


png(filename = 'C:/Users/pfern/research/gits/AnTracks/plots/Fig_SI_1.png', 4000, 4000, res = 400)
ggarrange(panel_a, panel_b, ncol = 1, nrow = 2)
dev.off()
## ANALYSIS 2 : EXPLORATION OF BORDER NODES ####

library(latex2exp)

nfd <- lapply(nf, function(i){
	i@food <- lapply(det[[1]]@food, function(x) x[, 1:2])
	i
})
nfd <- lapply(nfd, food_detection)

#
# ggplot(data = rbindlist(list(det = rbindlist(lapply(det, function(i){
#         data.table(mint = min(rbindlist(i@food)[['t']]) - i@data[1, Frame],
#                    exp = 'det')
# })), nfd = rbindlist(lapply(nfd, function(i){
#         data.table(mint = min(rbindlist(i@food)[['t']]) - i@data[1, Frame],
#                    exp = 'nfd')
# })))), aes(exp, fill = exp, 1/mint)) + geom_boxplot()

# ggplot(data = rbindlist(list(det = rbindlist(lapply(det, function(i){
# 	data.table(mint = min(rbindlist(i@food)[['t']]) - i@data[1, Frame],
# 		   exp = 'det')
# })), nfd = rbindlist(lapply(nfd[-c(1,2)], function(i){
# 	data.table(mint = min(rbindlist(i@food)[['t']]) - i@data[1, Frame],
# 		   exp = 'nfd')
# })))), aes(factor(exp), fill = exp, mint / 120)) +
# 	geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = FALSE)+
# 	scale_fill_manual('', values = c('mediumpurple', 'brown4'),
# 			  labels = c('DET', 'NFD')) +
# 	ylab('First passage time (min)') +
# 	scale_x_discrete('', labels = c('DET', 'NFD')) +
# 	geom_jitter(alpha = 0.25, size = 4, width = 0.25, show.legend = FALSE)



# ggplot(data = rbindlist(list(det = rbindlist(lapply(det, function(i){
# 	data.table(maxt = max(rbindlist(i@food)[['t']]) - i@data[1, Frame],
# 		   exp = 'det')
# })), nfd = rbindlist(lapply(nfd[-c(1,2)], function(i){
# 	data.table(maxt = max(rbindlist(i@food)[['t']]) - i@data[1, Frame],
# 		   exp = 'nfd')
# })))), aes(factor(exp), fill = exp, 1/(maxt/2))) +
# 	geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = FALSE)+
# 	scale_fill_manual('', values = c('mediumpurple', 'brown4'),
# 			  labels = c('DET', 'NFD')) +
# 	scale_y_continuous('Last detection time') +
# 	scale_x_discrete('', labels = c('DET', 'NFD')) +
# 	geom_jitter(alpha = 0.25, size = 4, width = 0.25, show.legend = FALSE)


range_exps <- rbindlist(list(det = rbindlist(lapply(det, function(i){
	t0 <- i@data[1, Frame]
	r <- (range(rbindlist(i@food)[['t']])-t0) /2
	data.table(mint = r[1], maxt = r[2], exp = 'det')
})), nfd = rbindlist(lapply(nfd[-c(1,2)], function(i){
	t0 <- i@data[1, Frame]
	r <- (range(rbindlist(i@food)[['t']]) - t0)/2
	data.table(mint = r[1], maxt = r[2], exp = 'nfd')
}))))

eff_exps <- rbindlist(list(det = rbindlist(lapply(det, get_eff)),
			   nfd = rbindlist(lapply(nfd, get_eff))), idcol = TRUE)

dt <- cbind(eff_exps[-c(11, 12)], mint = range_exps[['mint']], maxt = range_exps[['maxt']])


png('C:/Users/pfern/research/gits/AnTracks/plots/boxplots_effs_DET_NFD.png', 6000, 4000, res = 400)
ggplot(data = melt(dt, id.vars = '.id'), aes(factor(.id), fill = .id, 1/(value))) +
	geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = FALSE)+
	scale_fill_manual('', values = c('mediumpurple', 'brown4'),
			  labels = c('DET', 'NFD')) +
	scale_y_continuous(TeX('Exploration efficiency ($s^{-1}$)')) +
	scale_x_discrete('', labels = c('DET', 'NFD')) +
	geom_jitter(alpha = 0.25, size = 4, width = 0.25, show.legend = FALSE) +
	facet_wrap(~factor(variable, levels = c('mint', 'maxt',
						'tp1', 'tp2'), labels = c('First passage time',
									  'Last detection time',
									  'Exploration time',
									  'Exploitation time')), scales = 'free')

dev.off()

# ggplot(data = eff_exps[-c(11, 12)], aes(factor(.id), fill = .id, 1/(tp1))) +
# 	
# 	geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = FALSE)+
# 	scale_fill_manual('', values = c('mediumpurple', 'brown4'),
# 			  labels = c('DET', 'NFD')) +
# 	scale_y_continuous(TeX('Exploration efficiency ($s^{-1}$)')) +
# 	scale_x_discrete('', labels = c('DET', 'NFD')) +
# 	geom_jitter(alpha = 0.25, size = 4, width = 0.25, show.legend = FALSE)
# 
# ggplot(data = eff_exps[-c(11, 12)], aes(factor(.id), fill = .id, 1/(tp2))) +
# 	geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = FALSE)+
# 	scale_fill_manual('', values = c('mediumpurple', 'brown4'),
# 			  labels = c('DET', 'NFD')) +
# 	scale_y_continuous(TeX('Exploitation efficiency ($s^{-1}$)')) +
# 	scale_x_discrete('', labels = c('DET', 'NFD')) +
# 	geom_jitter(alpha = 0.25, size = 4, width = 0.25, show.legend = FALSE)




pos_det <- data.table(merge(hex, rbindlist(lapply(det, function(i){
	i@data[, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node'))
# pos_det[N < quantile(N, p = 0.75), 'N'] <- 0
pos_det[N > quantile(N, p = 0.75), 'N'] <- 0

draw_nodes(pos_det,
	   add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    		    add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))),
	   z = rank(pos_det[['N']]), size = 4) +
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = c(1, 620),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))







sbst_nodes <- data.table(rbind(hex[hex$y > 1000 & hex$y < 1100, c('x', 'y', 'node') ],
			       hex[hex$y > 1900, c('x', 'y', 'node') ],
			       hex[hex$x < 120 & hex$y > 1000, c('x', 'y', 'node') ],
			       hex[hex$y > 1000 & hex$x > 1850, c('x', 'y', 'node') ]))[, node]


sbst_nodes <- hex[hex$y > 1000 & hex$node %in% sbst_nodes, ]
sbst_nodes <- sbst_nodes[!sbst_nodes$node %in% pos_det[N == 0, node], ]
sbst_nodes <- sbst_nodes[-which(pdist(as.matrix(sbst_nodes[,c('x', 'y')]), 
				      as.matrix(hex[hex$node == 634, c('x', 'y')])) < 200), ]


draw_nodes(merge(hex, sbst_nodes)
	   , add = draw_hexagons())



all_det <- rbindlist(lapply(det, function(i){
	d <- i@data
	d[['Frame']] <- d[['Frame']] - d[1, Frame] + 1
	d
}), idcol = TRUE)[, c('Frame', 'node', '.id')]
all_det[['cat']] <- ifelse(all_det[['node']] %in% sbst_nodes[['node']], 1, 2)

ratio_det <- all_det[, .(ratio = sum(cat == 1)/sum(cat == 2),
			 proportion1 = sum(cat == 1) / sum(cat < 3),
			 cat1 = sum(cat == 1)), by = c('Frame', '.id')]

all_nfd <- rbindlist(lapply(nfd[-c(1, 2)], function(i){
	d <- i@data
	d[['Frame']] <- d[['Frame']] - d[1, Frame] + 1
	d
}), idcol = TRUE)[, c('Frame', 'node', '.id')]
all_nfd[['cat']] <- ifelse(all_nfd[['node']] %in% sbst_nodes[['node']], 1, 2)

ratio_nfd <- all_nfd[, .(ratio = sum(cat == 1)/sum(cat == 2),
			 proportion1 = sum(cat == 1) / sum(cat < 3),
			 cat1 = sum(cat == 1)), by = c('Frame', '.id')]

ggplot(data= ratio_det[is.finite(ratio)], aes(Frame, ratio))+ geom_line() +
	facet_wrap(~.id)
ggplot(data= ratio_det[is.finite(ratio), .(ratio = mean(ratio)), by = 'Frame'], aes(Frame, ratio))+ geom_line()
ggplot(data= ratio_nfd[is.finite(ratio), .(ratio = mean(ratio)), by = 'Frame'], aes(Frame, ratio))+ geom_line()

mixed_ratio <- rbindlist(list(det = ratio_det[is.finite(ratio), .(ratio = mean(ratio)), by = 'Frame'],
			 nfd = ratio_nfd[is.finite(ratio), .(ratio = mean(ratio)), by = 'Frame']), idcol = TRUE)
mixed_prop1 <- rbindlist(list(det = ratio_det[, .(N = mean(proportion1)), by = 'Frame'],
			     nfd = ratio_nfd[, .(N = mean(proportion1)), by = 'Frame']), idcol = TRUE)
mixed_cat1 <- rbindlist(list(det = ratio_det[, .(N = mean(cat1)), by = 'Frame'],
			      nfd = ratio_nfd[, .(N = mean(cat1)), by = 'Frame']), idcol = TRUE)

det_foods <- c(mean(dt[.id == 'det', mint*2]), mean(dt[.id == 'det', maxt*2]))

png('~/research/gits/AnTracks/plots/ocupation_of_border_nodes.png', 4000, 6000, res = 500)
ggarrange(
ggplot(data = mixed_ratio, aes(Frame, ratio, color = .id)) + geom_line(linewidth = 1)+
	scale_color_manual('', labels = c('DET', 'NFD'), values = c('mediumpurple', 'gold3'))+
	scale_x_continuous('Time (min)', breaks = seq(0, 21600, 50*120),
			   labels = seq(0, 180, 50))+
	scale_y_continuous('Ratio')+
	geom_vline(xintercept = det_foods, linetype = 2, color = 'grey35'),
ggplot(data = mixed_cat1, aes(Frame, N, color = .id)) + geom_line(linewidth = 1)+
	scale_color_manual('', labels = c('DET', 'NFD'), values = c('mediumpurple', 'gold3'))+
	scale_x_continuous('Time (min)', breaks = seq(0, 21600, 50*120),
			   labels = seq(0, 180, 50))+
	scale_y_continuous('Number of ants')+
	geom_vline(xintercept = det_foods, linetype = 2, color = 'grey35'),
ggplot(data = mixed_prop1, aes(Frame, N, color = .id)) + geom_line(linewidth = 1)+
	scale_color_manual('', labels = c('DET', 'NFD'), values = c('mediumpurple', 'gold3'))+
	scale_x_continuous('Time (min)', breaks = seq(0, 21600, 50*120),
			   labels = seq(0, 180, 50))+
	scale_y_continuous('Proportion')+
	geom_vline(xintercept = det_foods, linetype = 2, color = 'grey35')
, common.legend = TRUE, ncol = 1, nrow = 3, labels = c('Ocupation in border nodes', NA, NA),
font.label = list(face = 'plain'))
dev.off()

# NbyFrame <- megadata[, .N, by = c('Frame', '.id')]
# 
# 
# NbyCat <- megadata[, .N, by = c('Frame', '.id', 'node', 'cat')]

proportion <- NbyCat[, .(p = sum(N)), by = c('Frame', '.id', 'cat')]


ggplot(data= ratio[.id == 1], aes(Frame, ratio))+ geom_line()
ggplot(data= ratio[is.finite(ratio)], aes(Frame, ratio))+ geom_line() +
	facet_wrap(~.id)

ggplot(data= proportion[.id == 1], aes(Frame, p, color = factor(cat)))+ geom_line()

ggplot(data= proportion, aes(Frame, p, color = factor(cat)))+ geom_line() +
	facet_wrap(~.id)

## ANALYSIS 3 : DISTANCES TO FOOD, DISTANCES TO NEST ####
nest_ref <- as.matrix(hex[hex$node == 634, c('x', 'y')])
p1_ref <- t(as.matrix(apply(det[[1]]@food[[1]][, 1:2], 2, mean)))
p2_ref <- t(as.matrix(apply(det[[1]]@food[[2]][, 1:2], 2, mean)))
dnest_p <- pdist(nest_ref, p1_ref)


P <- nest_ref + 1
d1 <- pdist(P, p1_ref)
d2 <- pdist(P, p2_ref)
d3 <- pdist(P, nest_ref)
(d1/d3) - (d2/d3)



distances_nest_det <- rbindlist(lapply(det, function(i){
	data.table(Frame = i@data[['Frame']], d = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), nest_ref))
}), idcol = TRUE)
colnames(distances_nest_det) <- c('exp', 'Frame', 'd')
ggplot(data = distances_nest_det[, .(d = mean(d)), by = c('exp', 'Frame')], aes(Frame, d)) +
	geom_point() + geom_line() + 
	facet_wrap(~ exp)

distances_patch_det <- rbindlist(lapply(det, function(i){
	d1 = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), p1_ref)
	d2 = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), p2_ref)
	d = apply(cbind(d1, d2), 1, min)
	data.table(Frame = i@data[['Frame']], d = d)
}), idcol = TRUE)

trilateration_det <- rbindlist(lapply(det, function(i){
	d1 = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), p1_ref)
	d2 = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), p2_ref)
	d3 = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), nest_ref)
	d = (d1 / d3) - (d2 / d3)
	# d = apply(cbind(d1, d2), 1, min)
	data.table(Frame = i@data[['Frame']] - min(i@data[['Frame']]) +1, d = d)
}), idcol = TRUE)
colnames(trilateration_det) <- c('exp', 'Frame', 'd')
ggplot(data = trilateration_det[exp == 1, .(d = mean(d)), by = 'Frame'], aes(Frame, d)) + 
	geom_path() + geom_point() +	
	geom_vline(xintercept = c(mean(dt[.id == 'det', mint][1]*2), mean(dt[.id == 'det', maxt][1]*2)),
						linetype = 2, color = 'grey35')
ggplot(data = trilateration_det[exp == 1, .(d = median(d)), by = 'Frame'], aes(Frame, d)) + 
	geom_path() + geom_point() +	
	geom_vline(xintercept = c(mean(dt[.id == 'det', mint][1]*2), mean(dt[.id == 'det', maxt][1]*2)),
		   linetype = 2, color = 'grey35')
	

ggplot(data = det[[1]]@data[, .N, by = 'Frame'], aes(Frame-619, N)) + 
	geom_path()+	
	geom_vline(xintercept = c(mean(dt[.id == 'det', mint][1]*2), mean(dt[.id == 'det', maxt][1]*2)),
		   linetype = 2, color = 'grey35')
	



ggplot(data = distances_patch_det[, .(d = mean(d)), by = c('exp', 'Frame')], aes(Frame, d)) +
	geom_point() + geom_line() + 
	facet_wrap(~ exp)


distances_nest_nf <- rbindlist(lapply(nf[-c(1,2)], function(i){
	data.table(Frame = i@data[['Frame']], d = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), nest_ref))
}), idcol = TRUE)
colnames(distances_nest_nf) <- c('exp', 'Frame', 'd')
ggplot(data = distances_nest_nf[, .(d = mean(d)), by = c('exp', 'Frame')], aes(Frame, d)) +
	geom_point() + geom_line() + 
	facet_wrap(~ exp)

distances_patch_nf <- rbindlist(lapply(nf[-c(1,2)], function(i){
	d1 = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), p1_ref)
	d2 = pdist(as.matrix(i@data[, c('Xmm', 'Ymm')]), p2_ref)
	d = apply(cbind(d1, d2), 1, min)
	data.table(Frame = i@data[['Frame']], d = d)
}), idcol = TRUE)
colnames(distances_patch_nf) <- c('exp', 'Frame', 'd')
ggplot(data = distances_patch_nf[, .(d = mean(d)), by = c('exp', 'Frame')], aes(Frame, d)) +
	geom_point() + geom_line() + 
	facet_wrap(~ exp)

	

