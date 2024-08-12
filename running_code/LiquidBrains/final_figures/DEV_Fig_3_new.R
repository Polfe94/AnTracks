source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/results/mov_rho_R.RData')
source('~/research/gits/AnTracks/src/fit_functions.R')

library(arrow)

load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/nf.RData')



path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('default_\\d{1,2}', files))))

path_1 <- '/home/polfer/research/gits/AutomatAnts/results/2024/SUPPLEMENTARY/rho_1/'
files_1 <- list.files(path_1)
files_1 <- unique(unlist(regmatches(files_1,
				    gregexpr('rho_1.001_epsilon_1.001_\\d{1,2}', files_1))))

R_sim <- seq(2.5, 20, 1.25)
R_exp <- R_sim * 50


## couple minutes to process
result_sim <- lapply(files, function(f){
	print(which(files == f))
	do.call('gc', args = list(verbose = FALSE))
	data <- setDT(read_parquet(paste0(path, f, '_data.parquet')))
	pos <- data[, .(node = unique(parse_nodes(pos))), by = 'Frame']
	

	mdata <- merge(pos, hex_sim[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex_sim[hex_sim$node == '(0, 22)', c('x', 'y')]), 
		       as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_sim, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'sim')
})

result_sim_1 <- lapply(files_1, function(f){
	print(which(files_1 == f))
	do.call('gc', args = list(verbose = FALSE))
	data <- setDT(read_parquet(paste0(path_1, f, '_data.parquet')))
	pos <- data[, .(node = unique(parse_nodes(pos))), by = 'Frame']
	
	
	mdata <- merge(pos, hex_sim[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex_sim[hex_sim$node == '(0, 22)', c('x', 'y')]), 
		       as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_sim, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'sim_1')
})

estimated <- rbindlist(result_sim)
estimated_1 <- rbindlist(result_sim_1)

ggplot(data = rbind(rbindlist(list(det = det_result, nf = nf_result))[, .(t = mean(t)),
								      by = c('R', 'rho')],
		    estimated[, .(t = mean(t)), by = c('R', 'rho')]),
       aes(R/10, t / 120, fill = factor(rho, levels = c('det', 'nf', 'sim')),
           color = factor(rho, levels = c('det', 'nf', 'sim')),
           shape = factor(rho, levels = c('det', 'nf', 'sim'))))+
	
	
	
	scale_y_continuous('Mean first passage time (min)', breaks = seq(0, 30, 5))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'),
			  breaks = c('det', 'nf', 'sim'),
			  labels = c('DET', 'NFD', 'Model'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4'),
			   labels = c('DET', 'NFD', 'Model'))+
	scale_linetype_manual('', values = c(1, 1, 2),
			      breaks = c('det', 'nf', 'sim'),
			      labels = c('DET', 'NFD', 'Model'))+
	scale_shape_manual('', values = c(22, 22, 24),
			   breaks = c('det', 'nf', 'sim'),
			   labels = c('DET', 'NFD', 'Model'))+
	# geom_path(data = data_fits_combined, 
	# 	  aes(linetype = factor(rho, levels = c('det', 'nf', 'LR', 'SR'))),
	# 	  linewidth = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, color = 'black', alpha = 0.8)+
	
	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22),
	      aspect.ratio = 0.5) 


## combining with rho = 1
ggplot(data = rbind(rbindlist(list(det = det_result, nf = nf_result))[, .(t = mean(t)),
										 by = c('R', 'rho')],
			       estimated[, .(t = mean(t)), by = c('R', 'rho')],
		    estimated_1[, .(t = mean(t)), by = c('R', 'rho')]),
		  aes(R/10, t / 120, fill = factor(rho, levels = c('det', 'nf', 'sim', 'sim_1')),
		      color = factor(rho, levels = c('det', 'nf', 'sim', 'sim_1')),
		      shape = factor(rho, levels = c('det', 'nf', 'sim', 'sim_1'))))+



	scale_y_continuous('Mean first passage time (min)', breaks = seq(0, 30, 5))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'brown4'),
			  breaks = c('det', 'nf', 'sim', 'sim_1'),
			  labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4','brown4'),
			   labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	scale_linetype_manual('', values = c(1, 1, 2, 2),
			      breaks = c('det', 'nf', 'sim', 'sim_1'),
			      labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	scale_shape_manual('', values = c(22, 22, 24, 21),
			   breaks = c('det', 'nf', 'sim', 'sim_1'),
			   labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	# geom_path(data = data_fits_combined, 
	# 	  aes(linetype = factor(rho, levels = c('det', 'nf', 'LR', 'SR'))),
	# 	  linewidth = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, color = 'black', alpha = 0.8)+
	geom_path(aes(group = factor(rho)))+

	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22),
	      aspect.ratio = 0.5) 





fit_DET <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
					data = det_result[, .(t = mean(t)), by = 'R'],
					supp_errors = 'Y', iter = 3e3,
					start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
					start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))
fit_NFD <- lm(t ~ poly(R,2), data = nf_result[, .(t = mean(t)), by = 'R'])


R_exp_hd <- seq(2.5, 20, length.out = 150)*50

data_fits <- rbind(data.table(R = R_exp_hd, t = SSlogis(R_exp_hd, coef(fit_DET)[['Asym']],
							coef(fit_DET)[['xmid']], coef(fit_DET)[['scal']]),
							rho = 'det'),
		   data.table(R = R_exp_hd, t = predict(fit_NFD,
		   				     newdata = data.frame(R = R_exp_hd)),
		   	   rho = 'nf'))








# estimated <- mov_rho_R[rho == 0 | rho == 1, .(t = mean(N)), by = c('R', 'rho')]
# estimated[['rho']] <- ifelse(estimated[['rho']] == 0, 'SR', 'LR')
# estimated[['t']] <- estimated[['t']] * (50/10) * 2
# estimated[['R']] <- estimated[['R']] *50



fit_default <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
					    data = estimated[, .(t = mean(t)), by = 'R'],
					    supp_errors = 'Y', iter = 3e3,
					    start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
					    start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))
fit_1 <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
				      data = estimated_1[, .(t = mean(t)), by = 'R'],
				      supp_errors = 'Y', iter = 3e3,
				      start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
				      start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))


data_fits <- rbind(data.table(R = R_exp_hd, t = SSlogis(R_exp_hd, coef(fit_DET)[['Asym']],
							coef(fit_DET)[['xmid']], coef(fit_DET)[['scal']]),
							rho = 'det'),
		   data.table(R = R_exp_hd, t = predict(fit_NFD,
		   				     newdata = data.frame(R = R_exp_hd)),
		   	   rho = 'nf'))

data_fits_combined <- rbind(data_fits, data.table(R = sort(unique(estimated[['R']])),
						  t = SSlogis(sort(unique(estimated[['R']])), 
						  	  coef(fit_default)[['Asym']],
						  	  coef(fit_default)[['xmid']],
						  	  coef(fit_default)[['scal']]),
						  rho = 'sim'),
			    data.table(R = sort(unique(estimated_1[['R']])),
			    	   t = SSlogis(sort(unique(estimated_1[['R']])), 
			    	   	    coef(fit_1)[['Asym']],
			    	   	    coef(fit_1)[['xmid']],
			    	   	    coef(fit_1)[['scal']]),
			    	   rho = 'sim_1'))



ggplot(data = rbind(rbindlist(list(det = det_result, nf = nf_result))[, .(t = mean(t)),
								      by = c('R', 'rho')],
		    estimated[, .(t = mean(t)), by = c('R', 'rho')],
		    estimated_1[, .(t = mean(t)), by = c('R', 'rho')]),
       aes(R/10, t / 120, fill = factor(rho, levels = c('det', 'nf', 'sim', 'sim_1')),
           color = factor(rho, levels = c('det', 'nf', 'sim', 'sim_1')),
           shape = factor(rho, levels = c('det', 'nf', 'sim', 'sim_1'))))+
	
	
	
	scale_y_continuous('Mean first passage time (min)', breaks = seq(0, 30, 5))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'brown4'),
			  breaks = c('det', 'nf', 'sim', 'sim_1'),
			  labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4','brown4'),
			   labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	scale_linetype_manual('', values = c(1, 1, 2, 2),
			      breaks = c('det', 'nf', 'sim', 'sim_1'),
			      labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	scale_shape_manual('', values = c(22, 22, 24, 21),
			   breaks = c('det', 'nf', 'sim', 'sim_1'),
			   labels = c('DET', 'NFD', 'Model default', 'Model LR'))+
	# geom_path(data = data_fits_combined, 
	# 	  aes(linetype = factor(rho, levels = c('det', 'nf', 'LR', 'SR'))),
	# 	  linewidth = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, color = 'black', alpha = 0.8)+
	
	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22),
	      aspect.ratio = 0.5) +
geom_path(data = data_fits_combined, 
	  aes(linetype = factor(rho, levels = c('det', 'nf', 'sim', 'sim_1'))),
	  linewidth = 1.75, alpha = 0.6, show.legend = FALSE)
