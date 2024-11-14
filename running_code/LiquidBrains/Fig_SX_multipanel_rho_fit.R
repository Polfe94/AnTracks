#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
# library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')


## EXPERIMENTS
foods_det <- rbindlist(lapply(det[-7], function(i){
	a <- rbindlist(i@food)[['t']]
	data.table(mint = min(a)/2, maxt = max(a)/2)
}))[, lapply(.SD, mean)] # <--- median !

det_avg <- rbindlist(lapply(det[-7], function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	# if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	# }
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']

det_N <- rbindlist(lapply(det[-7], function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	# if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	# }
	x[, .(N = movingAverage(N, t = 60)[1:695], Frame = seq(31, 21545, 31))]
}), idcol = TRUE)

process_data <- function(path, 
			 vars = c('Frame', 'N')){
	ref <- data.table(Frame = 1:21600)
	t <- movingAverage(ref[['Frame']], 60, 0)
	
	files <- list.files(path)
	files <- files[!grepl('food', files) & !grepl('position', files) &
		       	!grepl('data', files) & !grepl('keys', files)]
	
	a <- lapply(files, function(i){
		
		y <- merge(ref, data.table(read_parquet(paste0(path, i)))[, ..vars], all = TRUE, by = 'Frame')
		y[['N']] <- fillVec(y[['N']])
		n <- movingAverage(y[['N']], 60L, 0L)
		z <- zscore(n)
		data.table(Frame = t, N = n, Z = z)
		
	})
	names(a) <- files
	a
}

path <- '/home/polfer/research/gits/AnTracks/results/sims/fig_S1/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('rho_\\d{1}.\\d{1,3}_epsilon_1.001_\\d{1,2}', files))))

#### POPULATION DYNAMICS ####

data <- process_data(path)

full_data <- rbindlist(lapply(files, function(i){
	cbind(data.table(read_parquet(paste0(path, i,'.parquet'))),
	      rho = round(as.numeric(strsplit(unlist(regmatches(i, gregexpr('rho_\\d{1}.\\d{1,3}', i))), split = '_')[[1]][2]), 2))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[Frame < 175*120, .(N = mean(N)), by = c('rho','Frame')] ## filtered

foods <- rbindlist(lapply(files, function(i){
	food <- data.table(read_parquet(paste0(path, i,'_food.parquet')))[['t']]
	data <- data.table(read_parquet(paste0(path, i,'.parquet')))
	mint <- data[N > 0, min(Frame)] /2
	r <- range(food) - mint
	data.table(mint = r[1], maxt = r[2],
		   rho = round(as.numeric(strsplit(unlist(regmatches(i, 
		   						  gregexpr('rho_\\d{1}.\\d{1,3}',
		   						  	 i))), split = '_')[[1]][2]), 2))
	
}), idcol = TRUE)[is.finite(mint),lapply(.SD, mean),
		  .SDcols = c('mint', 'maxt'), by = 'rho']


data_peak <- rbindlist(data, idcol = TRUE)
data_peak[['rho']] <- vapply(strsplit(unlist(regmatches(data_peak[['.id']],
						       gregexpr('rho_\\d{1}.\\d{1,3}',data_peak[['.id']]))),
				     split = '_'), function(i) round(as.numeric(i[2]), 2), numeric(1))


# comparison of population dynamics
a <- ggplot(data = data_peak[Frame < 175*120], aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.2, linewidth = 0.8)+
	geom_path(data = det_N[Frame < 175*120], aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg[Frame < 175*120], color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg[Frame < 175*120], aes(color = 'Experiments'), linewidth = 2, alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('N (number of ants in arena)')+
	scale_color_manual('',values = c('grey30', 'grey80'), 
			   labels = c('Experimental', 'Model'))+
	geom_vline(data = foods[, .(mint = mint * 2, maxt = maxt*2,
				    rho = rho)], linetype = 2, linewidth = 1,
		   aes(xintercept = mint))+
	geom_vline(data = foods[, .(mint = mint * 2, maxt = maxt*2,
				    rho = rho)], linetype = 2, linewidth = 1,
		   aes(xintercept = maxt))+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	# geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	# geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	theme(legend.position = c(0.9, 0.9),
	      legend.background = element_rect(fill = 'white', color = 'black'),
	      legend.title = element_blank(),plot.title = element_text(size = 22))+
	geom_text(data = data.frame(x = c(120*3, 63*120, 120*3, 63*120,
					  120*3, 55*120, 120*3, 55*120),
				    y = 43, label = c('TP1', 'TP2'),
					  rho = c(0, 0, 0.2, 0.2,
					  	0.8, 0.8, 1, 1)),
		  aes(x, y, label = label), size = 6)+
	facet_wrap(~ rho, labeller = as_labeller(function(x){
		paste0('Proportion of scouts = ', x)
	}))
