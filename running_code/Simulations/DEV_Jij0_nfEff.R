#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
# library(dtw)
# library(dtwclust)
library(latex2exp)
# source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')
load('~/research/gits/AnTracks/data/nf.RData')


foodEff <- function(data_path, label = NULL){
	files <- list.files(data_path)
	files <- files[grepl('food', files)]
	
	result <- rbindlist(lapply(seq_along(files), function(i){
		x <- read.csv(paste0(data_path, files[i]))[['t']]
		data.table(mint = min(x, na.rm = TRUE), 
			   maxt = max(x, na.rm = TRUE),
			   na = sum(is.na(x)),
			   inter_mean = mean(diff(sort(x))))
	}), idcol = TRUE)
	
	if(!is.null(label)){
		set(result, j = 'label', value = label)
	}
	result
}


f <- get_foodPatches('20180720M')
nf <- lapply(nf, function(i){
	i@food <- f[1:2]
	i <- food_detection(i)
	i
})

exp_Eff <- rbind(rbindlist(lapply(det, function(i){rbindlist(i@food)})),
		 rbindlist(lapply(sto, function(i){rbindlist(i@food)})))
exp_Eff[['.id']] <- unlist(lapply(1:20, function(i) rep(i, 12)))

nf_Eff <- rbindlist(lapply(nf, function(i) rbindlist(i@food)), idcol = TRUE)
nf_Eff_mlt <- reshape2::melt(nf_Eff[, .(eff_1 = 1/(min(t/2)),
					eff_2 = 1/(max(t/2)-min(t/2)),
					label = 'no_food'), by = .id], id.vars = c('.id', 'label'))


exp_Eff_mlt <- reshape2::melt(exp_Eff[, .(eff_1 = 1/ min(t/2), 
					  eff_2 = 1/(max(t/2)-min(t/2)),
					  label = 'experimental'), by = .id], id.vars = c('.id', 'label'))

ggplot(data = rbind(exp_Eff_mlt, nf_Eff_mlt), aes(factor(label), value)) +
	geom_boxplot() + facet_wrap(~ factor(variable, labels = c('Exploration efficiency',
								  'Exploitation efficiency')),
				    scales = 'free')+
	scale_x_discrete('')+ ylab('')+
	theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
	      axis.text = element_text(size = 12))

path_Jij0 <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/homogeneous_Jij/Jij_0/uniform/'
Jij0_Eff <- foodEff('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/homogeneous_Jij/Jij_0/uniform/', 'Jij = 0')

Jij0_mlt <- reshape2::melt(Jij0_Eff[, .(eff_1 = 1/mint, eff_2 = 1/(maxt-mint), label = label, .id = .id)],
			   id.vars = c('.id', 'label'))

ggplot(data = rbind(exp_Eff_mlt, nf_Eff_mlt, Jij0_mlt), aes(factor(label), value)) +
	geom_boxplot() + facet_wrap(~ factor(variable, labels = c('Exploration efficiency',
								  'Exploitation efficiency')),
				    scales = 'free')+
	scale_x_discrete('')+ ylab('')+
	theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
	      axis.text = element_text(size = 12))+
	coord_cartesian(ylim = c(0, 0.005))
