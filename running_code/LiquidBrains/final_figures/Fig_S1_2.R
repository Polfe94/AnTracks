source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')
# load('~/research/gits/AnTracks/results/mov_rho_R.RData')


maxt <- 20*120

library(trajr)
library(arrow)
library(egg)
library(latex2exp)

DET_scouts <- rbindlist(lapply(c(5, 10, 15, 20, 25) * 120,
			       function(t){
			       	
			       	rbindlist(lapply(det, function(p){
			       		data <- setDT(p@data)[Frame < t + min(Frame)]
			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
			       		set(data, j = 'd', value = as.numeric(dmatrix))
			       		data[, .(maxd = max(d), maxt = t), by = 'N_ind']
			       	}))
			       }))
DET_scouts[['tag']] <- 'DET'

NFD_scouts <- rbindlist(lapply(c(5, 10, 15, 20, 25) * 120,
			       function(t){
			       	
			       	rbindlist(lapply(nf[-c(1, 2)], function(p){
			       		data <- setDT(p@data)[Frame < maxt + min(Frame)]
			       		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
			       		set(data, j = 'd', value = as.numeric(dmatrix))
			       		data[, .(maxd = max(d), maxt = t), by = 'N_ind']
			       	}))
			       }))
NFD_scouts[['tag']] <- 'NFD'


png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_S1_2.png', 6000, 3000, res = 350)
ggplot(data = rbind(DET_scouts, NFD_scouts), aes(maxd, fill = tag)) + 
	geom_density(aes(fill = tag, y = after_stat(density)),
		     alpha = 0.9,  color = 'black')+
	geom_histogram(alpha = 0.5, color = 'black', aes(y = after_stat(density)),
		       bins = 25) +
	facet_wrap(~ maxt, scales = 'free',
		   labeller = as_labeller(function(i){
		   	paste0(
		   		as.numeric(i)/120,
		   		' min')
		   }))+
	ylab('Density') + 
	scale_x_continuous('Radius (mm)', breaks = seq(0, 1500, 250))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.position = c(0.8, 0.25))+
	geom_vline(xintercept = 650, linetype = 2)
dev.off()