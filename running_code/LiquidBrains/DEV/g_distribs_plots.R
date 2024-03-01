# required libraries
library(data.table)
library(ggplot2)

# generic paths
data_path <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameter_space/'
plot_path <- '/home/polfer/research/gits/AutomatAnts/plots/'

# generic functions
easyPlot <- function(data_path, plot_path){
	files <- list.files(data_path)
	files <- files[!grepl('food', files) & !grepl('positions', files)]
	

	
	for(i in seq_along(files)){
		dt <- data.table(read.csv(paste0(data_path, files[i])))
		png(paste0(plot_path, gsub('.csv', '', files[i]), '.png'),
		    1920, 1080, res = 120)
		print(ggplot(data = dt, aes(Frame, N)) + geom_path())
		dev.off()
	}
}

foodEff <- function(data_path, label = NULL){
	files <- list.files(data_path)
	files <- files[grepl('food', files)]
	
	result <- rbindlist(lapply(seq_along(files), function(i){
		x <- read.csv(paste0(data_path, files[i]))[, 2]
		data.table(mint = min(x), maxt = max(x))
	}), idcol = TRUE)
	
	if(!is.null(label)){
		set(result, j = 'label', value = label)
	}
	result
}

#### GAUSSIAN SIMULATIONS ####
easyPlot(paste0(data_path, 'gaussian/'), paste0(plot_path, 'gaussian_gains/'))

#### BIMODAL SIMULATIONS ####
easyPlot(paste0(data_path, 'bimodal/'), paste0(plot_path, 'bimodal_gains/'))

#### BIMODAL 50% SIMULATIONS ####
easyPlot(paste0(data_path, 'bimodal_50/'), paste0(plot_path, 'bimodal_50_gains/'))

#### BIMODAL 10% SIMULATIONS (LEFT-SKEWED) ####
easyPlot(paste0(data_path, 'bimodal_10/'), paste0(plot_path, 'bimodal_10_gains/'))

#### BIMODAL 90% SIMULATIONS (RIGHT-SKEWED) ####
easyPlot(paste0(data_path, 'bimodal_90/'), paste0(plot_path, 'bimodal_90_gains/'))


fef_gauss <- foodEff(paste0(data_path, 'gaussian/'), label = 'gaussian')
fef_bmd <- foodEff(paste0(data_path, 'bimodal/'), label = 'bimodal')
fef_bmd50 <- foodEff(paste0(data_path, 'bimodal_50/'), label = 'bimodal_50')
fef_bmd10 <- foodEff(paste0(data_path, 'bimodal_10/'), label = 'bimodal_10')
fef_bmd90 <- foodEff(paste0(data_path, 'bimodal_90/'), label = 'bimodal_90')

fef <- rbindlist(list(fef_gauss, fef_bmd, fef_bmd50, fef_bmd10, fef_bmd90))
# ggplot(data = reshape2::melt(fef, id.vars = c('.id', 'label')), 
#        aes(label, value, fill = variable))+
# 	geom_boxplot() + scale_fill_viridis_d()

ggplot(data = fef, aes(label, 1/mint)) + 
	geom_boxplot(outlier.shape = NA) + geom_jitter()+
	xlab('') + ylab('Food discovery efficiency')
ggplot(data = fef, aes(label, 1/(maxt-mint))) + 
	geom_boxplot(outlier.shape = NA) + geom_jitter()+
	xlab('')+ylab('Food retrieval efficiency')


ggplot(data = reshape2::melt(fef, id.vars = c('.id', 'label')), 
       aes(label, 1/value))+
	geom_boxplot()+facet_wrap(~ variable, scales = 'free') 
