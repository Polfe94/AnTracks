#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(latex2exp)
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')

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

n_data <- function(path, label){
	
	files <- list.files(path)
	files <- files[!grepl('food', files) & !grepl('position', files)]
	
	rbindlist(lapply(files, function(i){
		x <- data.table(read.csv(paste0(path, i)))[, c('Frame', 'N')]
		x[['label']] <- label
		x
	}), idcol = TRUE)
}



global_Eff <- rbind(foodEff('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_10/',
		      label = 'bimodal_10'),
		    foodEff('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_90/',
		    	label = 'bimodal_90'),
		    foodEff('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_50/',
		    	label = 'bimodal_50'),
		    foodEff('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal/',
		    	label = 'bimodal'),
		    foodEff('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gaussian/',
		    	label = 'gaussian'))

# takes 1 min
sims_N <- rbind(n_data('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_10/',
			label = 'bimodal_10'),
		n_data('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_90/',
			label = 'bimodal_90'),
		n_data('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal_50/',
			label = 'bimodal_50'),
		n_data('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/bimodal/',
			label = 'bimodal'),
		n_data('/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/gaussian/',
			label = 'gaussian'))




## COMPUTAR EFICIENCIA EXPERIMENTAL
det_Eff <- rbindlist(lapply(det, function(i){rbindlist(i@food)}), idcol = TRUE)
det_Eff_mlt <- reshape2::melt(det_Eff[, .(eff_1 = 1/ min(t/2), 
					     eff_2 = 1/(max(t/2)-min(t/2)),
					     eff_4 = 1 / mean(diff(sort(t / 2))),
					     no_found = 1,
					     label = 'determinist'), by = .id], id.vars = 'label')

sto_Eff <- rbindlist(lapply(sto, function(i){rbindlist(i@food)}), idcol = TRUE)
sto_Eff_mlt <- reshape2::melt(sto_Eff[, .(eff_1 = 1/ min(t/2), 
					  eff_2 = 1/(max(t/2)-min(t/2)),
					  eff_4 = 1 / mean(diff(sort(t/2))),
					  no_found = 1,
					  label = 'stochastic'), by = .id], id.vars = 'label')

exp_Eff <- rbind(rbindlist(lapply(det, function(i){rbindlist(i@food)})),
		 rbindlist(lapply(sto, function(i){rbindlist(i@food)})))
exp_Eff[['.id']] <- unlist(lapply(1:20, function(i) rep(i, 12)))



exp_Eff_mlt <- reshape2::melt(exp_Eff[, .(eff_1 = 1/ min(t/2), 
					  eff_2 = 1/(max(t/2)-min(t/2)),
					  eff_3 = 1/(max(t/2)-min(t/2)),
					  eff_4 = 1 / mean(diff(sort(t/2))),
					  no_found = 1,
					  label = 'experimental'), by = .id], id.vars = c('.id', 'label'))
ggplot(data = exp_Eff_mlt[exp_Eff_mlt$variable != 'eff_3',], aes(factor(label), value)) +
	geom_boxplot() + facet_wrap(~ factor(variable, labels = c('Discovery efficiency',
								  'Retrieval efficiency',
								  'Mean consecutive efficiency',
								  'Proportion of found food')),
				    scales = 'free')+
	scale_x_discrete('',labels = TeX)+ ylab('')+
	theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
	      axis.text = element_text(size = 12))


exp_Eff_xCapita <- rbind(rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- rbindlist(i@food)
	mint <- min(x[['t']])
	maxt <- max(x[['t']])
	x[['N1']] <- i@data[Frame == mint, .N]
	x[['N2']] <- mean(i@data[Frame > mint & Frame <= maxt, .N, by = 'Frame'][['N']])
	x
})),rbindlist(lapply(sto, function(i){
	setDT(i@data)
	x <- rbindlist(i@food)
	mint <- min(x[['t']])
	maxt <- max(x[['t']])
	x[['N1']] <- i@data[Frame == mint, .N]
	x[['N2']] <- mean(i@data[Frame > mint & Frame <= maxt, .N, by = 'Frame'][['N']])
	x
})))
exp_Eff_xCapita[['.id']] <- unlist(lapply(1:20, function(i) rep(i, 12)))

exp_EffxCapita_mlt <- reshape2::melt(exp_Eff_xCapita[, .(eff_1 = 1/ (min(t/2)*N1), 
					  eff_2 = 1/((max(t/2)-min(t/2))*N2),
					  eff_3 = 1/((max(t/2)-min(t/2))*N2),
					  eff_4 = 1 / mean(diff(sort(t/2))),
					  no_found = 1,
					  label = 'experimental'), by = .id], id.vars = c('.id', 'label'))

ggplot(data = exp_EffxCapita_mlt[exp_EffxCapita_mlt$variable != 'eff_3',], aes(factor(label), value)) +
	geom_boxplot() + facet_wrap(~ factor(variable, labels = c('Discovery efficiency',
								  'Retrieval efficiency',
								  'Mean consecutive efficiency',
								  'Proportion of found food')),
				    scales = 'free')+
	scale_x_discrete('',labels = TeX)+ ylab('')+
	theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
	      axis.text = element_text(size = 12))
## COMPUTAR EFICIENCIES PER CAPITA (TAN DE LES SIMUS COM EXPERIMENTALS)


global_Eff_mlt <- reshape2::melt(global_Eff[, .(eff_1 = 1/mint, 
						eff_2 = 1/(maxt-mint),
						eff_3 = 1/(ifelse(na > 0, 10800, maxt)-mint),
						eff_4 = 1 / inter_mean,
						no_found = round((12-na) / 12,2),
						label = label)], id.vars = 'label')

ggplot(data = global_Eff_mlt, aes(factor(label, levels = c('bimodal', 'gaussian','bimodal_50',
							   'bimodal_90', 'bimodal_10'),
					 # labels = c('$\\beta(0.5, 0.5)$',
					 # 	   '$\\beta(2, 2)$',
					 # 	   '$0.5 \\times \\beta(2, 6)+ 0.5 \\times \\beta(6, 2)$',
					 # 	   '$0.9 \\times\\beta(2, 6)+ 0.1 \\times \\beta(6, 2)$',
					 # 	   '$0.1 \\times\\beta(2, 6)+ 0.9 \\times \\beta(6, 2)$')
					 ),
				  value)) +
	geom_boxplot() + facet_wrap(~ factor(variable, labels = c('Discovery efficiency',
								  'Retrieval efficiency',
								  'Corrected retrieval efficiency',
								  'Mean consecutive efficiency',
								  'Proportion of found food')),
				    scales = 'free')+
	scale_x_discrete('',labels = TeX)+ ylab('')+
	theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
	      axis.text = element_text(size = 12))


##### EFFICIENCY INCLUDING EXPERIMENTAL VALUES
global_Eff[['Frame_min']] <- round(global_Eff[['mint']]*2)
global_Eff[['Frame_max']] <- round(global_Eff[['maxt']]*2)

global_Eff[['N1']] <- vapply(1:nrow(global_Eff), function(i){
	mint = global_Eff[i, Frame_min]
	maxt = global_Eff[i, Frame_max]
	mean(sims_N[.id == global_Eff[i, .id] &
	       	label == global_Eff[i, label] &
	       	Frame == mint, N])
}, numeric(1))

global_Eff[['N2']] <- vapply(1:nrow(global_Eff), function(i){
	mint = global_Eff[i, Frame_min]
	maxt = global_Eff[i, Frame_max]
	mean(sims_N[.id == global_Eff[i, .id] &
		    	label == global_Eff[i, label] &
		    	Frame > mint & Frame <= maxt, N])
}, numeric(1))

# sapply(1:nrow(global_Eff), function(i){
# 	mint = global_Eff[i, Frame_min]
# 	maxt = global_Eff[i, Frame_max]
# 	c(N1 = sims_N[.id == global_Eff[i, .id] &
# 	       	label == global_Eff[i, label] &
# 	       	Frame == mint, N],
# 	N2 = mean(sims_N[.id == global_Eff[i, .id] &
# 		     	label == global_Eff[i, label] &
# 		     	Frame > mint & Frame <= maxt, N]))
# })
global_EffxCapita_mlt <- reshape2::melt(global_Eff[, .(eff_1 = 1/(mint*N1), 
						eff_2 = 1/((maxt-mint)*N2),
						eff_3 = 1/((ifelse(na > 0, 10800, maxt)-mint)*N2),
						eff_4 = 1 / inter_mean,
						no_found = round((12-na) / 12,2),
						label = label)], id.vars = 'label')



###########

all_Eff <- rbind(global_Eff_mlt, exp_Eff_mlt[, c('label', 'variable', 'value')])
comparison_efficiencies <- ggplot(data = all_Eff, aes(factor(label, levels = c('experimental', 'bimodal', 'gaussian','bimodal_50',
							   'bimodal_90', 'bimodal_10'),
),
value)) +
	geom_boxplot() + facet_wrap(~ factor(variable, labels = c('Discovery efficiency',
								  'Retrieval efficiency',
								  'Corrected retrieval efficiency',
								  'Mean consecutive efficiency',
								  'Proportion of found food')),
				    scales = 'free')+
	scale_x_discrete('',labels = TeX)+ ylab('')+
	theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
	      axis.text = element_text(size = 12))

png('/home/polfer/research/gits/AnTracks/plots/efficiency_sims_with_recruitment.png', 1920, 1080, res = 100)
comparison_efficiencies
dev.off()

vars2test <- c('bimodal', 'gaussian','bimodal_50', 'bimodal_90', 'bimodal_10')

## DISCOVERY EFFICIENCY
for(i in vars2test){
	print(paste0('Experimental vs. ', i, ' Wilcoxon Test [DISCOVERY EFFICIENCY]'))
	print(wilcox.test(x = exp_Eff_mlt[exp_Eff_mlt$label == 'experimental' & exp_Eff_mlt$variable == 'eff_1', 'value'],
		     y = global_Eff_mlt[global_Eff_mlt$label == i & global_Eff_mlt$variable == 'eff_1', 'value'],
			  digits.rank = 7)$p.value)
}

## RETRIEVAL EFFICIENCY
for(i in vars2test){
	print(paste0('Experimental vs. ', i, ' Wilcoxon Test [RETRIEVAL EFFICIENCY]'))
	print(wilcox.test(x = exp_Eff_mlt[exp_Eff_mlt$label == 'experimental' & exp_Eff_mlt$variable == 'eff_2', 'value'],
			  y = global_Eff_mlt[global_Eff_mlt$label == i & global_Eff_mlt$variable == 'eff_2' &
			  		   	is.finite(global_Eff_mlt$value),
			  		   'value'],
			  digits.rank = 7)$p.value)
}

##### EFFICIENCY PER CAPITA INCLUDING EXPERIMENTAL VALUES

all_Eff_xCapita <- rbind(global_EffxCapita_mlt, exp_EffxCapita_mlt[, c('label', 'variable', 'value')])
filtered_ <- all_Eff_xCapita
filtered_$value[filtered_$variable == 'eff_1' & filtered_$value > 0.005] <- NA
filtered_$value[filtered_$variable == 'eff_2' & filtered_$value > 0.0005] <- NA
filtered_$value[filtered_$variable == 'eff_3' & filtered_$value > 3e-4] <- NA
comparison_efficiencies_xCapita <- ggplot(data = filtered_[filtered_$variable != 'no_found', ],
       aes(factor(label, levels = c('experimental', 'bimodal', 'gaussian','bimodal_50',
						    'bimodal_90', 'bimodal_10'),
),
value)) +
	geom_boxplot() + facet_wrap(~ factor(variable, 
					     labels = c('Discovery efficiency (per capita)',
					     	   'Retrieval efficiency (per capita)',
					     	   'Corrected retrieval efficiency (per capita)',
					     	   'Mean consecutive efficiency (per capita)'#,
					     	   #'Proportion of found food')
					     )),
				    scales = 'free_y')+
	scale_x_discrete('',labels = TeX)+ ylab('')+
	theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
	      axis.text = element_text(size = 12))

png('/home/polfer/research/gits/AnTracks/plots/efficiency_sims_with_recruitment_xCapita.png', 1920, 1080, res = 100)
comparison_efficiencies_xCapita
dev.off()


## DISCOVERY EFFICIENCY
for(i in vars2test){
	print(paste0('Experimental vs. ', i, ' Wilcoxon Test [DISCOVERY EFFICIENCY]'))
	print(wilcox.test(x = all_Eff_xCapita[all_Eff_xCapita$label == 'experimental' & all_Eff_xCapita$variable == 'eff_1', 'value'],
			  y = all_Eff_xCapita[all_Eff_xCapita$label == i & all_Eff_xCapita$variable == 'eff_1', 'value'],
			  digits.rank = 7)$p.value)
}

## RETRIEVAL EFFICIENCY
for(i in vars2test){
	print(paste0('Experimental vs. ', i, ' Wilcoxon Test [RETRIEVAL EFFICIENCY]'))
	print(wilcox.test(x = all_Eff_xCapita[all_Eff_xCapita$label == 'experimental' & all_Eff_xCapita$variable == 'eff_2', 'value'],
			  y = all_Eff_xCapita[all_Eff_xCapita$label == i & all_Eff_xCapita$variable == 'eff_2' &
			  		   	is.finite(all_Eff_xCapita$value),
			  		   'value'],
			  digits.rank = 7)$p.value)
}
