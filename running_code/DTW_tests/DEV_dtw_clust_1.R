source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/det.RData')

load("/home/polfer/research/gits/AutomatAnts/results/ref_seq_det.RData")
library(dtw)

## argument parser ...
# args <- commandArgs()
# ........
args <- argument_parser('')

path <- '~/research/gits/AutomatAnts/results/with_recruitment/food_conditions/det/'
files <- list.files(path)
files <- files[!grepl('food', files) & !grepl('positions', files)]


det_s <- lapply(seq_along(files), function(i){
	new('Simulation', data = read.csv(paste0(path, files[i])),
	    type = gsub('.csv', '', files[i]))
})

if(args[['plots']]){
	for(i in seq_along(det_s)){
		
		png(paste0(path, '../plots/det/', det_s[[i]]@type, '.png'), 1920, 1080, res = 120)
		print(ggplot(data = det_s[[i]]@data, aes(Frame, N)) + geom_path()+
		      	scale_x_continuous('Time (min)', breaks = seq(0, 21600, 120*15),
		      			   labels = seq(0, 180, 15))+
		      	scale_y_continuous('Population', breaks = seq(0, 100, 5)))
		dev.off()
	}
}


### DYNAMIC TIME WARPING CLUSTERIZATION ###
det_N <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- i@data[, .N, by = 'Frame']
	y <- merge.data.table(data.table(Frame = seq(1, 21600)), x, all = TRUE, by = 'Frame')
	y[is.na(N), 'N'] <- 0
	y[1:21600, .(N = mean(N)), by = .(Frame = ((Frame -1) %/% 60) * 60)]
}), idcol = TRUE)

det_I <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- i@data[Crossings != 0, .(I = .N), by = 'Frame']
	y <- merge.data.table(data.table(Frame = seq(1, 21600)), x, all = TRUE, by = 'Frame')
	y[is.na(I), 'I'] <- 0
	y[1:21600, .(I = mean(I)), by = .(Frame = ((Frame -1) %/% 60) * 60)]
}), idcol = TRUE)

detS_N <- rbindlist(lapply(det_s, function(i){
	y <- merge.data.table(data.table(Frame = seq(1, 21600)), i@data[, c('Frame', 'N')],
			      all = TRUE, by = 'Frame')
	y[is.na(N), 'N'] <- 0
	y[1:21600, .(N = mean(N)), by = .(Frame = ((Frame ) %/% 60) * 60)]
}), idcol = TRUE)

detS_I <- rbindlist(lapply(det_s, function(i){
	y <- merge.data.table(data.table(Frame = seq(1, 21600)), i@data[, c('Frame', 'I')],
			      all = TRUE, by = 'Frame')
	y[is.na(I), 'I'] <- 0
	y[1:21600, .(I = mean(I)), by = .(Frame = ((Frame ) %/% 60) * 60)]
}), idcol = TRUE)


dtw_results <- rbindlist(lapply(seq_len(max(detS_N[['.id']])), function(i){
	rbindlist(lapply(seq_len(max(det_N[['.id']])), function(j){
		data.table(nrm_d = dtw(detS_N[.id == i, N], det_N[.id == j, N], 
				       step.pattern = asymmetric,
				       distance.only = TRUE)[['normalizedDistance']],
			   .id = i, .id_exp = j)
	}))
}))

dtw_results_I <- rbindlist(lapply(seq_len(max(detS_I[['.id']])), function(i){
	rbindlist(lapply(seq_len(max(det_I[['.id']])), function(j){
		data.table(nrm_d = dtw(detS_I[.id == i, I], det_I[.id == j, I], 
				       step.pattern = asymmetric,
				       distance.only = TRUE)[['normalizedDistance']],
			   .id = i, .id_exp = j)
	}))
}))
dtw_classification <- dtw_results[, .(nrm_d = mean(nrm_d)), by = '.id']
dtw_class_I <- dtw_results_I[, .(nrm_d = mean(nrm_d)), by = '.id']
set(dtw_classification, j = 'k', 
    value = ifelse(dtw_classification[['nrm_d']] > 4.5, 0, 1))

K_IDX <- dtw_classification[k == 1, .id]

### FORAGING EFFICIENCY ###
food_files <- list.files(path)[grepl('food', list.files(path))]

detExpEff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)
	x[, .(mint = min(t/2), maxt = max(t/2), meant = mean(t/2), medt = median(t/2),
	      intt = (max(t/2) - min(t/2)))]
	# x[, .(mint = 1/ (min(t/2)), maxt = 1/ (max(t/2)), 
	#       meant = 1/ (mean(t/2)), medt = 1/(median(t/2)), intt = 1/(max(t/2) - min(t/2)))]
}))

detEff <- rbindlist(lapply(K_IDX, function(i){
	rfood <- setDT(read.csv(paste0(path, food_files[i])))
	rfood[, .(mint = min(t), maxt = max(t), meant = mean(t), medt = median(t),
		  intt = (max(t) - min(t)))]
	# rfood[, .(mint = 1/min(t), maxt = 1/max(t), meant = 1/mean(t), medt = 1/median(t))]
}))



### CROSSED ANALYSIS -- DTW ###
dddtw_x <- matrix(0, ncol = max(detS_N[['.id']]), nrow = max(detS_N[['.id']]))
for(i in seq_len(ncol(dddtw_x))){
	dddtw_x[, i] <- vapply(seq_len(nrow(dddtw_x)), function(j){
		dtw(detS_N[.id == i, N], detS_N[.id == j, N], 
	    step.pattern = asymmetric,
	    distance.only = TRUE)[['normalizedDistance']]}, numeric(1))
}
xhclust <- hclust(as.dist(dddtw_x))
xhclust <- hclust(as.dist(dddtw_x**2))
t <- cutree(xhclust, k = 3)
t <- cutree(xhclust, k = 4)
avghclust <- hclust(as.dist(dddtw_x), method = 'centroid')



dddtw_i <- matrix(0, ncol = max(detS_I[['.id']]), nrow = max(detS_I[['.id']]))
for(i in seq_len(ncol(dddtw_i))){
	dddtw_i[, i] <- vapply(seq_len(nrow(dddtw_i)), function(j){
		dtw(detS_I[.id == i, I], detS_I[.id == j, I], 
		    step.pattern = asymmetric,
		    distance.only = TRUE)[['normalizedDistance']]}, numeric(1))
}

whclust <- hclust(as.dist(dddtw_x * dddtw_i))
plot(xhclust)
plot(whclust)
t <- cutree(whclust, h = 10) 
# dtw_x <- rbindlist(lapply(1:(max(detS_N[['.id']])-1), function(i){
# 	rbindlist(lapply(2:(max(detS_N[['.id']])), function(j){
# 		data.table(nrm_d = dtw(detS_N[.id == i, N], detS_N[.id == j, N], 
# 				       step.pattern = asymmetric,
# 				       distance.only = TRUE)[['normalizedDistance']],
# 			   .id = i, .id_exp = j)
# 	}))
# }))
dtw_x <- dtw_x[nrm_d > 0]
dtw_y <- dcast(dtw_x, formula = .id ~ .id_exp, value.var = 'nrm_d', fill = 0)
dtw_y[['.id']] <- NULL
tsne <- bdm.bhtsne(dtw_y, 'dtw_det_sims', ppx = 100, is.distance = T, threads = 4)
tsne <- bdm.pakde(tsne, ppx = 5, threads = 4)
tsne <- bdm.wtt(tsne)
tsne <- bdm.merge.s2nr(as.matrix(dtw_y), tsne, k = 3, ret.merge = T, info = F)

ggplot(data = dtw_x, aes(.id, .id_exp, fill = nrm_d)) + geom_raster()+
	theme_void() + theme(aspect.ratio = 1) +
	scale_fill_viridis_c('Normalized DTW distance', direction = -1, 
			     option = 'magma', breaks = seq(0, 30, 5))



ids <- seq_along(t)[ t%in% c(1, 3)]
ggplot(data = detS_N[.id %in% ids], aes(Frame, N)) + facet_wrap(~factor(.id))+
	geom_path()


### EXTREME EVENTS ###
ids <- c(which(grepl('det_10', files)),
	 which(grepl('det_20', files)),
	 which(grepl('det_29', files)))
det10 <- detS_N[.id == ids[1]] # late-start
det20 <- detS_N[.id == ids[2]] # plateau
det29 <- detS_N[.id == ids[3]] # experiment
w <- c(rep(1000, nrow(det10)/2), rep(0.000001, nrow(det10)/2))
w <- 1:nrow(det20)

dtw_refs <- rbindlist(lapply(seq_len(max(detS_N[['.id']])), function(i){
	data.table(d_10 = dtw(detS_N[.id == i, N], det10[, N], 
			      step.pattern = asymmetric,
			      distance.only = TRUE,
			      weights = w)[['normalizedDistance']],
		   d_20 = dtw(detS_N[.id == i, N], det20[, N], 
		   	   step.pattern = asymmetric,
		   	   distance.only = TRUE, 
		   	   weights = w)[['normalizedDistance']],
		   d_29 = dtw(detS_N[.id == i, N], det29[, N], 
		   	   step.pattern = asymmetric,
		   	   distance.only = TRUE,
		   	   weights = w)[['normalizedDistance']],
		   .id = i)
}))


ggplot(data = dtw_refs, aes(d_10, d_20, color = d_29)) + geom_point(size = 5)







k <- tsclust(x, k = 3L, distance = 'dtw', type = 'fuzzy')
hk <- k@cluster

ggplot(data = as.data.frame(k@fcluster), aes(cluster_1, cluster_2, color = cluster_3))+
	geom_point(size = 5)


kp <- tsclust(x, k = 3L, distance = 'dtw', type = 'partitional', centroid = 'shape')
ggplot(data = as.data.frame(kp@fcluster), aes(cluster_1, cluster_2, color = cluster_3))+
	geom_point(size = 5)

kp2 <- tsclust(x, k = 4L, distance = 'dtw', type = 'partitional', centroid = 'shape')
ggplot(data = as.data.frame(kp@fcluster), aes(cluster_1, cluster_2, color = cluster_3))+
	geom_point(size = 5)

kf <- tsclust(x, k = 3L, type = 'fuzzy', shape = 'fcmdd')
plot(kf)

ndtw <- function(x, y, ...) {
	dtw::dtw(x, y, step.pattern = asymmetric,
		 distance.only = TRUE, ...)$normalizedDistance
}

if (!pr_DB$entry_exists("nDTW")){
	proxy::pr_DB$set_entry(FUN = ndtw, names=c("nDTW"),
			       loop = TRUE, type = "metric", distance = TRUE,
			       description = "Normalized asymmetric DTW")
}


kplan <- tsclust(x, k = 6L, distance = 'nDTW')

kplan_avg <- tsclust(x, k = 3L, distance = 'nDTW', centroid = 'shape', 
		  control = partitional_control(nrep = 100))
plot(kplan3)
kplan3



x <- lapply(det_s, function(i){
	y <- merge.data.table(data.table(Frame = seq(1, 21600)), i@data[, c('Frame', 'N')],
			      all = TRUE, by = 'Frame')
	y[is.na(N), 'N'] <- 0
	y[1:21600, .(N = mean(N)), by = .(Frame = ((Frame ) %/% 10) * 10)][['N']]
})
kx <- tsclust(x, type = 'partitional', k = 3L, distance = 'dtw_basic', preproc = zscore,
	centroid = 'shape')

