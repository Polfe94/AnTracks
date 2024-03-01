library(jsonlite)
source('~/research/gits/AnTracks/src/Experiment.R')


path2read1 <- "/home/polfer/research/gits/AnTracks/results/OLD_exps/20180717M/"
path2read2 <- "/home/polfer/research/gits/AnTracks/results/OLD_exps/20180717T/"
path2read3 <- "/home/polfer/research/gits/AnTracks/results/OLD_exps/20180717Mb/"

unveil_json <- function(json){
	require(data.table)
	rbindlist(lapply(seq_along(json), function(i){
		rbindlist(lapply(seq_along(json[[i]]), function(ii){
			a <- unlist(json[[i]][[ii]])
			data.frame(id = i, Frame = a[1], Xmm = a[2], Ymm = a[3], int = a[4])
		}))
	}))[order(Frame)]
}


tracks_1 <- read_json(paste0(path2read1, 'iTracks.json'))
tracks_2 <- read_json(paste0(path2read2, 'iTracks.json'))
tracks_3 <- read_json(paste0(path2read3, 'iTracks.json'))
data_20180717T <- unveil_json(tracks_1)
data_20180717M <- unveil_json(tracks_2)
data_20180717Mb <- unveil_json(tracks_3)

new_data <- new('Experiment', data = data_20180717T[Ymm > 1000], date = '20180717T')
new_data <- new('Experiment', data = data_20180717M[Ymm > 1000], date = '20180717M')
new_data <- new('Experiment', data = data_20180717Mb[Ymm > 1000], date = '20180717Mb')
new_data <- compute_nodes(new_data)
new_data <- food_detection(new_data)
geom_foodpatches(add = draw_hexagons(add = ggplot(data = new_data@data, aes(Xmm, Ymm))+geom_point()))
