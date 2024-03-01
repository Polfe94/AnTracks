#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
# library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')


path <- '/home/polfer/research/gits/AutomatAnts/results/2024/g_eff/'
conditions <- list.files(path)

result <- vector('list', length(conditions))
names(result) <- conditions

t0 <- Sys.time()
for(i in conditions){
	f <- list.files(paste0(path, i))
	files <- f[grepl('food', f)]
	l <- lapply(files, function(x){
		g <- as.numeric(gsub('g_', '', unlist(regmatches(x, gregexpr('g_\\d{1,2}\\.\\d{1,2}', x)))))
		data <- data.table(read_parquet(paste0(path, i, '/', x)))
		data.frame(tp1 = min(data[['t']]), tp2 = max(data[['t']]), g = g)
	})
	result[[i]] <- rbindlist(l, fill = TRUE)
}
Sys.time()-t0

ggplot(data = result[[3]][is.finite(tp1), .(tp1 = 1/mean(tp1)), by = 'g'], aes(g, tp1)) + geom_point()
# ggplot(data = result[[3]][is.finite(tp2), .(tp2 = mean(tp2)), by = 'g'], aes(g, tp2)) + geom_point()
ggplot(data = result[[3]][is.finite(tp1), .(exploit = 1/mean(tp2-tp1)), by = 'g'], aes(g, exploit)) + geom_point()

ggplot(data = rbindlist(result, idcol = TRUE)[!is.finite(tp1), tp1 := 0]
       [is.finite(tp1), .(tp1 = 1/mean(tp1)), by = c('.id', 'g')]
       [!is.finite(tp1), tp1 := 0], 
       aes(g, tp1, color = .id)) +
	geom_point(size = 5, shape = 21) + geom_line(linewidth = 1)+
	scale_color_viridis_d()

ggplot(data = rbindlist(result, idcol = TRUE)[!is.finite(tp2), tp2 := 0]
       [!is.finite(tp1), tp1 := 0]
       [is.finite(tp2), .(tp2 = 1/mean(tp2-tp1)), by = c('.id', 'g')]
       [!is.finite(tp2), tp2 := 0], 
       aes(g, tp2, color = .id)) +
	geom_point(size = 5, shape = 21) + geom_line(linewidth = 1)+
	scale_color_viridis_d()





#### WITH MEDIANS
ggarrange(
ggplot(data = rbindlist(result, idcol = TRUE)[!is.finite(tp1), tp1 := 0]
       [is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('.id', 'g')]
       [!is.finite(tp1), tp1 := 0], 
       aes(g, tp1, color = .id)) +
	geom_point(size = 5, shape = 21) + geom_line(linewidth = 1)+
	scale_color_viridis_d()+
	ylab('TP1')
,
ggplot(data = rbindlist(result, idcol = TRUE)[!is.finite(tp2), tp2 := 0]
       [!is.finite(tp1), tp1 := 0]
       [is.finite(tp2), .(tp2 = 1/median(tp2-tp1)), by = c('.id', 'g')]
       [!is.finite(tp2), tp2 := 0], 
       aes(g, tp2, color = .id)) +
	geom_point(size = 5, shape = 21) + geom_line(linewidth = 1)+
	scale_color_viridis_d()+
	ylab('TP2-TP1'),
ncol = 1)




scenario_collective_efficiency <- vector('list', length(conditions))

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/g_eff/'

lapply(seq_along(scenario_collective_efficiency), function(f){
	do.call('gc', args = list(verbose = FALSE))
	files <- list.files(paste0(path, conditions[f]))
	x <- eff(path, files[i])
	x[['g']] <- as.numeric(gsub('g_', '', unlist(regmatches(f[i], gregexpr('g_\\d{1}\\.\\d{1,2}', f[i])))))
	x
})
t0 <- Sys.time()
coll_eff <- lapply(seq_along(f), function(i){
	do.call('gc', args = list(verbose = FALSE))
	x <- eff(path, files[i])
	x[['g']] <- as.numeric(gsub('g_', '', unlist(regmatches(f[i], gregexpr('g_\\d{1}\\.\\d{1,2}', f[i])))))
	x
})
Sys.time()-t0
