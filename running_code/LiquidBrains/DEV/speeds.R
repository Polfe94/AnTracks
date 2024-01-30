spd <- det[[1]]@data[, .(dx = diff(Xmm), dy = diff(Ymm), dt = diff(Time_sec), t = Time_sec[-1]), by = 'N_ind']
spd[, .(v = sqrt(dx**2 + dy**2)/dt, t = t), by = 'N_ind']

fltrd_spd <- spd[, .(v = sqrt(dx**2 + dy**2)/dt, t = t), by = 'N_ind'][v <= quantile(v, p = 0.99)]

interval_width <- 1
time_width <- 60
fltrd_spd[, bin := round(v/interval_width)]
fltrd_spd[, time_bin := round(t/60)]

fltrd_spd[, .(v = .N), by = c('bin', 'time_bin')]


ggplot(data = fltrd_spd[, .(v = .N), by = c('bin', 'time_bin')][order(time_bin)], aes(x = bin/50, y = v)) + 
	geom_col(position = 'stack', aes(fill = time_bin)) +
	scico::scale_fill_scico('Time (min)', palette = 'vik', midpoint = 90)+
	scale_x_continuous('Speed (nodes / s)') +
	scale_y_continuous('Frequency', seq(0, 12000, length.out = 6))


sim_data <- sim_data <- data.table(read_parquet(paste0(path, files[1], '_data.parquet')))
sim_spd <- sim_data[, .(id = parse_ids(id_out), node = parse_nodes(pos)), by = 'T'][order(id, `T`)]
nds <- do.call('rbind', strsplit(gsub('[()]', '', sim_spd[['node']]), ', '))
set(sim_spd, j = 'x',value =  as.numeric(nds[, 1]))
set(sim_spd, j = 'y',value =  as.numeric(nds[, 2]))
sim_dts <- sim_spd[, .(dx = diff(x), dy = diff(y), dt = diff(`T`), t = `T`[-1]), by = 'id']
sim_dts[, group := c(0, cumsum(dx != 0 | dy != 0)[-length(dx)])]
# sim_dts[, group := (1+rleid(dx == 0 & dy == 0))%/% 2, by = 'id']
sim_dts <- sim_dts[, .(dx = sum(dx), dy = sum(dy), dt = sum(dt), t = tail(t, 1)), by = c('id', 'group')]
fltrd_sim_spd <- sim_dts[, .(v = sqrt(dx**2 + dy**2)/dt, t = t), by = c('id')][v <= quantile(v, p = 0.99)]


interval_width <- 0.25
time_width <- 60
fltrd_sim_spd[, bin := round(v/interval_width)]
fltrd_sim_spd[, time_bin := round(t/60)]



ggplot(data = fltrd_sim_spd[, .(v = .N), by = c('bin', 'time_bin')][order(time_bin)], aes(x = bin, y = v)) + 
	geom_col(position = 'stack', aes(fill = time_bin)) +
	scico::scale_fill_scico('Time (min)', palette = 'vik', midpoint = 90)+
	scale_x_continuous('Speed (nodes / s)') +
	coord_cartesian(xlim = c(0, 25))



mltd <- melt(fltrd_sim_spd, id.vars = 'id')

ggplot(data = mltd[mltd$variable == 'v', c('id', 'value')]) + geom_boxplot(aes(factor(id), value))+
	coord_cartesian(ylim = c(0, 5))









##### ALL SIMULATIONS #####
spd_det_p1 <- rbindlist(lapply(det, function(i){
	tp1 <- min(rbindlist(i@food)[['t']])
	spd <- i@data[Frame <= tp1, .(dx = diff(Xmm), dy = diff(Ymm), dt = diff(Time_sec), t = Time_sec[-1]), by = 'N_ind']
	spd[, .(v = sqrt(dx**2 + dy**2)/dt, t = t), by = 'N_ind'][v <= quantile(v, p = 0.99)]
}))

spd_det_p2 <- rbindlist(lapply(det, function(i){
	tp1 <- min(rbindlist(i@food)[['t']])
	tp2 <- max(rbindlist(i@food)[['t']])
	spd <- i@data[Frame > tp1 & Frame <= tp2,
		      .(dx = diff(Xmm), dy = diff(Ymm), dt = diff(Time_sec), t = Time_sec[-1]), by = 'N_ind']
	spd[, .(v = sqrt(dx**2 + dy**2)/dt, t = t), by = 'N_ind'][v <= quantile(v, p = 0.99)]
}))

spd_det_p3 <- rbindlist(lapply(det, function(i){
	tp2 <- max(rbindlist(i@food)[['t']])
	spd <- i@data[Frame > tp2,
		      .(dx = diff(Xmm), dy = diff(Ymm), dt = diff(Time_sec), t = Time_sec[-1]), by = 'N_ind']
	spd[, .(v = sqrt(dx**2 + dy**2)/dt, t = t), by = 'N_ind'][v <= quantile(v, p = 0.99)]
}))



ggplot(data = spd_det_p1, aes(v/50)) + 
	geom_histogram(color = 'black', fill = 'grey75')+
	scale_x_continuous('Speed (nodes / s)') +
	scale_y_continuous('Frequency')

ggplot(data = spd_det_p2, aes(v/50)) + 
	geom_histogram(color = 'black', fill = 'grey75')+
	scale_x_continuous('Speed (nodes / s)') +
	scale_y_continuous('Frequency')

ggplot(data = spd_det_p3, aes(v/50)) + 
	geom_histogram(color = 'black', fill = 'grey75')+
	scale_x_continuous('Speed (nodes / s)') +
	scale_y_continuous('Frequency')


ggplot(data = rbindlist(list(spd_det_p1, spd_det_p2, spd_det_p3), idcol = TRUE), 
       aes(v/50, fill = factor(.id, labels = c('Exploration', 'Exploitation', 'Relaxation')))) + 
	geom_histogram(color = 'black', alpha = 0.6, aes(y = after_stat(density)),
		       breaks = seq(0, 0.8, 0.01), position = 'identity')+
	scale_x_continuous('Speed (nodes / s)') +
	scale_y_continuous('Density') +
	scale_fill_manual('', values = c('gold3','mediumpurple', 'firebrick2'))






