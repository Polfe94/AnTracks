
KL_collective2 <- suppressMessages({
	rbindlist(list(coll_statictp1 = rbindlist(lapply(unique(coll_static[['rho']]), function(i){
		data.frame(rho = i, kl = KL.empirical(coll_static[rho == i, tp1], 
						      det_tps[['tp1']]), var = 'tp1', lbl = 'static')
	})),
	coll_statictp2 = rbindlist(lapply(unique(coll_static[['rho']]), function(i){
		
		data.frame(rho = i, kl = KL.empirical(coll_static[rho == i, tp2], 
						      det_tps[['tp2']]), var = 'tp2', lbl = 'static')
	})),
	
	
	coll_shutdowntp1 = rbindlist(lapply(unique(coll_shutdown[['rho']]), function(i){
		data.frame(rho = i, kl = KL.empirical(coll_shutdown[rho == i, tp1], 
						      det_tps[['tp1']]), var = 'tp1', lbl = 'shutdown')
	})),
	coll_shutdowntp2 = rbindlist(lapply(unique(coll_shutdown[['rho']]), function(i){
		
		data.frame(rho = i, kl = KL.empirical(coll_shutdown[rho == i, tp2], 
						      det_tps[['tp2']]), var = 'tp2', lbl = 'shutdown')
	}))))})
	


KL_collective <- suppressMessages({
	nbins <- 100
	dtp1 <- range(det_tps[['tp1']])
	dtp2 <- range(det_tps[['tp2']])
	rbindlist(list(coll_static = rbindlist(lapply(unique(coll_static[['rho']]), function(i){
		data <- coll_static[rho == i & is.finite(tp1) & is.finite(tp2) & tp1 > 0 & tp2 > 0]
		trange_x <- range(c(data[, tp1], dtp1))
		trange_y <- range(c(data[, tp2], dtp2))
		sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
		sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
		data[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
		ndata <- nrow(data)
		QX <- det_tps[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
		QX_tp1 <- as.numeric(table(QX[['bin_x']]) / nrow(QX))
		QX_tp2 <- as.numeric(table(QX[['bin_y']]) / nrow(QX))
		PX_tp1 <- as.numeric(table(data[['bin_x']]) / ndata)
		PX_tp2 <- as.numeric(table(data[['bin_y']]) / ndata)
		df_1 <- data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'static')
		df_2 <- data.frame(rho = i, kl = KL(rbind(PX_tp2, QX_tp2)), var = 'tp2', lbl = 'static')
		rbind(df_1, df_2)
	})), coll_shutdown = rbindlist(lapply(unique(coll_shutdown[['rho']]), function(i){
		data <- coll_shutdown[rho == i & is.finite(tp1) & is.finite(tp2) & tp1 > 0 & tp2 > 0]
		trange_x <- range(c(data[, tp1], dtp1))
		trange_y <- range(c(data[, tp2], dtp2))
		sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
		sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
		data[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
		ndata <- nrow(data)
		QX <- det_tps[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
		QX_tp1 <- as.numeric(table(QX[['bin_x']]) / nrow(QX))
		QX_tp2 <- as.numeric(table(QX[['bin_y']]) / nrow(QX))
		PX_tp1 <- as.numeric(table(data[['bin_x']]) / ndata)
		PX_tp2 <- as.numeric(table(data[['bin_y']]) / ndata)
		df_1 <- data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'shutdown')
		df_2 <- data.frame(rho = i, kl = KL(rbind(PX_tp2, QX_tp2)), var = 'tp2', lbl = 'shutdown')
		rbind(df_1, df_2)
	}))))})


KL_collective_log <- suppressMessages({
	nbins <- 100
	log_det_tps <- det_tps
	log_det_tps[['tp1']] = log(log_det_tps[['tp1']])
	log_det_tps[['tp2']] = log(log_det_tps[['tp2']])
	dtp1 <- range(log_det_tps[['tp1']])
	dtp2 <- range(log_det_tps[['tp2']])
	rbindlist(list(coll_static = rbindlist(lapply(unique(coll_static[['rho']]), function(i){
		data <- coll_static[rho == i & is.finite(tp1) & is.finite(tp2) & tp1 > 0 & tp2 > 0]
		trange_x <- range(c(data[, log(tp1)], dtp1))
		trange_y <- range(c(data[, log(tp2)], dtp2))
		sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
		sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
		data <- data[, bin_x := cut(log(tp1), breaks = sq_x)][, bin_y := cut(log(tp2), breaks = sq_y)]
		ndata <- nrow(data)
		QX <- log_det_tps[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
		QX_tp1 <- as.numeric(table(QX[['bin_x']]) / nrow(QX))
		QX_tp2 <- as.numeric(table(QX[['bin_y']]) / nrow(QX))
		PX_tp1 <- as.numeric(table(data[['bin_x']]) / ndata)
		PX_tp2 <- as.numeric(table(data[['bin_y']]) / ndata)
		df_1 <- data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'static')
		df_2 <- data.frame(rho = i, kl = KL(rbind(PX_tp2, QX_tp2)), var = 'tp2', lbl = 'static')
		rbind(df_1, df_2)
	})), coll_shutdown = rbindlist(lapply(unique(coll_shutdown[['rho']]), function(i){
		data <- coll_shutdown[rho == i & is.finite(tp1) & is.finite(tp2) & tp1 > 0 & tp2 > 0]
		trange_x <- range(c(data[, log(tp1)], dtp1))
		trange_y <- range(c(data[, log(tp2)], dtp2))
		sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
		sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
		data <- data[, bin_x := cut(log(tp1), breaks = sq_x)][, bin_y := cut(log(tp2), breaks = sq_y)]
		ndata <- nrow(data)
		QX <- log_det_tps[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
		QX_tp1 <- as.numeric(table(QX[['bin_x']]) / nrow(QX))
		QX_tp2 <- as.numeric(table(QX[['bin_y']]) / nrow(QX))
		PX_tp1 <- as.numeric(table(data[['bin_x']]) / ndata)
		PX_tp2 <- as.numeric(table(data[['bin_y']]) / ndata)
		df_1 <- data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'shutdown')
		df_2 <- data.frame(rho = i, kl = KL(rbind(PX_tp2, QX_tp2)), var = 'tp2', lbl = 'shutdown')
		rbind(df_1, df_2)
	}))))})

KL_ <- suppressMessages({
	nbins <- 100
	dmint <- range(det_global_eff[['mint']])
	dmaxt <- range(det_global_eff[['maxt']])
	rbindlist(list(static_ = rbindlist(lapply(unique(static_[['rho']]), function(i){
		trange_x <- range(c(static_[rho == i & is.finite(mint), mint], dmint))
		trange_x[1] <- trange_x[1] - .Machine$double.eps
		trange_x[2] <- trange_x[2] + .Machine$double.eps
		trange_y <- range(c(static_[rho == i & is.finite(maxt), maxt], dmaxt))
		trange_y[1] <- trange_y[1] - .Machine$double.eps
		trange_y[2] <- trange_y[2] + .Machine$double.eps
		sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
		sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
		data <- static_[rho == i][, bin_x := cut(mint, breaks = sq_x)][, bin_y := cut(maxt, breaks = sq_y)]
		ndata <- nrow(data[is.finite(mint) & is.finite(maxt)])
		QX <- det_global_eff[, bin_x := cut(mint, breaks = sq_x)][, bin_y := cut(maxt, breaks = sq_y)]
		QX_mint <- as.numeric(table(QX[['bin_x']]) / nrow(QX))
		QX_maxt <- as.numeric(table(QX[['bin_y']]) / nrow(QX))
		PX_mint <- as.numeric(table(data[['bin_x']]) / ndata)
		PX_maxt <- as.numeric(table(data[['bin_y']]) / ndata)
		df_1 <- data.frame(rho = i, kl = KL(rbind(PX_mint, QX_mint)), var = 'mint', lbl = 'static')
		df_2 <- data.frame(rho = i, kl = KL(rbind(PX_maxt, QX_maxt)), var = 'maxt', lbl = 'static')
		rbind(df_1, df_2)
	})), shutdown_ = rbindlist(lapply(unique(shutdown_[['rho']]), function(i){
		trange_x <- range(c(shutdown_[rho == i & is.finite(mint), mint], dmint))
		trange_x[1] <- trange_x[1] - .Machine$double.eps
		trange_x[2] <- trange_x[2] + .Machine$double.eps
		trange_y <- range(c(shutdown_[rho == i & is.finite(maxt), maxt], dmaxt))
		trange_y[1] <- trange_y[1] - .Machine$double.eps
		trange_y[2] <- trange_y[2] + .Machine$double.eps
		sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
		sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
		data <- shutdown_[rho == i][, bin_x := cut(mint, breaks = sq_x)][, bin_y := cut(maxt, breaks = sq_y)]
		ndata <- nrow(data[is.finite(mint) & is.finite(maxt)])
		QX <- det_global_eff[, bin_x := cut(mint, breaks = sq_x)][, bin_y := cut(maxt, breaks = sq_y)]
		QX_mint <- as.numeric(table(QX[['bin_x']]) / nrow(QX))
		QX_maxt <- as.numeric(table(QX[['bin_y']]) / nrow(QX))
		PX_mint <- as.numeric(table(data[['bin_x']]) / ndata)
		PX_maxt <- as.numeric(table(data[['bin_y']]) / ndata)
		df_1 <- data.frame(rho = i, kl = KL(rbind(PX_mint, QX_mint)), var = 'mint', lbl = 'shutdown')
		df_2 <- data.frame(rho = i, kl = KL(rbind(PX_maxt, QX_maxt)), var = 'maxt', lbl = 'shutdown')
		rbind(df_1, df_2)
	}))))})
	