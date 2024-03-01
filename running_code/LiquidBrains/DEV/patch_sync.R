path <- "/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.5/"
f <- list.files(path)



# first_patch_detection <- function(food){
# 	ifelse(which.min(food[['t']]) < 7, 1, 2)
# }


revert_node <- function(n){
	n <- do.call('rbind', n)
	apply(n, 1, function(i) paste0('(', paste0(i, collapse = ', '), ')'))
}

patch_detection_diff <- function(path){
	food <- data.table(read_parquet(paste0(path)))[, node := revert_node(node)]
	food <- merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y', 't')]
	food[['p']] <- c(rep(1, 6), rep(2, 6))
	food[, .(mint = min(t)), by = 'p'][, .(diff = abs(diff(mint)))][['diff']]
}

food_collection_diff <- function(path){
	food <- data.table(read_parquet(paste0(path)))[, node := revert_node(node)]
	food <- merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y', 't')]
	food[['p']] <- c(rep(1, 6), rep(2, 6))
	food[order(t), .(diff = abs(diff(t)))][['diff']]
}



path <- "/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.1/"
f <- list.files(path)
pdiff_rho01 <- as.numeric(unlist(lapply(f[grepl('food', f)], function(i){
	patch_detection_diff(paste0(path, i))
})))
cdiff_rho01 <- as.numeric(unlist(lapply(f[grepl('food', f)], function(i){
	food_collection_diff(paste0(path, i))
})))


path <- "/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.5/"
f <- list.files(path)
pdiff_rho05 <- as.numeric(unlist(lapply(f[grepl('food', f)], function(i){
	patch_detection_diff(paste0(path, i))
})))
cdiff_rho05 <- as.numeric(unlist(lapply(f[grepl('food', f)], function(i){
	food_collection_diff(paste0(path, i))
})))


path <- "/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.9/"
f <- list.files(path)
pdiff_rho09 <- as.numeric(unlist(lapply(f[grepl('food', f)], function(i){
	patch_detection_diff(paste0(path, i))
})))
cdiff_rho09 <- as.numeric(unlist(lapply(f[grepl('food', f)], function(i){
	food_collection_diff(paste0(path, i))
})))


boxplot(pdiff_rho01, pdiff_rho05, pdiff_rho09)
# boxplot(log(pdiff_rho01), log(pdiff_rho05), log(pdiff_rho09))
boxplot(log(cdiff_rho01), log(cdiff_rho05), log(cdiff_rho09))


path <- "/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.5/"
f <- list.files(path)
