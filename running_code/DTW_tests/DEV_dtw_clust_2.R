library(data.table)
library(dtw)
library(dtwclust)

standarize_variable <- function(x, var = 'N', tau = 30, frames = 21600){
	y <- merge.data.table(data.table(Frame = seq(1, frames)), x[, c('Frame', ..var)],
			 all = TRUE, by = 'Frame')
	y[is.na(get(var)), var] <- 0
	y[seq(frames), lapply(.SD, mean), .SDcols = var, by = .(Frame = ((Frame) %/% tau) * tau)]
}

pattern_0 <- data.table(read.csv('~/research/gits/AutomatAnts/results/example_3.csv')) ## replace with experimental determinist average
pattern_1 <- data.table(read.csv('~/research/gits/AutomatAnts/results/with_recruitment/food_conditions/det/det_87.csv'))
pattern_2 <- data.table(read.csv('~/research/gits/AutomatAnts/results/with_recruitment/food_conditions/det/det_78.csv'))

# one small test seems to reveal there is no difference between `tau = 30` and `tau = 60` 
p0 <- standarize_variable(pattern_0, tau = 60)
p1 <- standarize_variable(pattern_1, tau = 60)
p2 <- standarize_variable(pattern_2, tau = 60)

test <- data.table(read.csv('~/research/gits/AutomatAnts/results/with_recruitment/food_conditions/det/det_86.csv'))
test <- standarize_variable(test, tau = 60)


## zscore is used by `dtwclust` package to classify temporal series by 'shape' (then distance to centroid)
ds <- c(dtw(zscore(test[['N']]), zscore(p0[['N']]), step.pattern = asymmetric, distance.only = T)[['normalizedDistance']],
	dtw(zscore(test[['N']]), zscore(p1[['N']]), step.pattern = asymmetric, distance.only = T)[['normalizedDistance']],
	dtw(zscore(test[['N']]), zscore(p2[['N']]), step.pattern = asymmetric, distance.only = T)[['normalizedDistance']])
