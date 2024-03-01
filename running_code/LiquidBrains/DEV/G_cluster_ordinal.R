library(data.table)
library(ggplot2)
library(segmented)
library(latex2exp)


Rcpp::cppFunction('
CharacterVector fillCharVec(CharacterVector x) {
    for (int i = 1; i < x.size(); i++) {
        if (x[i] == NA_STRING) {
            x[i] = x[i - 1];
        }
    }
    return x;
}
')

Rcpp::cppFunction('
NumericVector fillVec(NumericVector x) {
    for (int i = 1; i < x.size(); i++) {
        if (NumericVector::is_na(x[i])) {
            x[i] = x[i - 1];
        }
    }
    return x;
}
')

path_results <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/Jij_0.4_5/'
files <- list.files(path_results)
d <- data.table(read.csv(paste0(path_results, files[1])))[-1, ]

ggplot(data = d[, .N, by = 'Frame'], aes(Frame, N)) + geom_path() 
ggplot(data = reshape2::melt(d[, c('gNest', 'gArena')])) + 
	geom_histogram(aes(x = value, fill = variable), color = 'black', breaks = seq(0, 1, 0.02)) +
	scale_fill_viridis_d(TeX('Sensitivity ($g$)'), labels = c('Nest', 'Arena'))


finite_mean <- function(x){
	x <- mean(x, na.rm = T)
	if(!is.finite(x) | is.na(x)){
		0
	} else {
		x
	}
}
ggplot(data = reshape2::melt(d[, .(gNest = finite_mean(gNest),
				   gArena = finite_mean(gArena)),
				   by = 'Frame'])) + 
	geom_histogram(aes(x = value, fill = variable), color = 'black', breaks = seq(0, 1, 0.02)) +
	scale_fill_viridis_d(TeX('Sensitivity ($g$)'), labels = c('Nest', 'Arena'))




d <- data.table(read.csv(paste0(path_results, files[11])))[-1, ]
d <- data.table(read.csv(paste0(path_results, files[30])))[-1, ]
ggplot(data = d[, .N, by = 'Frame'], aes(Frame, N)) + geom_path() 
ggplot(data = reshape2::melt(d[, .(gNest = finite_mean(gNest),
				   gArena = finite_mean(gArena)),
				   by = 'Frame'])) + 
	geom_histogram(aes(x = value, fill = variable), color = 'black', breaks = seq(0, 1, 0.02)) +
	scale_fill_viridis_d(TeX('Sensitivity ($g$)'), labels = c('Nest', 'Arena'))




d <- data.table(read.csv(paste0(path_results, files[11])))[-1, ]














# colnames(s)[1] <- 'Frame'
# s[['Frame']] <- s[['Frame']]*2
M <- matrix(nrow = 21600, ncol = 100)
M[1, ] <- 'Inactive'
for(i in unique(s[['Frame']])){
	ids <- s[Frame == i, c('id', 'states')]
	v <- unique(ids, by = 'id', fromlast = TRUE)
	M[i, 1+v[['id']]] <- v[['states']]
}




for(i in 1:100){
	M[, i] <- fillCharVec(M[, i])
}

data <- data.table(Frame = 1:21600, t(apply(M, 1, function(i){
	table(factor(i, levels = c('Active', 'Food', 'Inactive', 'Nest'))) / 100
})))
data[['Si']] <- merge(data.table(Frame = 1:21600), 
		      s[, .(Si = mean(si_nest)), by = 'Frame'], all = TRUE,
		      by = 'Frame')[['Si']]
data[['Si']][1] <- s[1, si_nest]
data[['Si']] <- fillVec(data[['Si']])

ggplot(data = reshape2::melt(data, id.vars = 'Frame'), aes(Frame, value, color = variable)) +
	geom_path() + scale_color_viridis_d('') + ylab('Proportion') +
	scale_x_continuous('Time (min)', breaks = seq(0, 180*120, 30*120),
			   labels = seq(0, 180, 30)) +
	theme_bw()

ggplot(data = data, aes(Si, Active)) +
	geom_point() + scale_color_viridis_d('') + ylab('Proportion') +
	theme_bw()

ggplot(data = data, aes(Si/Nest, Active)) +
	geom_point() + scale_color_viridis_d('') + ylab('Proportion') +
	theme_bw()

ggplot(data = data, aes(Si/ Inactive, Active)) +
	geom_point() + scale_color_viridis_d('') + ylab('Proportion') +
	theme_bw()

H <- data.table(apply(data[, 2:5] *100, 2, as.integer))

Hs <- apply(H, 1, infotheo::entropy)

dt <- c()

for(i in seq(0, 21600, 120)){
	dt <- rbind.data.frame(dt, c(Frame = mean(c(i, i+120))/2,
			  si = mean(data[i:(i+120), Si]),
			  ha = infotheo::entropy(H[i:(i+120), Active]),
			  hn = infotheo::entropy(H[i:(i+120), Nest])))	
}
setDT(dt)
colnames(dt) <- c('Frame', 'Si', 'Ha', 'Hn')

library(pracma)

### IDEAS: AVERAGE PEAK DURATION (SI), AVERAGE HEIGHT (SI), AVERAGE LAGG (SI - ACTIVE)
plot(findpeaks(data[['Si']]))
plot(data[, c('Frame', 'Si')], type = 'l')
points(findpeaks(loess(data[['Si']] ~ data[['Frame']], span = 0.05)$fitted, npeaks = 0,
		 minpeakheight = 0.002)[, 2:1])
abline(v = findpeaks(data[['Si']], npeaks = 0,
		     minpeakheight = 0.003)[, 2], lty = 'dashed')

my.seg <- segmented(lm(Active ~ Frame, data = data), seg.Z = ~ Frame, 
	  psi = list(Frame = c(5000, 10000)))

summary(my.seg)
ggplot(data = data) + geom_line(aes(Frame, Active)) + 
	geom_vline(xintercept = my.seg$psi[, 2], linetype = 'dashed')+
	 ylab('Proportion of active ants') +
	scale_x_continuous('Time (min)', breaks = seq(0, 180*120, 30*120),
			   labels = seq(0, 180, 30)) +
	theme_bw()


library(vars)
my.lm <- lm(Active ~ Si, data = data)
plot.ts(my.lm$residuals)
error.lm <- ur.df(my.lm$residuals, lags = 100, type = 'none')
summary(error.lm)

diff.Si <- diff(data$Si)[-1]
diff.Active <- diff(data$Active)[-1]
error.ecm1 <- my.lm$residuals[-1:-2]
diff.Active1 <- diff(data$Active)[-(nrow(data) - 1)]
diff.Si1 <- diff(data$Si)[-(nrow(data) - 1)]

my.lm_error <- lm(diff.Active ~ error.ecm1 + diff.Active1 + diff.Si1)
summary(my.lm_error)

my.lm_error <- lm(diff.Si ~ error.ecm1 + diff.Active1 + diff.Si1)
summary(my.lm_error)
