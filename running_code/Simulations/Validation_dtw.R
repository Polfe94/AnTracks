library(data.table)
library(ggplot2)
library(Rcpp)
path <- '~/research/gits/AutomatAnts/results/food_conditions/det_1000/'
load('~/research/gits/AnTracks/data/det_ref_seqs.RData')


#' C version of `moving_average`
Rcpp::cppFunction('
NumericVector movingAverage(NumericVector x, int t, int overlap = 0) {
    int slide = t / 2;
    
    if (overlap > slide) {
        Rcpp::warning("Overlap is bigger than window size. Capping overlap to window size.");
        overlap = slide;
    }
    
    int x0 = 1 + slide;
    int xn = x.size() - slide;
    int wsize = slide - overlap + 1;
    
    int N = (xn - x0) / wsize + 1;
    
    NumericVector result(N);
    
    for (int i = 0; i < N; i++) {
        int start = x0 + i * wsize - slide ;
        int end = start + 2*slide;
        
        double sum = 0.0;
        int count = 0;
        for (int j = start; j <= end; j++) {
            sum += x[j - 1];
            count += 1;
        }
        result[i] = sum / count;
    }
    
    return result;
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

moving_average <- function(x, t, overlap = 0){
	
	slide <- t %/% 2
	if(overlap > slide){
		warning('Overlap is bigger than window size. Capping overlap to window size.')
		overlap <- slide
	}
	x0 <- 1 + slide
	xn <- length(x) - slide
	
	sq <- seq(x0, xn, (slide - overlap + 1))
	v <- vapply(sq, function(i){
		mean(x[seq(i-slide, i+slide)])
	}, numeric(1))
	
	# fill in start and finish positions
	if(overlap == slide){
		c(x[1:slide], v, x[(xn+1):length(x)])
	} else {
		v
	}
	
}

arrange.dt <- function(path){
	f <- list.files(path, pattern = '.csv')
	results <- vector('list', length(f))
	for(i in seq_along(f)){
		x <- read.csv(paste0(path, f[i]))[, -1]
		setDT(x)
		x[['Frame']] <- x[['Frame']] * 2
		x <- merge(data.table(Frame = seq_len(21600)), x,
			   all.x = TRUE, by = 'Frame')
		cols <- c('N', 'I', 'SiOut')
		x[1, cols] <- 0
		for (ii in cols){
			x[[ii]] <- fillVec(x[[ii]])
		}
		results[[i]] <- x
	}
	rbindlist(results, idcol = TRUE)
}

## ~ 45 sec
t0 <- Sys.time()
dt <- arrange.dt('~/research/gits/AutomatAnts/results/food_conditions/det_1000/')
Sys.time()-t0

# for(i in seq_len(max(dt[['.id']]))){
# 	png(paste0('~/research/gits/AutomatAnts/plots/validation_1/det_', i, '.png'),
# 	    1920, 1080, res = 120)
# 	print(ggplot(data = dt[.id == i], aes(Frame, N)) + geom_line()+
# 		theme_bw()+
# 		scale_y_continuous('Number of Ants', breaks = seq(0,100, 5))+
# 		scale_x_continuous('Time (min)', breaks = seq(0, 21600, 15*120),
# 				   labels = seq(0, 180, 15))+
# 		theme(axis.title = element_text(size = 18),
# 		      axis.text = element_text(size = 16)))
# 	dev.off()
# }


#### ONLY AVERAGE
dt.dtw <- function(dt, ref){
	y <- moving_average(ref, 60)
	l <- max(dt[['.id']])
	result <- vector('list', l)
	for(i in seq_len(l)){
		x <- moving_average(dt[.id == i, N], 60)
		output <- unlist(dtw::dtw(x, y, 
					  distance.only = TRUE,
					  step.pattern = dtw::asymmetric)[c('distance',
					  				  'normalizedDistance')])
		result[[i]] <- data.frame(.id = i, d = output[1], nrm_d = output[2]) 
		gc()
	}
	rbindlist(result)
}

dtw_avg <- dt.dtw(dt[, c('.id', 'N')], det_NS[['avg']])

# 1 = Nice, 0 = Low N, 2 = Plateau, 3 = late-start, 4 = oscilation
k <- c(1, 1, 0, 2, 2, 3, 2, 2, 2, 0, 1, 3, 0, 2, 2, 1, 1, 2, 2, 1, 1, 1, 3,
       1, 1, 1, 3, 0, 1, 2, 2, 2, 0, 0, 1, 3, 1, 1, 1, 0, 1, 1, 3, 1, 3, 2,
       0, 1, 2, 1, 2, 1, 1, 1, 1, 2, 3, 0, 1, 1, 3, 1, 2, 3, 2, 0, 1, 2, 1,
       2, 1, 1, 3, 2, 1, 1, 1, 0, 3, 1, 3, 1, 1, 1, 1, 3, 3, 1, 1, 2, 1, 3, 
       1, 1, 3, 3, 2, 3, 1, 1, 2, 1, 0, 1, 3, 2, 1, 1, 2, 1, 3, 3, 1, 0, 3,
       3, 1, 3, 1, 3, 3, 3, 1, 1, 1, 1, 2, 1, 1, 3, 3, 1, 3, 1, 1, 3, 3, 1,
       3, 3, 1, 3, 3, 1, 2, 3, 1, 2, 2, 1, 3, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1,
       1, 1, 2, 0, 1, 1, 3, 3, 1, 1, 0, 1, 1, 2, 1, 2, 3, 2, 2, 3, 1, 1, 1,
       1, 1, 1, 3, 1, 1, 2, 1, 1, 1, 1, 2, 3, 3, 2, 2, 1, 3, 3, 1, 3, 1, 3,
       1, 1, 0, 1, 3, 1, 1, 1, 2, 3, 3, 0, 1, 2, 2, 1, 1, 3, 0, 3, 2, 2, 1,
       1, 1, 1, 1, 1, 3, 3, 1, 3, 1, 2, 1, 2, 1, 3, 1, 3, 3, 3, 1, 1, 2, 2, 
       1, 4, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 4, 1, 2, 1, 3, 2, 1, 1, 1, 
       1, 0, 0, 3, 1, 1, 3, 1, 3, 3, 0, 3, 1, 1, 2, 3, 1, 1, 1, 3, 1, 2, 1, 
       1, 1, 1, 0, 1, 2, 0, 3, 2, 1, 1, 0, 2, 1, 3, 1, 1, 1, 1, 0, 2, 1, 3,
       1, 3, 1, 1, 1, 3, 2, 0, 1, 3, 1, 1, 1, 1, 1, 1, 0, 0, 2, 3, 2, 1, 0,
       1, 2, 1, 1, 2, 0, 1, 2, 1, 3, 1, 3, 3, 1, 2, 1, 2, 3, 1, 2, 1, 1, 2, 
       2, 1, 1, 3, 1, 0, 0, 2, 3, 1, 2, 2, 1, 1, 1, 2, 3, 2, 2, 2, 1, 1, 1,
       3, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 3, 3, 1, 1, 1, 1, 3, 1, 3, 1, 
       3, 2, 1, 3, 3, 1, 1, 3, 3, 2, 0, 1, 3, 1, 2, 2, 1, 1, 3, 0, 1, 1, 3, 
       1, 3, 3, 3, 3, 2, 1, 0, 2, 2, 1, 2, 2, 1, 4, 1, 3, 3, 2, 1, 2, 2, 1, 
       1, 1, 1, 1, 3, 2, 1, 3, 1, 1, 3, 0, 3, 1, 2, 3, 1, 2, 0, 3, 0, 3, 3,
       3, 0, 1, 3, 2, 3, 3, 0, 3, 1, 0, 1, 1, 3, 2, 3, 2, 3, 2, 2, 1, 0, 3, 
       1, 3, 1, 1, 1, 2, 3, 1, 1, 0, 0, 1, 1, 3, 3, 1, 1, 0, 3, 1, 2, 1, 3,
       3, 3, 2, 1, 2, 3, 2, 1, 1, 1, 2, 4, 1, 3, 2, 1, 1, 2, 3, 2, 1, 2, 2,
       1, 1, 1, 2, 3, 2, 1, 1, 2, 1, 2, 3, 1, 1, 3, 4, 0, 1, 1, 1, 1, 1, 1, 
       3, 0, 1, 1, 0, 0, 0, 3, 3, 1, 1, 3, 2, 1, 4, 1, 2, 3, 1, 2, 0, 1, 3, 
       0, 0, 0, 0, 1, 1, 1, 2, 3, 1, 3, 3, 1, 1, 3, 1, 1, 2, 3, 2, 3, 2, 1,
       3, 2, 1, 1, 1, 2, 1, 2, 0, 1, 0, 3, 1, 2, 1, 1, 1, 1, 2, 3, 2, 0, 3, 
       3, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 0, 3, 3, 1, 4, 1, 1, 3, 1, 1, 3, 1,
       1, 0, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 0, 3, 1, 1, 3, 2, 1, 3, 2, 1, 
       1, 3, 3, 3, 2, 2, 1, 1, 1, 3, 2, 1, 2, 3, 1, 0, 1, 3, 2, 3, 1, 4, 1,
       2, 2, 2, 2, 3, 1, 2, 0, 2, 1, 1, 1, 1, 2, 3, 2, 0, 3, 2, 1, 1, 0, 1,
       3, 1, 1, 3, 3, 1, 1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 3, 1, 0, 1, 2, 3, 1, 
       1, 2, 3, 1, 2, 2, 1, 2, 3, 3, 1, 2, 1, 3, 3, 1, 1, 2, 2, 0, 2, 0, 2,
       1, 2, 1, 2, 1, 3, 1, 2, 1, 1, 1, 3, 1, 2, 2, 1, 2, 3, 0, 0, 0, 1, 3, 
       1, 1, 1, 2, 2, 2, 1, 3, 3, 1, 2, 3, 1, 3, 1, 2, 2, 1, 3, 2, 1, 2, 0, 
       1, 3, 3, 1, 1, 2, 0, 1, 1, 1, 2, 1, 1, 3, 2, 1, 3, 2, 1, 1, 1, 0, 2,
       1, 3, 3, 3, 3, 1, 1, 1, 2, 1, 1, 3, 1, 2, 3, 3, 3, 1, 1, 0, 3, 2, 0,
       1, 1, 1, 1, 0, 3, 1, 3, 1, 0, 3, 0, 2, 0, 1, 1, 1, 3, 1, 1, 2, 0, 1, 
       1, 1, 1, 1, 4, 1, 1, 1, 3, 3, 2, 2, 1, 0, 1, 1, 1, 2, 2, 3, 1, 3, 1, 
       1, 3, 2, 3, 2, 2, 1, 1, 0, 1, 1, 1, 2, 3, 2, 2, 1, 3, 3, 0, 0, 0, 2,
       3, 3, 0, 2, 1, 0, 3, 1, 1, 1, 3, 1, 0, 3, 1, 1, 3, 1, 1, 1, 0, 1, 2, 
       1, 0, 3, 1, 3, 1, 3, 3, 4, 1, 3, 0, 1, 1, 1, 1, 1, 2, 4, 3, 2, 3, 2, 
       1, 3, 1, 0, 1, 2, 1, 1, 1, 3, 1)

dtw_avg[['k']] <- k

png('/home/polfer/research/thesis/capitol_1/preliminary_results/simulation_classification/all_classes.png',
    1920, 1080, res = 120)
ggplot(data = dtw_avg, aes(.id, nrm_d, color = factor(k))) + geom_point(size = 3)+
	scale_color_viridis_d('Class', labels = c('Low N', 'Experiment',
						  'Plateau', 'Late-start', 'Oscillations')) + 
	theme_bw() +
	geom_hline(yintercept = 4.5, linetype = 2, linewidth = 1) +
	scale_y_continuous('Normalized DTW distance', breaks = c(4.5, seq(0, 15, 2.5)))+
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      axis.text.y = element_text(size = 16),
	      axis.title.y = element_text(size = 18),
	      legend.text = element_text(size = 16),
	      legend.title = element_text(size = 18))
dev.off()

M <- matrix(c(sum(k == 1 & dtw_avg[['nrm_d']] < 4.5),
	 sum(k != 1 & dtw_avg[['nrm_d']] >= 4.5),
	 sum(k != 1 & dtw_avg[['nrm_d']] < 4.5),
	 sum(k == 1 & dtw_avg[['nrm_d']] >= 4.5)), ncol = 2, nrow = 2,
       dimnames = list(c('Positive', 'Negative'), c('True', 'False')))

cnf_matrix <- c(Acc = c(sum(M[, 1]) / sum(M)),
  Sens = M[1, 1] / (M[1, 1] + M[2, 2]), 
  Spec = M[2, 1] / (M[2, 1] + M[1, 2]))

png('/home/polfer/research/thesis/capitol_1/preliminary_results/simulation_classification/binary_class_underThreshold.png',
    1920, 1080, res = 120)
ggplot(data = dtw_avg, aes(.id, nrm_d, color = factor(ifelse(k == 1, 0, 1)))) + geom_point(size = 3)+
	scale_color_viridis_d('Class', labels = c('Experiment', 'Others')) + 
	theme_bw() +
	scale_y_continuous('Normalized DTW distance', breaks = seq(0, 5, 1), limits = c(0, 4.5))+
	theme(axis.title.x = element_blank(),
	      axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      plot.title = element_text(size = 20),
	      axis.text.y = element_text(size = 16),
	      axis.title.y = element_text(size = 18),
	      legend.text = element_text(size = 16),
	      legend.title = element_text(size = 18))+
	ggtitle(paste0(names(cnf_matrix),'=', round(cnf_matrix, 3), collapse = '; '))
dev.off()
