exp_condition <- c('det', 'sto', 'nf')
source('G:/research/MoveAnts/code/config.R')
source('C:/Users/POL/Desktop/gants_current/code/dev/dev_functions.R')

hex <- hex90
colnames(hex) <- c('x', 'y')
h <- hex90[hex90$rotY > 980, ]
colnames(h) <- c('x', 'y')
source('G:/research/2022/AnTracks/src/generic.R')
source('G:/research/2022/AnTracks/src/coords.R')

det <- lapply(expL$det, function(i) init(i, refcoords = h,
                                         class = 'coords'))
sto <- lapply(expL$sto, function(i) init(i, refcoords = h,
                                         class = 'coords'))
nf <- lapply(expL$nf, function(i) init(i, refcoords = h,
                                         class = 'coords'))


det_N <- lapply(det, get_N)
sto_N <- lapply(sto, get_N)
nf_N <- lapply(nf, get_N)

for(i in seq_along(det_N)){
        x <- seq(0, max(det[[i]]$data$Frame))
        y <- det_N[[i]]
        write.table(cbind(x, y), 
                    file = paste('C:/Users/POL/Desktop/N_data_allExps/', names(det)[i], '.csv', sep = ''),
                    sep = ',')
}

for(i in seq_along(sto_N)){
        x <- seq(0, max(sto[[i]]$data$Frame))
        y <- sto_N[[i]]
        write.table(cbind(x, y), 
                    file = paste('C:/Users/POL/Desktop/N_data_allExps/', names(sto)[i], '.csv', sep = ''),
                    sep = ',')
}

for(i in seq_along(nf_N)){
        x <- seq(0, max(nf[[i]]$data$Frame))
        y <- nf_N[[i]]
        write.table(cbind(x, y), 
                    file = paste('C:/Users/POL/Desktop/N_data_allExps/', names(nf)[i], '.csv', sep = ''),
                    sep = ',')
}

